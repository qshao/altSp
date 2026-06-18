from __future__ import annotations
import json
import time
import urllib.error
import urllib.request
from typing import Iterable
from .models import IsoformProtein

REST_BASE = "https://grch37.rest.ensembl.org"


def _rest_get(path: str):
    req = urllib.request.Request(
        REST_BASE + path, headers={"Content-Type": "application/json"}
    )
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.load(r)


def resolve_rest(name: str) -> IsoformProtein:
    iso = IsoformProtein(name=name, source="rest", status="unresolved")
    try:
        xrefs = _rest_get(
            f"/xrefs/symbol/homo_sapiens/{name}?content-type=application/json"
        )
    except urllib.error.HTTPError:
        return iso
    tid = next(
        (x["id"] for x in xrefs if x.get("type") == "transcript"), None
    )
    if not tid:
        return iso
    iso.transcript_id = tid
    try:
        seq = _rest_get(
            f"/sequence/id/{tid}?type=protein;content-type=application/json"
        )
        iso.protein_seq = seq.get("seq")
        iso.protein_id = seq.get("id")
        iso.status = "protein" if iso.protein_seq else "noncoding"
    except urllib.error.HTTPError:
        iso.status = "noncoding"  # transcript exists but has no protein
    return iso


def resolve_all(
    names: Iterable[str],
    resolver_map: dict,
    use_rest: bool = True,
    sleep: float = 0.0,
) -> dict[str, IsoformProtein]:
    resolved: dict[str, IsoformProtein] = {}
    for n in sorted(names):
        if n in resolver_map:
            resolved[n] = resolver_map[n]
        elif use_rest:
            resolved[n] = resolve_rest(n)
            if sleep:
                time.sleep(sleep)
        else:
            resolved[n] = IsoformProtein(name=n, status="unresolved")
    return resolved
