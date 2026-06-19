"""Resolve before/after protein sequences for the UniProt PCa splicing set and
classify how each annotated isoform changes the canonical protein.

For every reviewed protein in uniprot_prad_splicing/proteins.json we pair the
canonical sequence ("before") with each annotated splice isoform ("after"),
using the actual isoform sequences in isoforms.fasta. For each pair we compute
the divergence point (common prefix / suffix), the length change, and a change
class, and we attach the curated VAR_SEQ descriptions, function and disease.

This is the UniProt-side analogue of analysis/sequence_impact.py (which worked on
the TCGA splice-out/splice-in peptides). Here "before" = canonical isoform,
"after" = an alternatively-spliced isoform — both curated by UniProt.

Output: analysis/uniprot_sequence_impact.json
"""
from __future__ import annotations
import json
import re

JSON_IN = "uniprot_prad_splicing/proteins.json"
FASTA = "uniprot_prad_splicing/isoforms.fasta"
OUT = "analysis/uniprot_sequence_impact.json"


def load_fasta(path: str) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name = None
    buf: list[str] = []
    for line in open(path):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            m = re.match(r">(?:sp|tr)\|([^|]+)\|", line)
            name = m.group(1) if m else line[1:].split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def common_prefix(a: str, b: str) -> int:
    n = min(len(a), len(b))
    i = 0
    while i < n and a[i] == b[i]:
        i += 1
    return i


def common_suffix(a: str, b: str, used: int) -> int:
    n = min(len(a), len(b)) - used
    i = 0
    while i < n and a[-1 - i] == b[-1 - i]:
        i += 1
    return i


def classify(before: str, after: str, pre: int, suf: int) -> tuple[str, str]:
    """Return (class, human description) for a canonical->isoform change."""
    db = len(before) - len(after)
    inner_b = len(before) - pre - suf
    inner_a = len(after) - pre - suf
    if after == before:
        return "identical", "Isoform sequence identical to canonical."
    if pre == len(after) and len(after) < len(before):
        return ("C-terminal truncation",
                f"Isoform stops early: shares the first {pre} residues with the "
                f"canonical protein, losing the C-terminal {len(before) - pre} aa.")
    if suf == len(after) and len(after) < len(before):
        return ("N-terminal truncation",
                f"Isoform lacks the N-terminal {len(before) - len(after)} aa; "
                f"shares the C-terminal {suf} residues (alternative start).")
    if inner_a == 0 and inner_b > 0:
        return ("internal deletion (in-frame)",
                f"Clean internal deletion of {inner_b} aa "
                f"(residues {pre + 1}-{pre + inner_b}); rest of the protein "
                "unchanged.")
    if inner_b == 0 and inner_a > 0:
        return ("internal insertion (in-frame)",
                f"In-frame insertion of {inner_a} aa after residue {pre}.")
    # both inner segments non-empty -> substituted stretch / frameshift-like
    return ("substituted segment",
            f"Residues {pre + 1}..(canonical {pre + inner_b}) are replaced: "
            f"{inner_b} aa swapped for {inner_a} aa before a shared "
            f"C-terminal tail of {suf} aa "
            f"({'frameshift-like / alt reading' if suf < 5 else 'cassette swap'}).")


def gene_of(rec: dict) -> str:
    g = rec.get("genes", [])
    return (g[0].get("geneName") or {}).get("value", "") if g else ""


def function_of(rec: dict) -> str:
    for c in rec.get("comments", []):
        if c.get("commentType") == "FUNCTION":
            t = c.get("texts", [])
            if t:
                return t[0].get("value", "")
    return ""


def diseases_of(rec: dict) -> list[dict]:
    out = []
    for c in rec.get("comments", []):
        if c.get("commentType") == "DISEASE":
            d = c.get("disease", {})
            note = ""
            dd = c.get("texts") or []
            if dd:
                note = dd[0].get("value", "")
            out.append({"id": d.get("diseaseId", ""),
                        "acc": d.get("diseaseAccession", ""),
                        "desc": (d.get("description") or note)})
    return out


def isoform_catalogue(rec: dict) -> list[dict]:
    """List isoforms with their name/synonyms and the VAR_SEQ ids they use."""
    out = []
    for c in rec.get("comments", []):
        if c.get("commentType") == "ALTERNATIVE PRODUCTS":
            for iso in c.get("isoforms", []):
                syn = [s.get("value", "") for s in iso.get("synonyms", [])]
                out.append({
                    "ids": iso.get("isoformIds", []),
                    "name": (iso.get("name") or {}).get("value", ""),
                    "synonyms": syn,
                    "status": iso.get("isoformSequenceStatus", ""),
                    "seqIds": iso.get("sequenceIds", []),
                })
    return out


def varseq_text(rec: dict) -> dict[str, str]:
    """featureId -> short description of each VAR_SEQ change."""
    out = {}
    for f in rec.get("features", []):
        if f.get("type") not in ("Alternative sequence", "VAR_SEQ"):
            continue
        loc = f.get("location", {})
        s = (loc.get("start") or {}).get("value")
        e = (loc.get("end") or {}).get("value")
        alt = f.get("alternativeSequence") or {}
        orig = alt.get("originalSequence")
        alts = alt.get("alternativeSequences")
        if not alt:
            change = f"residues {s}-{e} deleted"
        elif alts:
            change = (f"residues {s}-{e} "
                      f"({(orig or '')[:12]}…) -> {','.join(alts)[:18]}…")
        else:
            change = f"residues {s}-{e} altered"
        out[f.get("featureId", "")] = f"{f.get('description','')}: {change}"
    return out


def main() -> None:
    recs = json.load(open(JSON_IN))
    seqs = load_fasta(FASTA)
    results = []
    for rec in recs:
        acc = rec.get("primaryAccession", "")
        canonical = seqs.get(acc) or (rec.get("sequence") or {}).get("value", "")
        if not canonical:
            continue
        gene = gene_of(rec)
        fn = function_of(rec)
        dis = diseases_of(rec)
        is_pca = any("prostate" in (d["id"] + d["desc"]).lower() for d in dis)
        vseq = varseq_text(rec)
        cat = {tuple(i["ids"]): i for i in isoform_catalogue(rec)}
        iso_records = []
        # every isoform sequence present in the fasta for this accession
        for key in sorted(seqs):
            if not key.startswith(acc + "-"):
                continue
            after = seqs[key]
            pre = common_prefix(canonical, after)
            suf = common_suffix(canonical, after, pre)
            cls, desc = classify(canonical, after, pre, suf)
            meta = cat.get((key,), {})
            vs = [vseq.get(sid, sid) for sid in meta.get("seqIds", [])]
            iso_records.append({
                "isoform_id": key,
                "isoform_name": meta.get("name", ""),
                "synonyms": meta.get("synonyms", []),
                "status": meta.get("status", ""),
                "after_len": len(after),
                "after_seq": after,
                "identical_prefix": pre,
                "identical_suffix": suf,
                "len_change": len(after) - len(canonical),
                "change_class": cls,
                "change_desc": desc,
                "varseq": vs,
            })
        if not iso_records:
            continue
        results.append({
            "accession": acc,
            "gene": gene,
            "protein_name": ((rec.get("proteinDescription", {})
                              .get("recommendedName", {})
                              .get("fullName", {}) or {}).get("value", "")),
            "canonical_len": len(canonical),
            "canonical_seq": canonical,
            "function": fn,
            "diseases": dis,
            "is_prostate_disease": is_pca,
            "n_isoforms": len(iso_records) + 1,
            "isoforms": iso_records,
        })

    json.dump(results, open(OUT, "w"), indent=1)
    n_iso = sum(len(r["isoforms"]) for r in results)
    print(f"proteins with >=1 resolved isoform: {len(results)}")
    print(f"canonical->isoform pairs:           {n_iso}")
    print(f"  with curated prostate disease:    "
          f"{sum(1 for r in results if r['is_prostate_disease'])}")
    cls = {}
    for r in results:
        for i in r["isoforms"]:
            cls[i["change_class"]] = cls.get(i["change_class"], 0) + 1
    for k, v in sorted(cls.items(), key=lambda x: -x[1]):
        print(f"    {k:32s} {v}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
