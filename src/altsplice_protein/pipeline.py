from __future__ import annotations
import os
from .filtering import iter_significant_events, collect_unique_isoforms
from .ensembl_data import download, build_resolver_map, GTF_URL, PEP_URL
from .resolver import resolve_all
from .comparison import compare_event, representative


def write_fasta(path, resolved) -> None:
    with open(path, "w") as f:
        for name in sorted(resolved):
            iso = resolved[name]
            if not iso.protein_seq:
                continue
            f.write(
                f">{iso.name}|{iso.transcript_id}|{iso.protein_id}"
                f"|{iso.gene_symbol}\n"
            )
            s = iso.protein_seq
            for i in range(0, len(s), 60):
                f.write(s[i:i + 60] + "\n")


def _ids(names, resolved) -> str:
    if not names:
        return "."
    return ";".join(
        (resolved[n].protein_id or ".") if n in resolved else "."
        for n in names
    )


def write_events_tsv(path, events, resolved) -> None:
    cols = [
        "Splice_Event", "Gene_Symbol", "Splice_Type", "PSI_Difference",
        "FDR_Difference", "SpliceIn_isoforms", "SpliceIn_protein_ids",
        "SpliceOut_isoforms", "SpliceOut_protein_ids", "has_complete_pair",
        "before_len", "after_len", "comparison",
    ]
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for e in events:
            after_seqs = [
                resolved[n].protein_seq for n in e.splice_in if n in resolved
            ]
            before_seqs = [
                resolved[n].protein_seq for n in e.splice_out if n in resolved
            ]
            has_both = bool(e.splice_in) and bool(e.splice_out)
            comparison = compare_event(before_seqs, after_seqs, has_both)
            brep = representative(before_seqs)
            arep = representative(after_seqs)
            f.write("\t".join([
                e.splice_event, e.gene_symbol, e.splice_type,
                f"{e.psi_difference:g}", f"{e.fdr_difference:g}",
                ",".join(e.splice_in) or ".", _ids(e.splice_in, resolved),
                ",".join(e.splice_out) or ".", _ids(e.splice_out, resolved),
                "yes" if has_both else "no",
                str(len(brep)) if brep else ".",
                str(len(arep)) if arep else ".",
                comparison,
            ]) + "\n")


def write_unresolved(path, resolved) -> None:
    with open(path, "w") as f:
        f.write("isoform_name\ttranscript_id\tstatus\n")
        for name in sorted(resolved):
            iso = resolved[name]
            if iso.status != "protein":
                f.write(f"{name}\t{iso.transcript_id or '.'}\t{iso.status}\n")


def _stats(events, names, resolved) -> dict:
    return {
        "events": len(events),
        "unique_isoforms": len(names),
        "protein": sum(1 for i in resolved.values() if i.status == "protein"),
        "noncoding": sum(1 for i in resolved.values() if i.status == "noncoding"),
        "unresolved": sum(1 for i in resolved.values() if i.status == "unresolved"),
    }


def run_with_files(
    csv_path, gtf_path, pep_path, results_dir,
    fdr_max: float = 0.05, dpsi_min: float = 0.1, use_rest: bool = True,
) -> dict:
    os.makedirs(results_dir, exist_ok=True)
    resolver_map = build_resolver_map(gtf_path, pep_path)
    events = list(iter_significant_events(csv_path, fdr_max, dpsi_min))
    names = collect_unique_isoforms(events)
    resolved = resolve_all(names, resolver_map, use_rest=use_rest)
    write_fasta(os.path.join(results_dir, "proteins.fasta"), resolved)
    write_events_tsv(
        os.path.join(results_dir, "events_proteins.tsv"), events, resolved
    )
    write_unresolved(
        os.path.join(results_dir, "unresolved.txt"), resolved
    )
    return _stats(events, names, resolved)


def run(
    csv_path, data_dir: str = "data", results_dir: str = "results",
    fdr_max: float = 0.05, dpsi_min: float = 0.1, use_rest: bool = True,
) -> dict:
    os.makedirs(data_dir, exist_ok=True)
    gtf = download(GTF_URL, os.path.join(data_dir, "GRCh37.75.gtf.gz"))
    pep = download(PEP_URL, os.path.join(data_dir, "GRCh37.75.pep.all.fa.gz"))
    return run_with_files(
        csv_path, gtf, pep, results_dir, fdr_max, dpsi_min, use_rest
    )
