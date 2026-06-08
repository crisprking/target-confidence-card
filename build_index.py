"""
build_index.py — turn ClinVar's public bulk file into a per-gene index.

This is how the card goes from a 10-gene demo to every gene in ClinVar (~20k),
fully offline and reproducible. Run it once in YOUR environment (open network),
commit the snapshot date, and the app serves real counts for any gene.

WHY BULK, NOT THE LIVE API: parsing variant_summary.txt is far more robust and
reproducible than fighting NCBI eutils filter syntax, and a dated snapshot is
exactly what content-addressing wants.

DATA (commercially clean — ClinVar is public, NIH/NLM):
  https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz

Optional gnomAD constraint (LOEUF), also open:
  gnomAD v4 constraint metrics TSV with columns including a gene symbol and a
  LOEUF column (lof.oe_ci.upper). Pass with --constraint to enrich the flag.

USAGE:
  python build_index.py --clinvar variant_summary.txt.gz --out genes_index.json \
                        [--constraint gnomad_constraint.tsv] [--snapshot 2026-05]

CAVEAT (read this): variant_summary has no clean "missense" column, so missense
is inferred from the protein-change token in the Name field. It's a v1 heuristic
— spot-check a few genes against the ClinVar website before you publish a number.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import sys
from collections import defaultdict

# A missense protein change looks like (p.Arg175His). Exclude nonsense (Ter),
# synonymous (=), frameshift (fs), and indels (del/ins/dup/ext).
_PROT = re.compile(r"\(p\.[A-Za-z]{3}\d+[A-Za-z]{3}\)")
_NONMISSENSE = ("Ter", "fs", "del", "ins", "dup", "ext", "=")


def _is_missense_snv(vtype: str, name: str) -> bool:
    if vtype != "single nucleotide variant":
        return False
    m = _PROT.search(name)
    if not m:
        return False
    token = m.group(0)
    return not any(bad in token for bad in _NONMISSENSE)


def _stars(review_status: str) -> int:
    """ClinVar ReviewStatus -> star rating (0-4). >=2 stars is the review gate's 'solid'."""
    s = (review_status or "").strip().lower()
    if s == "practice guideline":
        return 4
    if s == "reviewed by expert panel":
        return 3
    if s == "criteria provided, multiple submitters, no conflicts":
        return 2
    if s.startswith("criteria provided"):   # single submitter, or conflicting interpretations
        return 1
    return 0


def _bucket(clinsig: str) -> str | None:
    s = clinsig.lower()
    if "conflicting" in s:
        return "conflicting"
    # Exclude pure drug-response / risk-factor / association rows from P-B buckets.
    if any(x in s for x in ("drug response", "risk factor", "association", "protective", "affects")):
        return None
    if "pathogenic" in s:  # covers Pathogenic, Likely pathogenic, Pathogenic/Likely pathogenic
        return "plp"
    if "benign" in s:  # covers Benign, Likely benign, Benign/Likely benign
        return "blb"
    if "uncertain significance" in s:
        return "vus"
    return None


def parse_variant_summary(path: str, assembly: str = "GRCh38", progress_every: int = 1_000_000) -> dict[str, dict]:
    """Stream the bulk file; return {GENE: {plp, blb, vus, conflicting}} for missense SNVs."""
    counts: dict[str, dict[str, int]] = defaultdict(
        lambda: {"plp": 0, "blb": 0, "vus": 0, "conflicting": 0, "ge2star": 0})
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        # ClinVar's header begins with '#AlleleID'; normalise that first key.
        if reader.fieldnames and reader.fieldnames[0].startswith("#"):
            reader.fieldnames[0] = reader.fieldnames[0].lstrip("#")
        n = 0
        for row in reader:
            n += 1
            if progress_every and n % progress_every == 0:
                print(f"  ...{n:,} rows scanned", flush=True)
            if row.get("Assembly") != assembly:
                continue  # avoid double-counting GRCh37 + GRCh38 rows
            gene = (row.get("GeneSymbol") or "").strip()
            if not gene or gene == "-" or ";" in gene:  # skip intergenic / multi-gene
                continue
            if not _is_missense_snv((row.get("Type") or "").strip(), row.get("Name") or ""):
                continue
            bucket = _bucket((row.get("ClinicalSignificance") or "").strip())
            if bucket:
                counts[gene][bucket] += 1
                # ge2star = clean (P/LP or B/LB) labels at >= 2 review stars
                if bucket in ("plp", "blb") and _stars(row.get("ReviewStatus") or "") >= 2:
                    counts[gene]["ge2star"] += 1
    return {g: dict(c) for g, c in counts.items()}


def merge_constraint(index: dict[str, dict], path: str) -> int:
    """Best-effort merge of gnomAD LOEUF. Returns number of genes enriched."""
    enriched = 0
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []
        gene_col = next((c for c in fields if c.lower() in ("gene", "gene_symbol", "symbol")), None)
        loeuf_col = next((c for c in fields if "oe_ci.upper" in c.lower() or c.lower() == "loeuf"), None)
        if not gene_col or not loeuf_col:
            print(f"  ! constraint file missing gene/LOEUF columns ({fields[:6]}…); skipping", file=sys.stderr)
            return 0
        for row in reader:
            g = (row.get(gene_col) or "").strip()
            val = (row.get(loeuf_col) or "").strip()
            if g in index and val not in ("", "NA"):
                try:
                    index[g]["loeuf"] = float(val)
                    enriched += 1
                except ValueError:
                    pass
    return enriched

def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Build a per-gene Target Confidence index from ClinVar's public bulk file.")
    ap.add_argument("--clinvar", required=True,
                    help="variant_summary.txt.gz (ClinVar tab-delimited bulk file)")
    ap.add_argument("--out", default="genes_index.json", help="output index path")
    ap.add_argument("--constraint", default=None,
                    help="optional gnomAD constraint TSV (adds LOEUF context)")
    ap.add_argument("--snapshot", default="unknown",
                    help="ClinVar snapshot label, e.g. 2026-06 (recorded for content-addressing)")
    ap.add_argument("--assembly", default="GRCh38",
                    help="genome assembly to keep (avoids double-counting GRCh37+GRCh38)")
    args = ap.parse_args(argv)

    print(f"parsing {args.clinvar}  (assembly {args.assembly}) ...", flush=True)
    index = parse_variant_summary(args.clinvar, assembly=args.assembly)
    if args.constraint:
        print(f"  enriched {merge_constraint(index, args.constraint)} genes with LOEUF")

    payload = {"snapshot": args.snapshot, "rule": "build_index-v1.1-reviewstars",
               "n_genes": len(index), "genes": index}
    with open(args.out, "w") as fh:
        json.dump(payload, fh)

    trust_gate = sum(1 for r in index.values()
                     if r["plp"] >= 1 and r["blb"] >= 1
                     and (r["plp"] + r["blb"]) >= 40
                     and r.get("ge2star", 0) / max(1, r["plp"] + r["blb"]) >= 0.50)
    print(f"wrote {args.out}: {len(index):,} genes  (~{trust_gate:,} meet the TRUST gate)")
    print(f"query it:  python reliability.py --index {args.out} BRCA1")


if __name__ == "__main__":
    main()
