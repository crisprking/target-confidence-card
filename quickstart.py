#!/usr/bin/env python3
"""
quickstart.py - see the Target Confidence Card work in ~5 seconds. No network, no GPU.

It demonstrates the two INDEPENDENT axes the tool keeps separate:

  1) VERDICT (predictor-INDEPENDENT): is a gene's ClinVar ground truth trustworthy
     enough to benchmark ANY predictor against?  ->  TRUST / CAUTION / REFUSE
  2) MEASURABILITY (predictor-DEPENDENT): given that truth, can we pin a predictor's
     per-gene AUROC tightly enough to be useful?  The SAME predictor is "measurable"
     on a large gene and not on a small one - the central finding, in miniature.

Run:
    python quickstart.py

Then try the CLI:
    python reliability.py BRCA1                          # one gene, built-in demo index
    python reliability.py --index genes_index.json TP53  # after running build_index.py
"""
import random

import reliability as tcc


def show_verdicts() -> None:
    print("=" * 72)
    print("1) VERDICTS  (predictor-independent; from the built-in demo index)")
    print("=" * 72)
    for gene in ["BRCA1", "PCSK9", "DEMO_ONESIDED", "DEMO_THIN", "NOTAGENE"]:
        v = tcc.assess_gene(gene)          # no index passed -> built-in DEMO_INDEX
        print(f"\n{gene:14} -> {v.tier.value}")
        print(f"   {v.headline}")
        print(f"   {v.reasons[0]}")
    print("\nTRUST = benchmarkable | CAUTION = usable with caveats |")
    print("REFUSE = not enough trustworthy truth (this includes 'absent from the index').")


def _labeled_scores(n_per_class: int, sep: float, seed: int):
    """A toy predictor: benign scores ~N(0,1), pathogenic ~N(sep,1). Higher sep = better predictor."""
    rng = random.Random(seed)
    scores = [rng.gauss(0.0, 1.0) for _ in range(n_per_class)] + \
             [rng.gauss(sep, 1.0) for _ in range(n_per_class)]
    labels = [0] * n_per_class + [1] * n_per_class
    return scores, labels


def show_measurability() -> None:
    print("\n" + "=" * 72)
    print("2) MEASURABILITY  (predictor-dependent; ONE predictor, three gene sizes)")
    print("=" * 72)
    sep = 0.8                               # a single fixed predictor skill (~AUROC 0.71)
    for n in (300, 60, 12):                 # clean labels available per class
        scores, labels = _labeled_scores(n, sep, seed=7)
        r = tcc.classify_measurability(scores, labels, seed=7)
        auc = "n/a" if r["auroc"] is None else f"{r['auroc']:.3f}"
        ciw = "n/a" if r["ci_width"] is None else f"{r['ci_width']:.3f}"
        print(f"\n{2 * n:4} clean labels ({n}/class):  AUROC {auc}   CI width {ciw}")
        print(f"   -> {r['klass'].upper()}: {r['why']}")
    print("\nSame predictor every time. Whether its accuracy is *measurable* on a gene is")
    print("gated by how much trustworthy truth that gene has - which is exactly why, on")
    print("real ClinVar, only 2 of 18,514 genes clear the bar. And 'measurable' is not")
    print("'good': the two that clear sit at AUROC 0.77 and 0.72. See docs/VALIDATED_RESULTS.md.")


if __name__ == "__main__":
    tcc.grantham_selfcheck()                # the engine's built-in regression check
    show_verdicts()
    show_measurability()
    print("\nNext: `python reliability.py BRCA1`, or build a real index with build_index.py.")
