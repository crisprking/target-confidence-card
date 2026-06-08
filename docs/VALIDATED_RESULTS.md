# Validated results

Every number below comes from running this repository's code on **real, public ClinVar**
(GRCh38, 2026-06 release) — not synthetic data. The synthetic fixtures in the test paths exist
only to prove plumbing and are clearly labelled as such.

- **Engine:** `reliability.py`, rule `tcc-verdict-v0.3-reviewstatus`
- **Source:** ClinVar `variant_summary.txt.gz`, 8,982,609 rows scanned
- **Universe:** **197,272** clean (≥1-star, two-sided-eligible) Pathogenic/Benign **missense**
  labels over **18,514** genes
- **Missense extraction:** HGVS protein-token heuristic (see *Caveats*)

---

## 1. The funnel

| Stage | Count | Gate |
|------|------:|------|
| Genes with curated missense | 18,514 | — |
| **TRUST** | **212** | two-sided · ≥40 clean labels · ≥50% at ≥2 review stars |
| eligible to estimate AUROC | 50 | ≥20 labels in **each** class |
| **MEASURABLE** | **2** | bootstrap 95% CI width ≤ 0.10 → **USH2A, BRCA1** |
| measurable at the ≥2-star gold standard | **0** | — |

Full tier breakdown: **212 TRUST / 3,903 CAUTION / 14,399 REFUSE**.

The two measurable genes are **mediocre, not good** — that is the point:

| Gene | clean labels | Grantham AUROC | 95% CI width |
|------|------:|:---:|:---:|
| USH2A | 450 | 0.766 | ≤ 0.10 |
| BRCA1 | 591 | 0.717 | ≤ 0.10 |

For ~18,512 other genes we cannot even *state* a predictor's per-gene accuracy with a usable
confidence interval, because the trustworthy ground truth to check it against does not exist.

---

## 2. The five supporting findings

**GEM 1 — real in aggregate, invisible per gene.**
Pooled Grantham AUROC over all 197,272 labels = **0.6714** (95% CI [0.6687, 0.6741]) — pinned to
±0.003. Yet per gene the signal is unmeasurable for all but two. The uncomfortable twist: on the
**higher-quality ≥2-star labels it gets *worse*** (0.6714 → **0.6494**). Same predictor, opposite
verdict, purely from where the labels sit.

**GEM 2 — the empirical substitution map (the standalone showpiece).**
A 20×20 matrix of pathogenic-fraction per amino-acid substitution. Spearman **ρ = +0.66** vs
Grantham distance. Cysteine and tryptophan dominate the pathogenic corner (C→W 86%, C→F 81%,
W→C / W→G / W→S ≈ 80%) — structurally load-bearing residues (disulfides, core packing). The
falsifiable disagreement cells: Grantham **under-warns** W→L (distance only 61, yet 78% pathogenic)
and F↔I (distance 22, ~50–58% pathogenic); it **over-warns** S→L, R→C, Y→C. "Chemistry says mild,
curation says severe" — checkable biology.

**GEM 3 — the wall is part sample size, part predictor skill.**
CI width ≈ **1.31 / √(labels per class)** at Grantham's accuracy; hitting CI ≤ 0.10 implies
**~171 labels per class**. Only **1** gene has that many — yet **2** are measurable, because USH2A
clears at n < 171 where its variance happens to be favorable. So a stronger predictor *might* move
the wall by reaching a higher-AUROC/lower-variance regime — that is the ESM-C question, set up but
not yet answered.

**GEM 4 — benign is the majority *and* the binding constraint, because it's hoarded.**
Global clean labels: **59,234 P/LP : 138,038 B/LB** — benign are the **70%** majority (2.3×). Yet
they are **concentrated**: **Gini 0.70** across 18,514 genes; the top **1,244** genes hold 50% and
the top **5,177** hold 80% of all benign labels; only **15,014** genes have *any* benign label
(≈3,500 have zero). Per gene, missing-benign is still the more common eligibility blocker (**98** vs
**64**), and the median benign-blocked near-TRUST gene carries **67 P/LP but only 8 B/LB**. It is a
*distribution* problem (saturation-screened genes hoard benign while the typical clinical gene is
pathogenic-rich and benign-poor), not aggregate scarcity.

> Methods note: an early pass concluded "benign is under-submitted." The data says the opposite —
> it is concentrated. The corrected, fully data-driven finding is the one above.

**GEM 5 — the quality of the "truth."**
Review stars across all clean labels: **1★ 134,459 (68.2%) · 2★ 58,091 (29.4%) · 3★ 4,714 (2.4%) ·
4★ 8 (0.0%)**. Only **32%** reach ≥2 stars, and expert-panel review is vanishingly rare. The
`REVIEW_BAR = 50%` gate is the TRUST-vs-CAUTION divider among well-labelled genes — which is why the
strict gate yields **212** TRUST genes rather than the ~844 a count-only gate would pass.

Figures for these are produced by `tcc_article_figures.py` → `tcc_funnel.png`,
`tcc_substitution_matrix.png`, `tcc_benign_concentration.png`.

---

## Status: ESM-C sequel — apparatus verified, result NOT yet produced

`tcc_esm_swap.py` is the falsifiable follow-up: replace the funnel's `score()` with ESM-C and ask
whether a real model moves the measurability wall.

What is **validated** (standin, on real ClinVar):

- The Grantham control arm **reproduces the published two genes** ({USH2A, BRCA1}) with matching
  per-gene AUROCs (USH2A 0.767, BRCA1 0.720; global 0.681) — so the sequel's baseline *is* the Two
  Genes baseline, and any gene ESM-C adds is cleanly attributable.
- The read-off-true standin (no label signal) correctly sits **at chance** (median per-gene AUROC
  **0.496**), the resumable cache is per-scorer (a standin run cannot contaminate a later GPU run),
  and the position/offset/vocab read-off and `seq[pos-1]==ref` gate all execute.

What is **NOT** validated:

- **There are no real ESM-C numbers.** Running the model requires a GPU + the `esm` package; the
  module otherwise falls back to the no-signal standin above. The standin's table is plumbing proof,
  **not a result** — do not cite its chance-level figures.
- An older, separate **ESM-2** pass on a different base hinted the honest per-gene median is ~0.87,
  which sets a rough expectation but is not this pipeline's validated output.

**Open prediction (falsifiable, not yet checked):** among eligible genes, TP53 (Grantham AUROC 0.668,
CI width ~0.11 — one hair over the line) is the marquee candidate to become the *third* measurable
gene if ESM-C lifts its ranking toward ~0.85. Benign-starved genes (LDLR, NF1, GAA, MYH7) and
high-AUROC-but-tiny-n genes (FGFR3, RET, MYOC) should **not** move no matter how good the model — they
are gated by class imbalance and sample size, not predictor skill.
