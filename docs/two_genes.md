# Two Genes

### Of every gene in ClinVar, only two hold enough trustworthy truth to tell you whether a missense predictor actually works there. Here's what that taught me about the data we benchmark on.

*Figures referenced below live in `../figures/`. For Substack, drop the three PNGs inline where each one is marked.*

We are not short on missense predictors. AlphaMissense put a number on ~71 million variants; protein language models read sequence directly; EVE, REVEL, and a dozen others crowd the leaderboards. The standard move is to benchmark them against ClinVar and report an AUROC. I wanted to ask a narrower question first — one that comes *before* "how good is the predictor":

> For a given gene, does ClinVar even hold enough trustworthy, two-sided ground truth to know whether *any* predictor can be benchmarked there at all?

So I built a small, deliberately stubborn tool to answer it. It is refusal-first: it would rather say "I can't tell" than hand you a confident number it can't stand behind. It never rescores a predictor — it only judges the *truth* a predictor would be scored against — so its verdict doesn't move when you swap models. And it is content-addressed: every verdict is hashed over the evidence plus the ClinVar release, so a number I cite today is reproducible and any change to the underlying counts moves the hash.

I pointed it at the whole thing: the public ClinVar release (GRCh38, 2026-06), 8.98 million rows, every gene with curated missense variants — **18,514** of them.

The answer was **two**.

## Two axes that should never be conflated

The tool keeps two questions strictly apart, because collapsing them is how benchmarks lie to you.

The first is the **verdict**, and it does *not* depend on the predictor: is this gene's ClinVar ground truth trustworthy enough to benchmark against? That means two-sided (you need both pathogenic and benign labels, or AUROC is undefined), deep enough (≥40 clean labels), and well-reviewed (at least half of those labels at ClinVar's ≥2-star "multiple submitters, no conflicts" tier). Pass all three and the gene earns **TRUST**; fail softly and it's **CAUTION**; fail hard and it's **REFUSE** — including *refuse by absence*, because a gene with no curated truth gets no opinion.

The second is **measurability**, and it *does* depend on the predictor: given trustworthy truth, can we actually pin the predictor's per-gene AUROC tightly enough to be useful? My bar is a bootstrap 95% confidence interval no wider than 0.10. A gene can be TRUST and still not measurable, if it simply doesn't have enough labels in each class to estimate anything stable.

## The funnel

Here is the collapse, stage by stage:

*(figure: `../figures/tcc_funnel.png`)*

- **18,514** genes with curated missense
- **212** earn TRUST — two-sided, ≥40 clean labels, ≥50% at ≥2 stars
- **50** are even *eligible* to estimate an AUROC (≥20 labels in **each** class)
- **2** are **measurable** — confidence interval ≤ 0.10: **USH2A** and **BRCA1**
- **0** survive if you raise the bar to ≥2-star "gold standard" labels only

Two genes out of eighteen thousand. And the most important thing about those two is that they are **not impressive**. USH2A lands at AUROC 0.766; BRCA1 at 0.717. These are mediocre numbers. That is the entire point: for ~18,512 genes we cannot even *state* how good a predictor is with a straight face, and the two places we can, the honest answer is "so-so." The headline isn't "predictors are good on two genes." It's "for almost every gene, the ground truth to check them simply isn't there."

## Then I went digging

The funnel is the thesis. Underneath it, five findings I didn't expect.

**1. The signal is real in aggregate and invisible per gene.** Pool all 197,272 clean labels and a dumb baseline — Grantham's 1974 amino-acid distance — scores AUROC **0.6714**, with a 95% CI of [0.6687, 0.6741]. Pinned to ±0.003 across 197k variants: unambiguously real. Yet per gene that same signal is unmeasurable for all but two. Same predictor, opposite verdict, purely from where the labels sit. The twist nobody scripted: on the *higher-quality* ≥2-star labels, Grantham gets **worse** (0.6714 → **0.6494**). The cleaner the truth, the less the cheap heuristic explains — a small, slightly uncomfortable result that I think is true and worth sitting with.

**2. "Measurable" is not "good."** Worth repeating as its own point, because it's the one most likely to be misread from a leaderboard. The two genes that clear the bar are the two we can *characterize*, not the two that are *well predicted*. Measurability is a statement about the data's resolving power, not the model's quality.

**3. Benign labels are the majority — and the binding constraint — because they're hoarded.** This is the structural finding, and my first pass got it backwards. Globally the clean labels run **59,234 pathogenic : 138,038 benign** — benign are the **70%** majority, 2.3× the pathogenic count. So benign can't be the bottleneck, right? But per gene, "too few benign" is still the *more common* reason a near-TRUST gene fails eligibility (98 genes vs 64 blocked by too-few pathogenic). Both facts are true only if benign is **concentrated** — and it is, hard. The Gini coefficient across genes is **0.70**; the top ~7% of genes hold half of all benign labels; roughly **3,500** genes have *zero* benign labels at all. The median gene that just misses TRUST on the benign side carries **67 pathogenic but only 8 benign**.

*(figure: `../figures/tcc_benign_concentration.png`)*

A handful of saturation-screened genes (BRCA1 and friends) hoard the benign labels, while the typical clinically-defined gene is pathogenic-rich and benign-starved. It's a *distribution* problem, not a scarcity problem — and that distinction is the difference between "submit more benign variants" and "the few genes we can measure are the few that got systematically screened."

**4. Chemistry says mild; curation says severe — and you can check where.** Instead of trusting Grantham, I built the empirical version: a 20×20 matrix of the observed pathogenic fraction for every amino-acid substitution in the curated set.

*(figure: `../figures/tcc_substitution_matrix.png`)*

It correlates with Grantham at Spearman **ρ = +0.66** — close enough to validate the old chemistry, loose enough to be interesting. Cysteine and tryptophan dominate the pathogenic corner (C→W 86%, C→F 81%, W→C / W→G / W→S all ≈ 80%) — exactly the structurally load-bearing residues you'd predict from disulfides and core packing. The gold is in the **disagreement cells**, the falsifiable biology a reader can go check: Grantham *under-warns* W→L (its distance is only 61, yet 78% of those variants are pathogenic) and the F↔I swap (distance 22, but ~50–58% pathogenic); it *over-warns* S→L, R→C, and Y→C. "Cheap chemistry and curated outcomes disagree here" is a testable claim, not a vibe.

**5. The "truth" is mostly one-star.** The quality breakdown of all clean labels: **1-star 68.2%**, 2-star 29.4%, 3-star 2.4%, and exactly **8 variants** at the 4-star "practice guideline" tier. Only **32%** of curated labels reach ≥2 stars; expert-panel review is vanishingly rare. This is why the review gate matters so much — and why a strict ≥50%-at-≥2-stars rule yields **212** TRUST genes where a naive count-only rule would wave through ~844. Three-quarters of the genes that *look* benchmarkable are leaning on single-submitter labels.

## What I'm not claiming

I want to be precise about the limits, because the strong version of this post would overreach.

This is **a floor, not a verdict on AlphaMissense.** The baseline here is Grantham distance — a deliberately weak, license-clean predictor — so "2 measurable" is a floor for *this* predictor. A stronger model could, in principle, push more genes over the line by reaching a higher-AUROC, lower-variance regime. (That's the sequel.)

The missense definition is **a v1 heuristic.** I infer missense from the HGVS protein token, which over-calls by ~10% versus a consequence-based extraction. The *shape* of the funnel is robust to this; the third decimal of any single AUROC is not.

"Measurable on a gene" is **not "validated across it."** A tight per-gene CI means the AUROC is pin-down-able for the variants ClinVar happens to hold — not that the predictor is right everywhere along the protein.

And the deepest caveat: **ClinVar labels are not predictor-independent.** ACMG's PP3/BP4 criteria explicitly let in-silico predictions count as evidence toward a classification. So when you benchmark a predictor against labels that were themselves partly informed by predictors, you inflate apparent accuracy. This is not a footnote — it's the cleanest justification for the whole design. Asking "is this score correct?" is contaminated. Asking "is there enough trustworthy truth here to benchmark at all?" is the question you can still answer honestly.

## The sequel

All of which sets up one falsifiable question. If the wall is part sample size and part predictor mediocrity — and the math says it's both — then a genuinely strong model should move it. Not by a lot: a better predictor can't manufacture benign labels for the 3,500 genes that have none, and it can't fix the genes gated by tiny sample sizes. But for the handful of genes that are *mediocrity-limited* — healthy on both classes, deep enough, held back only by a so-so ranking — a real model might tighten the interval under 0.10.

The marquee candidate is **TP53**: 240 pathogenic, 151 benign, Grantham AUROC 0.668, CI width 0.11 — a single hair over the line. If a protein language model lifts its ranking toward ~0.85, TP53 becomes the third measurable gene. If it doesn't, the wall was never about the model.

That test — swapping the baseline for ESM-C and asking whether the wall moves — is the next post. The tool, the data, and the five findings here are what it stands on.

---

*The tool, the validated numbers, and the code to reproduce every figure are open: [github.com/crisprking/target-confidence-card](https://github.com/crisprking/target-confidence-card). ClinVar is NIH/NLM public domain; Grantham distances are from Grantham, Science 1974.*
