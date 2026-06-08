# Target Confidence Card (TCC)

**A refusal-first reliability prior for missense variant-effect predictors.**

Before you trust *any* missense predictor (AlphaMissense, ESM, EVE, REVEL, …) on a gene, TCC
answers a prior question it actually *can* answer from public data:

> Does this gene have enough trustworthy, two-sided ClinVar ground truth to know whether
> any predictor can be benchmarked here at all?

It returns **TRUST / CAUTION / REFUSE**, and it **never rescores a predictor** — so the verdict
does not move when you swap models. Every verdict is content-addressed (SHA-256 over the evidence
+ ClinVar release + rule version), so a cited verdict is reproducible and any count change moves
the hash.

---

## The headline (validated on real ClinVar — GRCh38, 2026-06 release)

Of **18,514** genes with curated missense variants, the funnel collapses to almost nothing:

```
18,514  genes with curated missense
   ->  212  TRUST       two-sided, >=40 clean labels, >=50% of them at >=2 review stars
   ->   50  eligible    >=20 labels in EACH class (enough to even estimate an AUROC)
   ->    2  MEASURABLE  bootstrap 95% CI width <= 0.10   ->   USH2A, BRCA1
   ->    0  at the >=2-star "gold standard"
```

**"Measurable" is not "good."** The two genes that clear sit at AUROC **0.77** and **0.72**. The
point is not that they are well predicted — it is that for ~18,512 other genes we cannot even *say*
how good a predictor is, because the trustworthy ground truth to check it against isn't there.

Full numbers and the five supporting findings: [`docs/VALIDATED_RESULTS.md`](docs/VALIDATED_RESULTS.md).
The write-up: [`docs/two_genes.md`](docs/two_genes.md).

---

## Two axes it never conflates

| Axis | Depends on the predictor? | The question |
|------|:---:|------|
| **Verdict** — `assess_gene()` | **No** | Is this gene's ClinVar truth trustworthy enough to benchmark against? |
| **Measurability** — `classify_measurability()` | **Yes** | Given that truth, can we pin the predictor's per-gene AUROC tightly (CI width ≤ 0.10)? |

A gene can be TRUST (good truth) yet still not measurable for a given predictor if there aren't
enough labels per class — and the same predictor flips between measurable and not purely on gene size.

---

## Install

Python 3.10+. The engine (`reliability.py`) and `quickstart.py` are **pure standard library** —
nothing to install to get a verdict.

```bash
git clone https://github.com/crisprking/target-confidence-card
cd target-confidence-card
python quickstart.py          # ~5 seconds, no network, no GPU
```

`pip install -r requirements.txt` (numpy, matplotlib) is needed only to reproduce the figures
or run the ESM-C sequel.

---

## Use it

**On the built-in demo index (zero setup):**

```bash
python reliability.py BRCA1 TP53 PCSK9      # TRUST/CAUTION/REFUSE + reasons + hash
python reliability.py --json BRCA1          # machine-readable
```

**On real, current ClinVar (any gene):** build the index once, then query offline.

```bash
# 1) ClinVar's public bulk file (NIH/NLM, public domain), ~150 MB
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz

# 2) build a per-gene index: counts P/LP, B/LB, VUS, conflicting, and >=2-star labels
python build_index.py --clinvar variant_summary.txt.gz --out genes_index.json --snapshot 2026-06

# 3) query any gene against your dated snapshot
python reliability.py --index genes_index.json SCN1A
```

**In Python:**

```python
import reliability as tcc
idx = tcc.load_index("genes_index.json")    # or None to use the built-in demo index
v = tcc.assess_gene("BRCA1", idx)
print(v.tier.value, v.hash)                  # e.g. 'TRUST' '3f0a…'
print(v.as_dict())                           # gene, tier, headline, reasons, lof_intolerant, hash
```

---

## Reproduce the paper figures

```bash
python tcc_article_figures.py     # rebuilds the index from ClinVar, writes the three figures
```

Writes `tcc_funnel.png`, `tcc_substitution_matrix.png`, `tcc_benign_concentration.png` (move them
into [`figures/`](figures/)). The script is self-bootstrapping — it downloads and parses ClinVar if
the parsed files aren't present — so run order can't starve it.

---

## The sequel — apparatus, not yet a result

`tcc_esm_swap.py` swaps the funnel's `score()` for a real protein language model (ESM-C) to ask
whether a stronger predictor *moves the measurability wall* (tightens CIs at fixed n by reaching a
higher-AUROC, lower-variance regime — or whether the wall is pure sample size).

**No real ESM-C numbers exist yet.** Running the model needs a GPU and the EvolutionaryScale `esm`
package; without them the module auto-falls-back to a **no-signal standin** whose output sits at
chance by construction. Do **not** read the standin's table as a result. Details and status:
[`docs/VALIDATED_RESULTS.md`](docs/VALIDATED_RESULTS.md).

---

## What this is — and what it is not

- **A floor, not a verdict on AlphaMissense.** The bundled baseline is Grantham distance (1974) — a
  deliberately weak, license-clean predictor — so "2 measurable" is a floor *for this predictor*, not
  a judgment of any specific modern model.
- **Missense extraction is a v1 heuristic.** `build_index.py` infers missense from the HGVS protein
  token, which over-calls ~10% versus a consequence-based extraction. The funnel *shape* is robust; the
  third decimal isn't.
- **"Measurable on a gene" ≠ "validated across it."** A tight per-gene CI means the AUROC is
  pin-down-able, not that the predictor is correct everywhere in that gene.
- **ClinVar labels are not predictor-independent.** ACMG's PP3/BP4 criteria let in-silico predictions
  count as evidence toward a classification, so the labels are partly predictor-informed. This is
  precisely why the tool measures *whether you can benchmark* rather than asserting a predictor is right.

---

## Data & credit

ClinVar (Landrum et al., *Nucleic Acids Research*) — NIH/NLM, public domain. Grantham distances —
Grantham, *Science*, 1974. Optional constraint: gnomAD LOEUF (open). Review-status star ratings follow
ClinVar's own scheme. Code: **MIT** (see [`LICENSE`](LICENSE)).
