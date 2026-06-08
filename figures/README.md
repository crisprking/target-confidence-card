# Figures

These three figures are referenced by `../README.md` and `../docs/two_genes.md`.
They are **reproducible from public ClinVar**, so they are not committed as binaries
that could drift from the snapshot.

Generate them (Internet on for the first run to fetch ClinVar; ~2-3 min total, almost
all of it the funnel's per-gene bootstrap):

```bash
python ../tcc_article_figures.py     # writes the three PNGs to the working directory
```

then place the outputs here:

- `tcc_funnel.png` — the 18,514 -> 212 -> 50 -> 2 measurability funnel (0 at the >=2-star gold standard)
- `tcc_substitution_matrix.png` — the empirical 20x20 pathogenic-fraction-per-substitution map
- `tcc_benign_concentration.png` — Lorenz curve + per-gene class-count scatter (Gini 0.70)

If you have already generated them, just drop the PNGs into this folder before committing.
