"""
reliability.py — stdlib-only Target Confidence engine (TCC), single-file edition.

Two INDEPENDENT axes, never conflated:
  - Verdict       (predictor-INDEPENDENT): is this gene's ClinVar ground truth
                  trustworthy enough to benchmark a predictor against?
  - Measurability (predictor-DEPENDENT):  given that truth, can we estimate the
                  predictor's per-gene AUROC tightly enough to be useful?

The gate is ClinVar's own review-star quality, NOT the conflicting fraction
(the conflicting-fraction gate was refuted empirically -- it tracked attention,
not reliability, and demoted the expert-panel crown jewels). Do not reintroduce it.

Content-addressed: verdict_hash() is sha256 over (evidence + release + rule version),
so a cited verdict is reproducible and a number change moves the hash.

No third-party deps. No network. No GPU.
"""
from __future__ import annotations

import hashlib
import json
import random
from dataclasses import dataclass, asdict
from enum import Enum
from typing import Optional, Sequence

RULE_VERSION = "tcc-verdict-v0.3-reviewstatus"

# --- locked verdict thresholds (named on purpose) ---
MIN_LABELED_RATE = 10      # < this many clean labels -> REFUSE (nothing to test against)
MIN_LABELED_TRUST = 40     # two-sided & >= this & well-reviewed -> TRUST
REVIEW_BAR = 0.50          # >= this fraction of clean labels must be >= 2-star
LOEUF_INTOLERANT = 0.35    # gnomAD LOEUF < this => LoF-intolerant (CONTEXT hint, not a tier input)

# --- measurability thresholds ---
MIN_CLASS = 20             # need >= this many in EACH class to estimate at all
TARGET_CI = 0.10           # measurable iff stratified bootstrap CI width <= this and not separable


class Tier(str, Enum):
    TRUST = "TRUST"
    CAUTION = "CAUTION"
    REFUSE = "REFUSE"


@dataclass(frozen=True)
class GeneEvidence:
    """Per-gene ClinVar evidence counts. `ge2star` = clean labels at >= 2 review stars."""
    gene: str
    plp: int = 0            # pathogenic / likely pathogenic (clean)
    blb: int = 0            # benign / likely benign (clean)
    vus: int = 0            # uncertain significance: registers the gene, never a label
    conflicting: int = 0    # conflicting classifications: registered, NEVER a quality signal
    ge2star: int = 0        # clean labels at >= 2 stars
    release: str = "unknown"
    loeuf: Optional[float] = None  # gnomAD LOEUF, optional context only

    @property
    def clean(self) -> int:
        return self.plp + self.blb

    @property
    def q2(self) -> float:
        return (self.ge2star / self.clean) if self.clean else 0.0


@dataclass(frozen=True)
class Verdict:
    gene: str
    tier: Tier
    headline: str
    reasons: tuple
    lof_intolerant: Optional[bool]
    hash: str

    def as_dict(self) -> dict:
        d = asdict(self)
        d["tier"] = self.tier.value
        d["reasons"] = list(self.reasons)
        return d


def verdict_hash(ev: GeneEvidence) -> str:
    """Content key over the evidence that drives the tier + the rule version.
    Identical inputs -> identical hash; any count or release change moves it."""
    payload = {
        "gene": ev.gene.upper(), "plp": ev.plp, "blb": ev.blb, "vus": ev.vus,
        "conflicting": ev.conflicting, "ge2star": ev.ge2star,
        "release": ev.release, "rule": RULE_VERSION,
    }
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode()).hexdigest()


def assess(ev: GeneEvidence) -> Verdict:
    """Tier a gene from its evidence. Refusal-first: REFUSE -> CAUTION -> TRUST.
    The tier is predictor-independent: it does NOT move when you swap predictors."""
    clean = ev.clean
    q2 = ev.q2
    reasons: list[str] = []

    if clean < MIN_LABELED_RATE:
        tier = Tier.REFUSE
        reasons.append(f"only {clean} clean labels (< {MIN_LABELED_RATE}); nothing to test against")
    else:
        one_sided = (ev.plp == 0 or ev.blb == 0)
        sparse = clean < MIN_LABELED_TRUST
        weak = q2 < REVIEW_BAR
        if one_sided or sparse or weak:
            tier = Tier.CAUTION
            if one_sided:
                reasons.append(f"one-sided labels (P/LP={ev.plp}, B/LB={ev.blb}) -> AUROC undefined")
            if sparse:
                reasons.append(f"sparse: {clean} clean labels (< {MIN_LABELED_TRUST})")
            if weak:
                reasons.append(f"weakly reviewed: {q2:.0%} of clean labels >= 2-star (< {REVIEW_BAR:.0%})")
        else:
            tier = Tier.TRUST
            reasons.append(f"two-sided, {clean} clean labels, {q2:.0%} >= 2-star (ClinVar 'solid')")

    lof = None if ev.loeuf is None else (ev.loeuf < LOEUF_INTOLERANT)
    if lof:
        reasons.append(f"LoF-intolerant (LOEUF={ev.loeuf} < {LOEUF_INTOLERANT}) -- raises prior cost of a wrong call")

    headline = {
        Tier.TRUST: "ground truth is trustworthy enough to benchmark against",
        Tier.CAUTION: "usable with caveats -- not a clean benchmark",
        Tier.REFUSE: "insufficient trustworthy truth -- refuse",
    }[tier]

    return Verdict(ev.gene.upper(), tier, headline, tuple(reasons), lof, verdict_hash(ev))


# ----------------------------------------------------------------------------
# Index lookup with REFUSE-by-absence
# ----------------------------------------------------------------------------

# A tiny ILLUSTRATIVE index. Real runs load genes_index.json (the pinned
# ClinVar build). Counts here are representative, not authoritative.
DEMO_INDEX: dict[str, dict] = {
    "BRCA1": dict(plp=224, blb=376, vus=3100, conflicting=2639, ge2star=342, loeuf=0.36),
    "TP53":  dict(plp=236, blb=149, vus=900,  conflicting=444,  ge2star=277, loeuf=0.74),
    "LDLR":  dict(plp=731, blb=40,  vus=1200, conflicting=226,  ge2star=470, loeuf=0.51),
    "PTEN":  dict(plp=248, blb=7,   vus=300,  conflicting=129,  ge2star=176, loeuf=0.12),
    "MLH1":  dict(plp=188, blb=36,  vus=700,  conflicting=337,  ge2star=168, loeuf=0.55),
    "PCSK9": dict(plp=22,  blb=18,  vus=260,  conflicting=70,   ge2star=21,  loeuf=0.62),
    "USH2A": dict(plp=232, blb=232, vus=2400, conflicting=335,  ge2star=236, loeuf=0.99),
    "SCN5A": dict(plp=117, blb=31,  vus=1100, conflicting=226,  ge2star=75,  loeuf=0.21),
    # --- synthetic entries (fake symbols) so triage/gate visibly exercise CAUTION + REFUSE ---
    "DEMO_ONESIDED": dict(plp=120, blb=0,  vus=300, conflicting=40, ge2star=90),  # one-sided -> CAUTION
    "DEMO_SPARSE":   dict(plp=12,  blb=10, vus=80,  conflicting=9,  ge2star=14),  # <40 clean -> CAUTION
    "DEMO_THIN":     dict(plp=3,   blb=2,  vus=20,  conflicting=2,  ge2star=4),   # <10 clean -> REFUSE
}


def load_index(path: str = "genes_index.json") -> Optional[dict]:
    """Load a real per-gene index if present in the cwd (Kaggle: /kaggle/working).
    Returns a {GENE: {plp,blb,vus,conflicting,ge2star,loeuf}} dict, or None if absent."""
    try:
        with open(path) as fh:
            raw = json.load(fh)
    except (FileNotFoundError, json.JSONDecodeError):
        return None
    genes = raw.get("genes", raw)  # tolerate {"genes": {...}} or a bare {...}
    out: dict[str, dict] = {}
    for g, rec in genes.items():
        if not isinstance(rec, dict):
            continue
        out[g.upper()] = dict(
            plp=int(rec.get("plp", 0)), blb=int(rec.get("blb", 0)),
            vus=int(rec.get("vus", 0)), conflicting=int(rec.get("conflicting", 0)),
            ge2star=int(rec.get("ge2star", rec.get("ge2", 0))),
            loeuf=rec.get("loeuf"),
        )
    return out or None


def assess_gene(gene: str, index: Optional[dict] = None, release: str = "DEMO") -> Verdict:
    """Look a gene up in `index` and tier it. Absent -> REFUSE BY ABSENCE:
    the engine never guesses about truth that isn't there."""
    idx = index if index is not None else DEMO_INDEX
    key = gene.upper()
    if key not in idx:
        ev = GeneEvidence(key, release=release)
        return Verdict(
            key, Tier.REFUSE, "absent from index -- refuse by absence",
            (f"{key} not in the index ({len(idx)} genes); no truth to vouch for",),
            None, verdict_hash(ev),
        )
    ev = GeneEvidence(gene=key, release=release, **idx[key])
    return assess(ev)


# ----------------------------------------------------------------------------
# Measurability math (pure Python: tie-aware AUROC + stratified bootstrap)
# ----------------------------------------------------------------------------

def auroc(scores: Sequence[float], labels: Sequence[int]) -> Optional[float]:
    """Tie-aware AUROC via the rank-based Mann-Whitney form. Ties -> 0.5 credit.
    Returns None if a class is empty."""
    pairs = sorted(zip(scores, labels), key=lambda t: t[0])
    n = len(pairs)
    # average ranks (1-based), tie-aware
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and pairs[j + 1][0] == pairs[i][0]:
            j += 1
        avg = (i + j) / 2 + 1  # average of 1-based ranks i+1..j+1
        for k in range(i, j + 1):
            ranks[k] = avg
        i = j + 1
    n_pos = sum(1 for _, y in pairs if y == 1)
    n_neg = n - n_pos
    if n_pos == 0 or n_neg == 0:
        return None
    sum_pos = sum(r for r, (_, y) in zip(ranks, pairs) if y == 1)
    u = sum_pos - n_pos * (n_pos + 1) / 2
    return u / (n_pos * n_neg)


def is_separable(scores: Sequence[float], labels: Sequence[int]) -> bool:
    """Perfect separation: every positive scores strictly above every negative
    (or vice versa). Bootstrap CIs collapse to ~0 here -> exclude, don't trust."""
    pos = [s for s, y in zip(scores, labels) if y == 1]
    neg = [s for s, y in zip(scores, labels) if y == 0]
    if not pos or not neg:
        return True
    return (max(neg) < min(pos)) or (max(pos) < min(neg))


def _percentile(sorted_vals: list[float], p: float) -> float:
    if not sorted_vals:
        return float("nan")
    if len(sorted_vals) == 1:
        return sorted_vals[0]
    idx = p * (len(sorted_vals) - 1)
    lo = int(idx)
    frac = idx - lo
    if lo + 1 >= len(sorted_vals):
        return sorted_vals[-1]
    return sorted_vals[lo] * (1 - frac) + sorted_vals[lo + 1] * frac


def stratified_ci(scores: Sequence[float], labels: Sequence[int],
                  B: int = 2000, level: float = 0.95, seed: int = 0):
    """STRATIFIED percentile bootstrap: resample WITHIN each class so class sizes
    are preserved (pooled resampling lets the minority count drift -> biased AND
    degenerate). Returns (lo, hi, width) or (None, None, None)."""
    rng = random.Random(seed)
    pos = [s for s, y in zip(scores, labels) if y == 1]
    neg = [s for s, y in zip(scores, labels) if y == 0]
    if not pos or not neg:
        return (None, None, None)
    aucs: list[float] = []
    np_, nn_ = len(pos), len(neg)
    for _ in range(B):
        ps = [pos[rng.randrange(np_)] for _ in range(np_)]
        ns = [neg[rng.randrange(nn_)] for _ in range(nn_)]
        a = auroc(ps + ns, [1] * np_ + [0] * nn_)
        if a is not None:
            aucs.append(a)
    aucs.sort()
    alpha = (1 - level) / 2
    lo = _percentile(aucs, alpha)
    hi = _percentile(aucs, 1 - alpha)
    return (round(lo, 4), round(hi, 4), round(hi - lo, 4))


def classify_measurability(scores: Sequence[float], labels: Sequence[int],
                           seed: int = 0) -> dict:
    """The v2 measurability bridge. Excludes separable, then 3-way splits the rest:

      separable_excluded  perfect separation; CI is fake-tight
      unmeasurable        < MIN_CLASS in some class; AUROC is undefined-by-noise
      measurable          CI width <= TARGET_CI               -> trust the number
      label_limited       AUROC high but CI wide               -> predictor works, need more labels
      borderline          AUROC moderate, CI wide              -> ambiguous
      predictor_limited   AUROC near chance, CI wide           -> the predictor itself is weak
    """
    n_pos = sum(1 for y in labels if y == 1)
    n_neg = len(labels) - n_pos
    n_min = min(n_pos, n_neg)
    base = dict(n_pos=n_pos, n_neg=n_neg, n_min=n_min,
                auroc=None, ci=None, ci_width=None)

    if n_min < MIN_CLASS:
        return {**base, "klass": "unmeasurable",
                "why": f"< {MIN_CLASS} in some class (min={n_min}); AUROC is noise"}
    if is_separable(scores, labels):
        return {**base, "klass": "separable_excluded",
                "why": "perfect separation -> bootstrap CI collapses to ~0 (not real stability)"}

    a = auroc(scores, labels)
    lo, hi, w = stratified_ci(scores, labels, seed=seed)
    out = {**base, "auroc": round(a, 4), "ci": [lo, hi], "ci_width": w}

    if w is not None and w <= TARGET_CI:
        out["klass"], out["why"] = "measurable", f"CI width {w} <= {TARGET_CI}"
    elif a >= 0.85:
        out["klass"], out["why"] = "label_limited", f"AUROC {a:.3f} high but CI width {w} > {TARGET_CI}"
    elif a >= 0.70:
        out["klass"], out["why"] = "borderline", f"AUROC {a:.3f} moderate, CI width {w} > {TARGET_CI}"
    else:
        out["klass"], out["why"] = "predictor_limited", f"AUROC {a:.3f} near chance"
    return out


def _draw(n_pos: int, n_neg: int, sep: float, noise: float, seed: int):
    """Synthetic predictor scores for a gene: positives ~ N(sep, noise),
    negatives ~ N(0, noise). Larger sep -> higher AUROC; smaller n -> wider CI.
    A stand-in so the file runs with no GPU/ESM; swap in real {gene:(labels,scores)}."""
    rng = random.Random(seed)
    pos = [rng.gauss(sep, noise) for _ in range(n_pos)]
    neg = [rng.gauss(0.0, noise) for _ in range(n_neg)]
    return pos + neg, [1] * n_pos + [0] * n_neg


def grantham_selfcheck() -> None:
    """Placeholder so a notebook self-check line has something to call."""
    # The real engine ships a full Grantham matrix; the single-file edition
    # uses synthetic predictor scores (_draw), so there is nothing to verify here.
    return None


__all__ = [
    "Tier", "GeneEvidence", "Verdict", "assess", "verdict_hash",
    "assess_gene", "load_index", "DEMO_INDEX",
    "auroc", "is_separable", "stratified_ci", "classify_measurability", "_draw",
    "RULE_VERSION", "MIN_LABELED_RATE", "MIN_LABELED_TRUST", "REVIEW_BAR",
    "MIN_CLASS", "TARGET_CI",
]



# ---------------------------------------------------------------------------
# CLI  —  python reliability.py [GENE ...] [--index genes_index.json] [--json]
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(
        description="Target Confidence Card: is a gene's ClinVar ground truth trustworthy "
                    "enough to benchmark a variant predictor against?  TRUST / CAUTION / REFUSE.")
    ap.add_argument("genes", nargs="*",
                    help="gene symbol(s) to assess; omit to tour the built-in demo index")
    ap.add_argument("--index", default=None,
                    help="genes_index.json from build_index.py (omit to use the built-in demo)")
    ap.add_argument("--json", action="store_true", help="emit machine-readable JSON")
    args = ap.parse_args()

    idx = load_index(args.index) if args.index else None
    if args.index and idx is None:
        print(f"(index {args.index!r} not found -- using the built-in demo index)\n")
    release = "DEMO" if idx is None else args.index
    source = idx if idx is not None else DEMO_INDEX
    targets = args.genes if args.genes else sorted(source)

    verdicts = [assess_gene(g, idx, release=release) for g in targets]
    if args.json:
        print(json.dumps([v.as_dict() for v in verdicts], indent=2))
    else:
        for v in verdicts:
            print(f"{v.gene:14} {v.tier.value}")
            print(f"   {v.headline}")
            for r in v.reasons:
                print(f"     - {r}")
            print(f"   hash {v.hash[:16]}\n")
