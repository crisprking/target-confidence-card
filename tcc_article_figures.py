# ===========================================================================
# ARTICLE FIGURES -- one self-bootstrapping cell that writes ALL THREE figures
# for the post: tcc_funnel.png, tcc_substitution_matrix.png,
# tcc_benign_concentration.png.
#
# Order can't starve it: uses clinvar_variants.jsonl + genes_index.json if present
# (cwd or /kaggle/input), else locates-or-downloads ClinVar and parses. The funnel
# panel runs a stratified bootstrap over the ~50 eligible genes, so budget ~1-3 min
# on the real data; the other two figures are fast.
# Needs matplotlib + numpy (preinstalled on Kaggle); stdlib otherwise.
# ===========================================================================
# ===========================================================================
# ARTICLE CELL (self-bootstrapping) -- the funnel + figure that ARE the post.
#
# One cell, no ordering, can't be starved by a session reset:
#   - if clinvar_variants.jsonl + genes_index.json are on disk (cwd OR mounted
#     as a /kaggle/input dataset) it uses them;
#   - else it locates-or-downloads ClinVar's variant_summary.txt.gz and parses
#     it (Internet ON for the first download, ~439 MB, ~2 min);
#   - then tiers every gene, runs the two-axis measurability funnel at the
#     >=1-star floor AND the >=2-star gold standard, writes tcc_funnel.png
#     (the article image) + tcc_provenance.json (the citable record).
#
# Defines its own engine inline -- imports nothing, depends on no prior cell.
# Needs matplotlib (preinstalled on Kaggle); stdlib otherwise.
# ===========================================================================
import csv, gzip, json, os, re, glob, math, random, hashlib, datetime, urllib.request, time
from collections import defaultdict

# ---- config ---------------------------------------------------------------
VARIANTS_PATH = "clinvar_variants.jsonl"
INDEX_PATH    = "genes_index.json"
FIG_PATH      = "tcc_funnel.png"
PROV_PATH     = "tcc_provenance.json"
CLINVAR_URL   = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
ASSEMBLY      = "GRCh38"
PREDICTOR_NAME = "grantham-v1"
RULE_VERSION   = "tcc-verdict-v0.3-reviewstatus"
FORCE_REBUILD  = False     # set True to re-parse ClinVar even if parsed files exist

# ---- locked thresholds (identical to the engine) --------------------------
MIN_LABELED_RATE  = 10
MIN_LABELED_TRUST = 40
REVIEW_BAR        = 0.50
MIN_CLASS         = 20
TARGET_CI         = 0.10
BOOTSTRAP_B       = 2000
STAR_FLOORS       = [1, 2]

AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
       "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
       "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
       "Tyr": "Y", "Val": "V"}
_PTOK = re.compile(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})")

# ===========================================================================
# PART 1 -- data layer: ensure clinvar_variants.jsonl + genes_index.json exist
# ===========================================================================
def _stars(review_status):
    rs = (review_status or "").lower()
    if "practice guideline" in rs: return 4
    if "expert panel" in rs: return 3
    if "no assertion" in rs or "no classification" in rs or "no interpretation" in rs: return 0
    if "criteria provided" in rs:
        if "multiple submitters" in rs and "no conflict" in rs: return 2
        return 1
    return 0

def _bucket(clinsig):
    s = (clinsig or "").lower()
    if "conflicting" in s: return "conflicting"      # check first (contains 'pathogenicity')
    if "pathogenic" in s:  return "plp"
    if "benign" in s:      return "blb"
    if "uncertain" in s:   return "vus"
    return None

def _locate_or_download():
    for p in (glob.glob("/kaggle/input/**/variant_summary.txt.gz", recursive=True)
              + ["variant_summary.txt.gz", "/kaggle/working/variant_summary.txt.gz"]):
        if os.path.exists(p):
            print(f"using raw ClinVar at {p}")
            return p
    print(f"downloading ClinVar (Internet must be ON)...\n  {CLINVAR_URL}")
    urllib.request.urlretrieve(CLINVAR_URL, "variant_summary.txt.gz")
    print(f"  done: {os.path.getsize('variant_summary.txt.gz') / 1e6:.1f} MB")
    return "variant_summary.txt.gz"

def _build_from_clinvar():
    src = _locate_or_download()
    counts = defaultdict(lambda: dict(plp=0, blb=0, vus=0, conflicting=0, ge2star=0))
    audit = defaultdict(int)
    n = emitted = 0; t0 = time.time()
    with gzip.open(src, "rt") as fh, open(VARIANTS_PATH, "w") as out:
        for row in csv.DictReader(fh, delimiter="\t"):
            n += 1
            if n % 1_000_000 == 0: print(f"  ...{n:,} rows scanned", flush=True)
            if (row.get("Assembly") or "") != ASSEMBLY: audit["assembly"] += 1; continue
            if (row.get("Type") or "") != "single nucleotide variant": audit["not_snv"] += 1; continue
            g = (row.get("GeneSymbol") or "").strip()
            if not g or ";" in g or g == "-": audit["intergenic/multigene"] += 1; continue
            m = _PTOK.search(row.get("Name") or "")
            if not m: audit["no_missense_token"] += 1; continue
            ref, alt = AA3.get(m.group(1)), AA3.get(m.group(3))
            if ref is None or alt is None: audit["nonstandard_aa"] += 1; continue
            if ref == alt: audit["synonymous"] += 1; continue
            bk = _bucket(row.get("ClinicalSignificance"))
            if bk is None: audit["clinsig_other"] += 1; continue
            st = _stars(row.get("ReviewStatus"))
            if st < 1: audit["zero_star"] += 1; continue          # >=1-star benchmark floor
            G = g.upper()
            counts[G][bk] += 1
            if bk in ("plp", "blb"):
                if st >= 2: counts[G]["ge2star"] += 1
                out.write(json.dumps({"gene": G, "ref_aa": ref, "alt_aa": alt,
                                      "label": 1 if bk == "plp" else 0, "stars": st}) + "\n")
                emitted += 1
    with open(INDEX_PATH, "w") as fh:
        json.dump({"release": "rebuilt", "assembly": ASSEMBLY,
                   "genes": {g: dict(c) for g, c in counts.items()}}, fh)
    print(f"\nparsed {n:,} rows in {time.time() - t0:.0f}s")
    print(f"wrote {VARIANTS_PATH}  ({emitted:,} scoreable >=1-star P/B missense records)")
    print(f"wrote {INDEX_PATH}     ({len(counts):,} genes)")
    print("parse audit (dropped, by reason):")
    for k, v in sorted(audit.items(), key=lambda x: -x[1]):
        print(f"    {k:22} {v:,}")
    return VARIANTS_PATH, INDEX_PATH

def _find(name):
    if os.path.exists(name):
        return name
    hits = sorted(glob.glob(f"/kaggle/input/**/{name}", recursive=True))
    return hits[0] if hits else None

def ensure_data():
    """Return (variants_path, index_path), building from ClinVar only if needed."""
    if not FORCE_REBUILD:
        v, i = _find(VARIANTS_PATH), _find(INDEX_PATH)
        if v and i:
            print(f"using cached parsed files:\n   {v}\n   {i}")
            return v, i
    print("parsed files not found -> building from ClinVar this session.\n")
    return _build_from_clinvar()

# ===========================================================================
# PART 2 -- the engine (verdict gate + measurability math + Grantham)
# ===========================================================================
def tier_from_counts(plp, blb, ge2star):
    clean = plp + blb
    if clean < MIN_LABELED_RATE:
        return "REFUSE"
    q2 = (ge2star / clean) if clean else 0.0
    if (plp == 0 or blb == 0) or clean < MIN_LABELED_TRUST or q2 < REVIEW_BAR:
        return "CAUTION"
    return "TRUST"

def auroc(scores, labels):
    pairs = sorted(zip(scores, labels), key=lambda t: t[0])
    n = len(pairs); ranks = [0.0] * n; i = 0
    while i < n:
        j = i
        while j + 1 < n and pairs[j + 1][0] == pairs[i][0]:
            j += 1
        avg = (i + j) / 2 + 1
        for k in range(i, j + 1):
            ranks[k] = avg
        i = j + 1
    n_pos = sum(1 for _, y in pairs if y == 1); n_neg = n - n_pos
    if n_pos == 0 or n_neg == 0:
        return None
    sum_pos = sum(r for r, (_, y) in zip(ranks, pairs) if y == 1)
    u = sum_pos - n_pos * (n_pos + 1) / 2
    return u / (n_pos * n_neg)

def is_separable(scores, labels):
    pos = [s for s, y in zip(scores, labels) if y == 1]
    neg = [s for s, y in zip(scores, labels) if y == 0]
    if not pos or not neg:
        return True
    return (max(neg) < min(pos)) or (max(pos) < min(neg))

def _pct(v, p):
    if not v: return float("nan")
    if len(v) == 1: return v[0]
    idx = p * (len(v) - 1); lo = int(idx); frac = idx - lo
    if lo + 1 >= len(v): return v[-1]
    return v[lo] * (1 - frac) + v[lo + 1] * frac

def stratified_ci(scores, labels, B=BOOTSTRAP_B, level=0.95, seed=0):
    rng = random.Random(seed)
    pos = [s for s, y in zip(scores, labels) if y == 1]
    neg = [s for s, y in zip(scores, labels) if y == 0]
    if not pos or not neg:
        return (None, None, None)
    aucs = []; np_, nn_ = len(pos), len(neg)
    for _ in range(B):
        ps = [pos[rng.randrange(np_)] for _ in range(np_)]
        ns = [neg[rng.randrange(nn_)] for _ in range(nn_)]
        a = auroc(ps + ns, [1] * np_ + [0] * nn_)
        if a is not None:
            aucs.append(a)
    aucs.sort(); alpha = (1 - level) / 2
    lo = _pct(aucs, alpha); hi = _pct(aucs, 1 - alpha)
    return (round(lo, 4), round(hi, 4), round(hi - lo, 4))

def classify(scores, labels, seed=0):
    n_pos = sum(1 for y in labels if y == 1); n_neg = len(labels) - n_pos
    n_min = min(n_pos, n_neg)
    if n_min < MIN_CLASS:
        return dict(klass="unmeasurable", auroc=None, ci_width=None, n_min=n_min)
    if is_separable(scores, labels):
        return dict(klass="separable_excluded", auroc=None, ci_width=None, n_min=n_min)
    a = auroc(scores, labels)
    lo, hi, w = stratified_ci(scores, labels, seed=seed)
    if w is not None and w <= TARGET_CI: k = "measurable"
    elif a >= 0.85: k = "label_limited"
    elif a >= 0.70: k = "borderline"
    else: k = "predictor_limited"
    return dict(klass=k, auroc=round(a, 4), ci_width=w, n_min=n_min)

_PROP = {
    "S": (1.42, 9.2, 32),   "R": (0.65, 10.5, 124), "L": (0.0, 4.9, 111),
    "P": (0.39, 8.0, 32.5), "T": (0.71, 8.6, 61),   "A": (0.0, 8.1, 31),
    "V": (0.0, 5.9, 84),    "G": (0.74, 9.0, 3),    "I": (0.0, 5.2, 111),
    "F": (0.0, 5.0, 132),   "Y": (0.20, 6.2, 136),  "C": (2.75, 5.5, 55),
    "H": (0.58, 10.4, 96),  "Q": (0.89, 10.5, 85),  "N": (1.33, 11.6, 56),
    "K": (0.33, 11.3, 119), "D": (1.38, 13.0, 54),  "E": (0.92, 12.3, 83),
    "M": (0.0, 5.7, 105),   "W": (0.13, 5.4, 170),
}
_ALPHA, _BETA, _GAMMA, _RHO = 1.833, 0.1018, 0.000399, 50.723

def grantham(a, b):
    ca, pa, va = _PROP[a]; cb, pb, vb = _PROP[b]
    return _RHO * math.sqrt(_ALPHA * (ca - cb) ** 2 + _BETA * (pa - pb) ** 2 + _GAMMA * (va - vb) ** 2)

def score(ref, alt):
    return grantham(ref, alt)

# ===========================================================================
# PART 3 -- load + the two-axis funnel
# ===========================================================================
def load_variants(path):
    gv = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            r = json.loads(line)
            gv[str(r["gene"]).upper()].append(
                (r["ref_aa"], r["alt_aa"], int(r["label"]), int(r.get("stars", 1)))
            )
    return gv

def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def build_funnel(vpath, ipath):
    gv = load_variants(vpath)
    with open(ipath) as fh:
        idx = json.load(fh)
    genes_idx = idx.get("genes", idx)
    n_all = len(genes_idx)

    def tier_of(g):
        rec = genes_idx.get(g)
        if rec:
            return tier_from_counts(int(rec.get("plp", 0)), int(rec.get("blb", 0)),
                                    int(rec.get("ge2star", rec.get("ge2", 0))))
        v = gv.get(g, [])
        plp = sum(1 for _, _, y, _ in v if y == 1)
        blb = sum(1 for _, _, y, _ in v if y == 0)
        g2 = sum(1 for _, _, _, s in v if s >= 2)
        return tier_from_counts(plp, blb, g2)

    tiers = {g: tier_of(g) for g in genes_idx}
    n_tier = defaultdict(int)
    for t in tiers.values():
        n_tier[t] += 1
    trust = sorted(g for g, t in tiers.items() if t == "TRUST")

    funnel = {}
    for floor in STAR_FLOORS:
        rows, kc = [], defaultdict(int)
        for g in trust:
            v = [(r, a, y) for r, a, y, s in gv.get(g, []) if s >= floor]
            n_pos = sum(1 for _, _, y in v if y == 1)
            n_neg = sum(1 for _, _, y in v if y == 0)
            if n_pos < MIN_CLASS or n_neg < MIN_CLASS:
                continue
            res = classify([score(r, a) for r, a, _ in v], [y for _, _, y in v], seed=0)
            res["gene"] = g
            kc[res["klass"]] += 1
            rows.append(res)
        funnel[floor] = dict(eligible=len(rows), measurable=kc.get("measurable", 0),
                             kc=dict(kc), rows=rows)

    return dict(n_all=n_all, n_tier=dict(n_tier), n_trust=len(trust),
                release=idx.get("release", "unknown"), assembly=idx.get("assembly", "unknown"),
                funnel=funnel, vpath=vpath)

# ===========================================================================
# PART 4 -- report + provenance + figure
# ===========================================================================
KLASS_ORDER = ["measurable", "label_limited", "borderline", "predictor_limited"]

def report(R):
    sha = sha256_of(R["vpath"]); n_records = sum(1 for _ in open(R["vpath"]))
    today = datetime.date.today().isoformat()
    print("\n" + "=" * 74)
    print("TARGET CONFIDENCE CARD -- measurability funnel (citable)")
    print("=" * 74)
    print(f"ClinVar build : release={R['release']}  assembly={R['assembly']}  computed={today}")
    print(f"benchmark file: {R['vpath']}")
    print(f"   sha256     : {sha}")
    print(f"   records    : {n_records:,} (>=1-star clean P/B missense)")
    print(f"predictor     : {PREDICTOR_NAME}   rule: {RULE_VERSION}")
    print(f"thresholds    : TRUST>={MIN_LABELED_TRUST} clean & >={REVIEW_BAR:.0%} 2-star; "
          f"eligible>={MIN_CLASS}/class; measurable CI<={TARGET_CI}")
    print()
    print(f"verdict axis  : {R['n_all']:,} genes -> "
          f"TRUST {R['n_tier'].get('TRUST', 0):,} / "
          f"CAUTION {R['n_tier'].get('CAUTION', 0):,} / "
          f"REFUSE {R['n_tier'].get('REFUSE', 0):,}")
    for floor in STAR_FLOORS:
        f = R["funnel"][floor]; tag = "gold standard" if floor == 2 else "permissive   "
        print(f"\n>= {floor}-star ({tag}):  TRUST {R['n_trust']:,} "
              f"-> eligible {f['eligible']} -> measurable {f['measurable']}")
        for k in KLASS_ORDER:
            if f["kc"].get(k):
                print(f"      {k:18} {f['kc'][k]}")
    m1 = R["funnel"][1]["measurable"]; m2 = R["funnel"][2]["measurable"]
    print("\n" + "-" * 74)
    print(f"HEADLINE: under {PREDICTOR_NAME}, {m1} of {R['n_all']:,} genes "
          f"({m1 / R['n_all'] * 100:.3f}%) are *measurable* at the >=1-star floor; "
          f"{m2} at the >=2-star gold standard.")
    print("-" * 74)
    meas = sorted((r for r in R["funnel"][1]["rows"] if r["klass"] == "measurable"),
                  key=lambda r: -(r["auroc"] or 0))
    if meas:
        print("\nmeasurable genes (>=1-star) -- the ones you can actually benchmark on:")
        print(f"   {'gene':14} {'n_min':>5} {'AUROC':>6} {'CI_width':>9}")
        for r in meas:
            print(f"   {r['gene']:14} {r['n_min']:>5} {r['auroc']:>6.3f} {r['ci_width']:>9.3f}")
    prov = dict(
        tool="target-confidence-card", computed=today,
        clinvar=dict(release=R["release"], assembly=R["assembly"]),
        benchmark_file=dict(path=R["vpath"], sha256=sha, records=n_records),
        predictor=PREDICTOR_NAME, rule_version=RULE_VERSION,
        thresholds=dict(MIN_LABELED_TRUST=MIN_LABELED_TRUST, REVIEW_BAR=REVIEW_BAR,
                        MIN_CLASS=MIN_CLASS, TARGET_CI=TARGET_CI, BOOTSTRAP_B=BOOTSTRAP_B),
        verdict_axis=R["n_tier"],
        funnel={str(fl): dict(trust=R["n_trust"], eligible=R["funnel"][fl]["eligible"],
                              measurable=R["funnel"][fl]["measurable"], classes=R["funnel"][fl]["kc"])
                for fl in STAR_FLOORS},
        measurable_genes=[r["gene"] for r in meas],
    )
    with open(PROV_PATH, "w") as fh:
        json.dump(prov, fh, indent=2)
    print(f"\nwrote provenance -> {PROV_PATH}")

def make_figure(R):
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    INK, MUTE = "#1d2433", "#5b667a"
    C1, C2 = "#2f6f8f", "#16384a"
    KCOL = {"measurable": "#2e8b57", "label_limited": "#d99000",
            "borderline": "#8a6db0", "predictor_limited": "#c0392b"}
    plt.rcParams.update({
        "font.size": 11, "axes.edgecolor": MUTE, "axes.labelcolor": INK,
        "xtick.color": MUTE, "ytick.color": MUTE, "text.color": INK,
        "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 130,
    })
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.5, 5.4),
                                   gridspec_kw=dict(width_ratios=[1.05, 1.0]))

    f1, f2 = R["funnel"][1], R["funnel"][2]
    stages = [
        ("All genes in build", R["n_all"], None),
        ("Trustworthy truth\n(TRUST verdict)", R["n_trust"], None),
        (f"Eligible to benchmark\n(>={MIN_CLASS}/class)", f1["eligible"], f2["eligible"]),
        (f"Measurable\n(CI <= {TARGET_CI})", f1["measurable"], f2["measurable"]),
    ]
    base = 0.7
    y = list(range(len(stages)))[::-1]
    for yi, (lab, v1, v2) in zip(y, stages):
        if v2 is None:
            axA.barh(yi, max(v1 - base, 1e-6), left=base, height=0.62, color=C1, zorder=3)
            axA.text(max(v1, base + 0.2) * 1.18, yi, f"{v1:,}", va="center", ha="left",
                     fontsize=11, color=INK, fontweight="bold")
        else:
            axA.barh(yi + 0.16, max(v1 - base, 1e-6), left=base, height=0.30, color=C1, zorder=3)
            axA.barh(yi - 0.16, max(v2 - base, 1e-6), left=base, height=0.30, color=C2, zorder=3)
            axA.text(max(v1, base + 0.2) * 1.18, yi + 0.16, f"{v1:,}", va="center",
                     ha="left", fontsize=10, color=C1, fontweight="bold")
            axA.text(max(v2, base + 0.2) * 1.18, yi - 0.16, f"{v2:,}", va="center",
                     ha="left", fontsize=10, color=C2, fontweight="bold")
    axA.set_xscale("log")
    axA.set_xlim(base, R["n_all"] * 3)
    axA.set_yticks(y); axA.set_yticklabels([s[0] for s in stages], fontsize=10.5)
    axA.set_xlabel("number of genes  (log scale)")
    axA.set_title("How many genes can actually validate a variant predictor?",
                  fontsize=12.5, fontweight="bold", loc="left", pad=10)
    axA.grid(axis="x", color="#e6e9ee", zorder=0)
    axA.legend(handles=[Patch(color=C1, label=">= 1-star labels"),
                        Patch(color=C2, label=">= 2-star (gold standard)")],
               loc="lower right", frameon=False, fontsize=9.5)

    rows = [r for r in f1["rows"] if r["auroc"] is not None and r["ci_width"] is not None]
    axB.axvspan(0, TARGET_CI, color="#eaf3ee", zorder=0)
    axB.axvline(TARGET_CI, color=KCOL["measurable"], ls="--", lw=1.3, zorder=2)
    axB.axhline(0.5, color=MUTE, ls=":", lw=1.0, zorder=1)
    for k in KLASS_ORDER:
        pts = [(r["ci_width"], r["auroc"]) for r in rows if r["klass"] == k]
        if pts:
            xs, ys = zip(*pts)
            axB.scatter(xs, ys, s=34, color=KCOL[k], alpha=0.82, edgecolor="white",
                        linewidth=0.5, label=k.replace("_", " "), zorder=3)
    for r in rows:
        if r["klass"] == "measurable":
            axB.annotate(r["gene"], (r["ci_width"], r["auroc"]), fontsize=8.5,
                         xytext=(4, 4), textcoords="offset points", color=INK)
    xmax = max([r["ci_width"] for r in rows] + [TARGET_CI * 2]) * 1.08
    axB.set_xlim(0, xmax); axB.set_ylim(0.45, 1.02)
    axB.set_xlabel("bootstrap 95% CI width on AUROC  (narrower = more certain)")
    axB.set_ylabel("per-gene AUROC")
    axB.set_title(f"Why: of {f1['eligible']} eligible genes, {f1['measurable']} land in the\n"
                  f"measurable corner (AUROC certain to +/- {TARGET_CI})",
                  fontsize=12.5, fontweight="bold", loc="left", pad=10)
    axB.text(TARGET_CI / 2, 0.47, "measurable", color=KCOL["measurable"],
             fontsize=9, ha="center", style="italic")
    if rows:
        axB.legend(loc="upper right", frameon=True, framealpha=0.9,
                   edgecolor="#e6e9ee", fontsize=9)

    fig.suptitle(f"Most genes have no ground truth to test a predictor against -- and where "
                 f"they do, {PREDICTOR_NAME} usually can't be pinned down",
                 fontsize=10, color=MUTE, y=1.005)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(FIG_PATH, dpi=150, bbox_inches="tight", facecolor="white")
    print(f"wrote figure     -> {FIG_PATH}")
    try:
        plt.show()
    except Exception:
        pass

MATRIX_FIG = "tcc_substitution_matrix.png"
MIN_CELL   = 25   # min variants in a substitution cell to trust its pathogenic fraction
_AA_ORDER  = ['A','V','L','I','M','F','Y','W','G','P','C','S','T','N','Q','K','R','H','D','E']

def make_matrix_figure(gv):
    from collections import defaultdict as _dd
    tot = _dd(int); pat = _dd(int); n_tot = 0; n_path = 0
    for g in gv:
        for (r, a, y, s) in gv[g]:
            tot[(r, a)] += 1; pat[(r, a)] += y; n_tot += 1; n_path += y
    base = (n_path / n_tot) if n_tot else 0.0
    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.colors import TwoSlopeNorm
    M = np.full((20, 20), np.nan)
    pos = {aa: i for i, aa in enumerate(_AA_ORDER)}
    for (r, a), t in tot.items():
        if t >= MIN_CELL and r in pos and a in pos:
            M[pos[r], pos[a]] = pat[(r, a)] / t
    fig, ax = plt.subplots(figsize=(8.6, 7.4))
    cmap = plt.cm.RdBu_r.copy(); cmap.set_bad("#eceff3")
    vmax = float(np.nanmax(M)) if np.isfinite(M).any() else 1.0
    norm = TwoSlopeNorm(vmin=0.0, vcenter=min(max(base, 1e-3), vmax - 1e-3), vmax=vmax)
    im = ax.imshow(M, cmap=cmap, norm=norm, aspect="equal")
    ax.set_xticks(range(20)); ax.set_xticklabels(_AA_ORDER, fontsize=9)
    ax.set_yticks(range(20)); ax.set_yticklabels(_AA_ORDER, fontsize=9)
    ax.set_xlabel("substituted to (alt)"); ax.set_ylabel("from (ref)")
    ax.set_title(f"Empirical pathogenicity of amino-acid substitutions in ClinVar\n"
                 f"(red = more pathogenic than the {base:.0%} base rate; grey = < {MIN_CELL} variants)",
                 fontsize=11, fontweight="bold")
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("fraction of clean labels that are pathogenic")
    fig.tight_layout()
    fig.savefig(MATRIX_FIG, dpi=150, bbox_inches="tight", facecolor="white")
    print(f"wrote figure -> {MATRIX_FIG}  (base rate {base:.0%}, "
          f"{sum(1 for v in tot.values() if v >= MIN_CELL)} cells >= {MIN_CELL})")
    try: plt.show()
    except Exception: pass

CONC_FIG = "tcc_benign_concentration.png"

def make_concentration_figure(gv, idx):
    def _gini(xs):
        xs = sorted(x for x in xs if x >= 0); n = len(xs); s = sum(xs)
        if n == 0 or s == 0: return float("nan")
        return sum((2 * i - n - 1) * x for i, x in enumerate(xs, 1)) / (n * s)

    def tier_of(g):
        rec = idx.get(g)
        if rec:
            return tier_from_counts(int(rec.get("plp", 0)), int(rec.get("blb", 0)),
                                    int(rec.get("ge2star", rec.get("ge2", 0))))
        v = gv.get(g, [])
        p = sum(1 for _, _, y, _ in v if y == 1); b = sum(1 for _, _, y, _ in v if y == 0)
        return tier_from_counts(p, b, 0)

    blb_all = sorted((int(idx[g].get("blb", 0)) for g in idx))
    n_genes = len(blb_all); tot_b = sum(blb_all)
    Gv = _gini(blb_all)
    cum_x = [i / n_genes for i in range(n_genes + 1)]
    run = 0; cum_y = [0.0]
    for x in blb_all:
        run += x; cum_y.append(run / tot_b if tot_b else 0.0)
    desc = sorted((int(idx[g].get("blb", 0)) for g in idx), reverse=True)
    run = 0; k50 = n_genes
    for i, x in enumerate(desc, 1):
        run += x
        if run >= 0.5 * tot_b: k50 = i; break
    share50 = k50 / n_genes

    trust_pts = []
    for g in idx:
        if tier_of(g) != "TRUST": continue
        v = gv.get(g, [])
        p = sum(1 for _, _, y, _ in v if y == 1); b = sum(1 for _, _, y, _ in v if y == 0)
        if p + b == 0: continue
        cls = "eligible" if (p >= MIN_CLASS and b >= MIN_CLASS) else ("benign_blocked" if b < MIN_CLASS else "path_blocked")
        trust_pts.append((max(p, 1), max(b, 1), cls))

    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    INK, MUTE = "#1d2433", "#5b667a"
    ELIG, BBLK, PBLK = "#2e8b57", "#c0392b", "#2f6f8f"
    plt.rcParams.update({"font.size": 11, "axes.edgecolor": MUTE, "axes.labelcolor": INK,
        "xtick.color": MUTE, "ytick.color": MUTE, "text.color": INK,
        "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 130})
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.5, 5.4))

    axA.plot([0, 1], [0, 1], color=MUTE, ls="--", lw=1.2, zorder=2, label="perfect equality")
    axA.plot(cum_x, cum_y, color=BBLK, lw=2.2, zorder=3, label="benign labels")
    axA.fill_between(cum_x, cum_y, cum_x, color=BBLK, alpha=0.10, zorder=1)
    axA.set_xlim(0, 1); axA.set_ylim(0, 1)
    axA.set_xlabel("cumulative fraction of genes (fewest-benign first)")
    axA.set_ylabel("cumulative fraction of all benign labels")
    axA.set_title(f"Benign labels are hoarded: Gini = {Gv:.2f}", fontsize=12.5,
                  fontweight="bold", loc="left", pad=10)
    axA.annotate(f"top {share50:.0%} of genes\nhold 50% of benign",
                 xy=(1 - share50, 0.5), xytext=(0.30, 0.62), fontsize=9.5, color=INK,
                 arrowprops=dict(arrowstyle="->", color=MUTE, lw=1.0))
    axA.axhline(0.5, color=MUTE, ls=":", lw=0.8, zorder=1)
    axA.legend(loc="upper left", frameon=False, fontsize=9.5)

    order = {"eligible": ("eligible (>=20/class)", ELIG),
             "benign_blocked": ("blocked: too few benign", BBLK),
             "path_blocked": ("blocked: too few pathogenic", PBLK)}
    xmax = max((p for p, _, _ in trust_pts), default=10) * 1.5
    ymax = max((b for _, b, _ in trust_pts), default=10) * 1.5
    axB.set_xscale("log"); axB.set_yscale("log")
    axB.set_xlim(1, xmax); axB.set_ylim(1, ymax)
    axB.add_patch(Rectangle((MIN_CLASS, MIN_CLASS), xmax - MIN_CLASS, ymax - MIN_CLASS,
                            facecolor="#eef4f0", edgecolor="none", zorder=0))
    axB.axvline(MIN_CLASS, color=MUTE, ls="--", lw=1.1, zorder=2)
    axB.axhline(MIN_CLASS, color=MUTE, ls="--", lw=1.1, zorder=2)
    for cls, (lab, col) in order.items():
        pts = [(p, b) for p, b, c in trust_pts if c == cls]
        if pts:
            xs, ys = zip(*pts)
            axB.scatter(xs, ys, s=26, color=col, alpha=0.7, edgecolor="white",
                        linewidth=0.4, label=f"{lab}  (n={len(pts)})", zorder=3)
    axB.set_xlabel("pathogenic labels in gene  (log)")
    axB.set_ylabel("benign labels in gene  (log)")
    axB.set_title("Most TRUST genes have the pathogenic side but not the benign",
                  fontsize=12.5, fontweight="bold", loc="left", pad=10)
    axB.text(MIN_CLASS * 1.15, ymax / 1.25, "eligible\n(>=20 each)", fontsize=8.5,
             color=ELIG, style="italic", va="top", zorder=4)
    axB.legend(loc="lower left", frameon=True, framealpha=0.92, edgecolor="#e6e9ee", fontsize=8.5)

    fig.suptitle("Benign missense is the majority of labels, but it pools in a few genes -- "
                 "so per-gene it is the binding constraint", fontsize=10, color=MUTE, y=1.005)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(CONC_FIG, dpi=150, bbox_inches="tight", facecolor="white")
    print(f"wrote figure -> {CONC_FIG}  (Gini {Gv:.2f}, top {share50:.0%} of genes hold 50% of benign)")
    try: plt.show()
    except Exception: pass

# ===========================================================================
# run -- build/find data once, then write all three article figures
# ===========================================================================
vpath, ipath = ensure_data()

R = build_funnel(vpath, ipath)          # Figure 1: the funnel (+ bootstrap)
make_figure(R)                          # writes tcc_funnel.png

gv = load_variants(vpath)               # load once for figures 2 and 3
with open(ipath) as fh:
    _raw = json.load(fh)
idx = _raw.get("genes", _raw)
make_matrix_figure(gv)                  # Figure 2: tcc_substitution_matrix.png
make_concentration_figure(gv, idx)      # Figure 3: tcc_benign_concentration.png

m1 = R["funnel"][1]["measurable"]; m2 = R["funnel"][2]["measurable"]
print(f"\nall three article figures written. headline: {m1} of {R['n_all']:,} genes "
      f"measurable at >=1-star, {m2} at the >=2-star gold standard.")
