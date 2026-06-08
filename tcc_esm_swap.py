# =============================================================================
# TCC SEQUEL -- does a real model move the measurability wall?
#   score() swap (Grantham -> ESM-C) + resumable cached scoring pass,
#   wired into the IDENTICAL per-gene AUROC + bootstrap-CI funnel.
#
# Self-bootstrapping, same anti-starvation contract as the deep-dive cell:
#   - finds-or-downloads-or-parses ClinVar (cwd or a mounted /kaggle/input ds)
#   - eligibility/TRUST are LABEL-ONLY (predictor-independent) -> computed first
#   - ESM-C is scored ONLY on the eligible genes (the genes where the
#     measurability question actually lives) -- ~50 proteins, not 18,514
#   - scoring is RESUMABLE: a session reset loses at most the in-flight gene
#
# MODE is auto-detected:
#   - "esm"      : torch + a CUDA GPU + the `esm` package are present  -> real
#   - "standin"  : none of the above                                   -> mock
# The mock REPLACES ONLY the GPU forward pass; the offset/log-softmax/vocab/LLR
# read-off, the cache, and the funnel are the SAME code that runs for real.
# Anything the mock prints is plumbing proof, NOT a result.
# =============================================================================
import os, sys, json, gzip, math, time, hashlib, urllib.request, random
import numpy as np

# ----------------------------------------------------------------------------- 0. CONFIG
RNG_SEED          = 1234
N_BOOT            = 1000          # bootstrap resamples per gene (matches funnel)
CI_MAX            = 0.10          # "measurable" iff bootstrap 95% CI width <= this
ELIGIBLE_PER_CLASS= 20            # >=20 P and >=20 B to benchmark a gene
TRUST_MIN_CLEAN   = 40            # >=40 clean labels for TRUST (label-only gate)
TRUST_REVIEW_BAR  = 0.50          # >=50% of clean labels at >=2-star for TRUST
REVIEW_FLOOR      = 1             # keep variants with >= this many review stars
MODEL_NAME        = "esmc_300m"   # or "esmc_600m"
SCORING           = "wt_marginal" # "wt_marginal" (1 pass/protein) | "masked_marginal" (L passes)
TCC_BENCH_MODE    = "trust_nested" # "trust_nested" reproduces the published 2-gene funnel | "all_eligible" = wider >=20/class base
WORKDIR           = os.environ.get("TCC_WORKDIR", ".")

CLINVAR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
F_RAW   = os.path.join(WORKDIR, "variant_summary.txt.gz")
F_VARS  = os.path.join(WORKDIR, "clinvar_variants.jsonl")
F_IDX   = os.path.join(WORKDIR, "genes_index.json")
F_SEQS  = os.path.join(WORKDIR, "esm_sequences.json")     # {gene: {seq, source, acc}}
F_SCORES= os.path.join(WORKDIR, "esm_scores.jsonl")       # append-only {gene,key,llr,score}
F_DONE  = os.path.join(WORKDIR, "esm_done.json")          # ["GENE", ...] genes fully scored
F_FIG   = os.path.join(WORKDIR, "tcc_esm_vs_grantham.png")

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_SET = set(AA)

# Grantham (1974) distances -- the license-free baseline arm.
_GRANTHAM = {  # symmetric upper triangle; filled both ways below
 ('S','R'):110,('S','L'):145,('S','P'):74,('S','T'):58,('S','A'):99,('S','V'):124,
 ('S','G'):56,('S','I'):142,('S','F'):155,('S','Y'):144,('S','C'):112,('S','H'):89,
 ('S','Q'):68,('S','N'):46,('S','K'):121,('S','D'):65,('S','E'):80,('S','M'):135,('S','W'):177,
 ('R','L'):102,('R','P'):103,('R','T'):71,('R','A'):112,('R','V'):96,('R','G'):125,('R','I'):97,
 ('R','F'):97,('R','Y'):77,('R','C'):180,('R','H'):29,('R','Q'):43,('R','N'):86,('R','K'):26,
 ('R','D'):96,('R','E'):54,('R','M'):91,('R','W'):101,
 ('L','P'):98,('L','T'):92,('L','A'):96,('L','V'):32,('L','G'):138,('L','I'):5,('L','F'):22,
 ('L','Y'):36,('L','C'):198,('L','H'):99,('L','Q'):113,('L','N'):153,('L','K'):107,('L','D'):172,
 ('L','E'):138,('L','M'):15,('L','W'):61,
 ('P','T'):38,('P','A'):27,('P','V'):68,('P','G'):42,('P','I'):95,('P','F'):114,('P','Y'):110,
 ('P','C'):169,('P','H'):77,('P','Q'):76,('P','N'):91,('P','K'):103,('P','D'):108,('P','E'):93,
 ('P','M'):87,('P','W'):147,
 ('T','A'):58,('T','V'):69,('T','G'):59,('T','I'):89,('T','F'):103,('T','Y'):92,('T','C'):149,
 ('T','H'):47,('T','Q'):42,('T','N'):65,('T','K'):78,('T','D'):85,('T','E'):65,('T','M'):81,('T','W'):128,
 ('A','V'):64,('A','G'):60,('A','I'):94,('A','F'):113,('A','Y'):112,('A','C'):195,('A','H'):86,
 ('A','Q'):91,('A','N'):111,('A','K'):106,('A','D'):126,('A','E'):107,('A','M'):84,('A','W'):148,
 ('V','G'):109,('V','I'):29,('V','F'):50,('V','Y'):55,('V','C'):192,('V','H'):84,('V','Q'):96,
 ('V','N'):133,('V','K'):97,('V','D'):152,('V','E'):121,('V','M'):21,('V','W'):88,
 ('G','I'):135,('G','F'):153,('G','Y'):147,('G','C'):159,('G','H'):98,('G','Q'):87,('G','N'):80,
 ('G','K'):127,('G','D'):94,('G','E'):98,('G','M'):127,('G','W'):184,
 ('I','F'):21,('I','Y'):33,('I','C'):198,('I','H'):94,('I','Q'):109,('I','N'):149,('I','K'):102,
 ('I','D'):168,('I','E'):134,('I','M'):10,('I','W'):61,
 ('F','Y'):22,('F','C'):205,('F','H'):100,('F','Q'):116,('F','N'):158,('F','K'):102,('F','D'):177,
 ('F','E'):140,('F','M'):28,('F','W'):40,
 ('Y','C'):194,('Y','H'):83,('Y','Q'):99,('Y','N'):143,('Y','K'):85,('Y','D'):160,('Y','E'):122,
 ('Y','M'):36,('Y','W'):37,
 ('C','H'):174,('C','Q'):154,('C','N'):139,('C','K'):202,('C','D'):154,('C','E'):170,('C','M'):196,('C','W'):215,
 ('H','Q'):24,('H','N'):68,('H','K'):32,('H','D'):81,('H','E'):40,('H','M'):87,('H','W'):115,
 ('Q','N'):46,('Q','K'):53,('Q','D'):61,('Q','E'):29,('Q','M'):101,('Q','W'):130,
 ('N','K'):94,('N','D'):23,('N','E'):42,('N','M'):142,('N','W'):174,
 ('K','D'):101,('K','E'):56,('K','M'):95,('K','W'):110,
 ('D','E'):45,('D','M'):160,('D','W'):181,
 ('E','M'):126,('E','W'):152,
 ('M','W'):67,
}
def grantham(a, b):
    if a == b: return 0.0
    return float(_GRANTHAM.get((a, b)) or _GRANTHAM.get((b, a)) or 100.0)

# ----------------------------------------------------------------------------- 1. DATA BOOTSTRAP
def _find(name):
    """Locate a parsed artefact in cwd or any mounted /kaggle/input dataset."""
    if os.path.exists(name):
        return name
    base = os.path.basename(name)
    for root in ("/kaggle/input", WORKDIR):
        if os.path.isdir(root):
            for dp, _, fs in os.walk(root):
                if base in fs:
                    return os.path.join(dp, base)
    return None

_STARS = {  # ReviewStatus -> gold-star count (ClinVar convention)
 "practice guideline": 4,
 "reviewed by expert panel": 3,
 "criteria provided, multiple submitters, no conflicts": 2,
 "criteria provided, single submitter": 1,
 "criteria provided, conflicting classifications": 1,
 "criteria provided, conflicting interpretations": 1,
 "no assertion criteria provided": 0,
 "no assertion provided": 0,
 "no classification provided": 0,
}
def _stars(review_status):
    return _STARS.get((review_status or "").strip().lower(), 0)

import re
_PMAP3 = {  # 3-letter -> 1-letter
 'ala':'A','arg':'R','asn':'N','asp':'D','cys':'C','gln':'Q','glu':'E','gly':'G',
 'his':'H','ile':'I','leu':'L','lys':'K','met':'M','phe':'F','pro':'P','ser':'S',
 'thr':'T','trp':'W','tyr':'Y','val':'V'}
_PROT_RE = re.compile(r'p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})')
_ACC_RE  = re.compile(r'(NP_\d+\.\d+)')

def _parse_missense(name):
    """From a ClinVar Name string -> (ref1, pos, alt1, np_accession) or None.
    Keeps only single-aa missense (excludes synonymous, nonsense, fs, del/ins)."""
    m = _PROT_RE.search(name or "")
    if not m: return None
    r3, pos, a3 = m.group(1).lower(), int(m.group(2)), m.group(3).lower()
    if a3 in ("ter", "*"): return None          # nonsense
    r1, a1 = _PMAP3.get(r3), _PMAP3.get(a3)
    if not r1 or not a1 or r1 == a1: return None # unknown or synonymous
    acc = _ACC_RE.search(name or "")
    return r1, pos, a1, (acc.group(1) if acc else None)

def _sig_to_y(clinsig):
    s = (clinsig or "").lower()
    if "conflicting" in s: return None
    path   = ("pathogenic" in s) and ("benign" not in s)
    benign = ("benign" in s) and ("pathogenic" not in s)
    if path and not benign: return 1
    if benign and not path: return 0
    return None

def build_clinvar_from_raw():
    """Cold path: download variant_summary.txt.gz and parse -> clinvar_variants.jsonl."""
    if not os.path.exists(F_RAW):
        print(f"downloading ClinVar (Internet must be ON)...\n  {CLINVAR_URL}")
        urllib.request.urlretrieve(CLINVAR_URL, F_RAW)
        print(f"done: {os.path.getsize(F_RAW)/1e6:.1f} MB")
    print("parsing variant_summary -> clinvar_variants.jsonl ...")
    n_in = n_out = 0
    cols = None
    with gzip.open(F_RAW, "rt", errors="replace") as fh, open(F_VARS, "w") as out:
        for line in fh:
            if cols is None:
                cols = line.rstrip("\n").split("\t")
                ci = {c: i for i, c in enumerate(cols)}
                need = ("GeneSymbol","Name","ClinicalSignificance","ReviewStatus","Assembly")
                idxs = {k: ci.get(k) for k in need}
                continue
            n_in += 1
            if n_in % 1_000_000 == 0:
                print(f"  ...{n_in:,} rows scanned")
            f = line.rstrip("\n").split("\t")
            if idxs["Assembly"] is not None and f[idxs["Assembly"]] != "GRCh38":
                continue
            y = _sig_to_y(f[idxs["ClinicalSignificance"]])
            if y is None: continue
            pm = _parse_missense(f[idxs["Name"]])
            if pm is None: continue
            ref, pos, alt, acc = pm
            gene = (f[idxs["GeneSymbol"]] or "").split(";")[0].strip()
            if not gene or gene == "-": continue
            st = _stars(f[idxs["ReviewStatus"]])
            if st < REVIEW_FLOOR: continue
            out.write(json.dumps({"gene":gene,"ref":ref,"alt":alt,"pos":pos,
                                  "y":y,"stars":st,"acc":acc}) + "\n")
            n_out += 1
    print(f"parsed {n_out:,} clean >= {REVIEW_FLOOR}-star missense from {n_in:,} rows")

def normalize_record(d):
    """Accept whatever schema is on disk (mine or an earlier cell's) -> canonical."""
    g = d.get("gene") or d.get("GeneSymbol") or d.get("symbol")
    ref = d.get("ref") or d.get("wt") or d.get("from") or d.get("aa_ref")
    alt = d.get("alt") or d.get("mut") or d.get("to")  or d.get("aa_alt")
    pos = d.get("pos") or d.get("position") or d.get("aa_pos") or d.get("protein_pos")
    y   = d.get("y");  y = d.get("label") if y is None else y
    y   = d.get("clin") if y is None else y
    st  = d.get("stars"); st = d.get("review_stars") if st is None else st
    acc = d.get("acc") or d.get("np") or d.get("protein_acc")
    if isinstance(y, str):
        y = _sig_to_y(y)
    try:
        ref, alt = (ref or "").upper(), (alt or "").upper()
        pos = int(pos); y = int(y); st = int(st if st is not None else 0)
    except (TypeError, ValueError):
        return None
    if g and ref in AA_SET and alt in AA_SET and y in (0,1) and pos > 0:
        return {"gene":g,"ref":ref,"alt":alt,"pos":pos,"y":y,"stars":st,"acc":acc}
    return None

def load_clinvar():
    """Warm: read cached jsonl (any schema). Cold: build then read."""
    src = _find(F_VARS)
    if src is None:
        print("parsed files not found -> building from ClinVar this session.")
        build_clinvar_from_raw(); src = F_VARS
    else:
        print(f"using cached parsed file: {src}")
    recs = []
    with open(src) as fh:
        for line in fh:
            line = line.strip()
            if not line: continue
            r = normalize_record(json.loads(line))
            if r: recs.append(r)
    print(f"loaded {len(recs):,} clean missense records over "
          f"{len({r['gene'] for r in recs}):,} genes")
    return recs

# ----------------------------------------------------------------------------- 2. METRICS (tie-correct)
def auroc(scores, labels):
    """Mann-Whitney AUROC with average ranks (tie-correct -- Grantham ties matter)."""
    s = np.asarray(scores, float); y = np.asarray(labels, int)
    n1 = int(y.sum()); n0 = len(y) - n1
    if n1 == 0 or n0 == 0: return float("nan")
    order = np.argsort(s, kind="mergesort")
    ranks = np.empty(len(s), float)
    sr = s[order]; i = 0
    while i < len(sr):
        j = i
        while j + 1 < len(sr) and sr[j+1] == sr[i]:
            j += 1
        ranks[order[i:j+1]] = (i + j) / 2.0 + 1.0  # average rank, 1-based
        i = j + 1
    r1 = ranks[y == 1].sum()
    return (r1 - n1 * (n1 + 1) / 2.0) / (n1 * n0)

def bootstrap_ci(scores, labels, n_boot=N_BOOT, seed=RNG_SEED):
    s = np.asarray(scores, float); y = np.asarray(labels, int)
    pos = np.where(y == 1)[0]; neg = np.where(y == 0)[0]
    rng = np.random.default_rng(seed)
    out = np.empty(n_boot)
    for b in range(n_boot):
        ip = rng.choice(pos, len(pos), replace=True)
        ineg = rng.choice(neg, len(neg), replace=True)
        idx = np.concatenate([ip, ineg])
        out[b] = auroc(s[idx], y[idx])
    lo, hi = np.nanpercentile(out, [2.5, 97.5])
    return float(lo), float(hi), float(hi - lo)

# ----------------------------------------------------------------------------- 3. SCORERS
# A scorer maps a list of unique (gene,pos,ref,alt) keys -> {key: pathogenicity score},
# higher = more pathogenic. The funnel never sees which scorer produced the numbers.

def variant_key(r): return f"{r['gene']}:{r['pos']}{r['ref']}>{r['alt']}"

class GranthamScorer:
    name = "grantham-v1"
    needs_sequence = False
    def score_keys(self, keys, recs_by_key, seqs=None):
        return {k: grantham(recs_by_key[k]['ref'], recs_by_key[k]['alt']) for k in keys}

class ESMCScorer:
    """Real ESM-C scorer. The numeric read-off (_llr_from_logits) is pure numpy
    and unit-tested; only _forward touches torch/GPU/weights. The mock twin
    overrides _forward alone, so this exact read-off path is what gets verified."""
    needs_sequence = True
    def __init__(self, model_name=MODEL_NAME, scoring=SCORING, device="cuda"):
        self.model_name = model_name; self.scoring = scoring; self.device = device
        self.name = f"esmc:{model_name}:{scoring}"
        self._model = None; self._tok = None
        self.aa_to_idx = None; self.bos_offset = 1   # ESMProtein.encode prepends BOS

    def _lazy_load(self):
        if self._model is not None: return
        import torch
        from esm.models.esmc import ESMC
        from esm.tokenization import EsmSequenceTokenizer
        self._torch = torch
        self._model = ESMC.from_pretrained(self.model_name).to(self.device).eval()
        self._tok = EsmSequenceTokenizer()
        self.aa_to_idx = {a: self._tok.convert_tokens_to_ids(a) for a in AA}

    def _forward(self, seq):
        """seq (str) -> log-prob matrix [L+special, vocab] as numpy (BOS at row 0)."""
        self._lazy_load()
        import torch
        from esm.sdk.api import ESMProtein, LogitsConfig
        prot = ESMProtein(sequence=seq)
        with torch.no_grad():
            t = self._model.encode(prot)
            out = self._model.logits(t, LogitsConfig(sequence=True))
            logits = out.logits.sequence[0]                  # [L+special, vocab]
            logp = torch.log_softmax(logits.float(), dim=-1)
        return logp.cpu().numpy()

    def _llr_from_logits(self, logp, pos, ref, alt):
        """Pure: pathogenicity score = -(logP(alt) - logP(ref)) at residue `pos`.
        Higher => model finds the substitution less likely => more pathogenic."""
        row = pos - 1 + self.bos_offset
        ri, ai = self.aa_to_idx[ref], self.aa_to_idx[alt]
        llr = float(logp[row, ai] - logp[row, ri])   # negative when alt less likely
        return -llr                                   # flip so higher = more pathogenic

    def _assert_offset(self, seq, logp):
        """Guard against a silent off-by-one: at most residues the model's own
        top token should *tend* to be the WT residue. We require the mean logP of
        WT residues to beat a uniform baseline -- a weak but catch-all sanity gate."""
        if self.aa_to_idx is None: return
        idxs = [self.aa_to_idx[c] for c in seq if c in self.aa_to_idx]
        rows = [i - 1 + self.bos_offset for i, c in enumerate(seq, 1) if c in self.aa_to_idx]
        wt_lp = float(np.mean([logp[r, j] for r, j in zip(rows, idxs)]))
        if wt_lp < math.log(1.0 / logp.shape[-1]):
            raise RuntimeError(f"offset/vocab sanity failed (mean WT logP {wt_lp:.2f}); "
                               "check bos_offset / aa_to_idx mapping before trusting scores.")

    def score_keys(self, keys, recs_by_key, seqs):
        """Group keys by gene; ONE forward pass per protein (wt_marginal), read off
        every variant in that protein. masked_marginal re-runs per focal position."""
        by_gene = {}
        for k in keys:
            by_gene.setdefault(recs_by_key[k]['gene'], []).append(k)
        out = {}
        for gi, (gene, gkeys) in enumerate(by_gene.items(), 1):
            seq = seqs.get(gene)
            if not seq:
                continue
            if self.scoring == "wt_marginal":
                logp = self._forward(seq)
                if gi == 1: self._assert_offset(seq, logp)
                for k in gkeys:
                    r = recs_by_key[k]
                    if r['pos'] - 1 < len(seq) and seq[r['pos']-1] == r['ref']:
                        out[k] = self._llr_from_logits(logp, r['pos'], r['ref'], r['alt'])
            else:  # masked_marginal: mask the focal residue, then read it off
                positions = sorted({recs_by_key[k]['pos'] for k in gkeys})
                logp_at = {}
                for p in positions:
                    if p - 1 >= len(seq): continue
                    masked = seq[:p-1] + "<mask>" + seq[p:]  # tokenizer handles <mask>
                    logp_at[p] = self._forward(masked)
                for k in gkeys:
                    r = recs_by_key[k]
                    lp = logp_at.get(r['pos'])
                    if lp is not None and r['pos']-1 < len(seq) and seq[r['pos']-1] == r['ref']:
                        out[k] = self._llr_from_logits(lp, r['pos'], r['ref'], r['alt'])
        return out

# --------------------------------------------------------------------- 3b. OFFLINE STANDIN / SELF-TEST SCORERS
# Exist so the GPU read-off path runs with NO GPU and NO weights. None is a result.
#   MockForwardESMC_ReadoffTrue -- overrides _forward ONLY; inherits the real
#     score_keys/_llr_from_logits/_assert_offset; carries NO label signal, so a
#     correct standin run leaves the wall UNMOVED (median AUROC ~0.5). Its job is
#     to prove offset/sign/vocab and the seq[pos-1]==ref gate execute offline.
#   OracleScorer -- LEAKS the label on purpose; proves only that the FUNNEL can
#     detect a moved wall. A moved wall here says nothing about ESM-C.
class MockForwardESMC_ReadoffTrue(ESMCScorer):
    """Synthetic _forward only -> real read-off, no label signal -> AUROC ~0.5."""
    def __init__(self, **kw):
        super().__init__(**kw)
        self.name = "STANDIN-readoff(real-readoff,no-signal)"
        self.aa_to_idx = {a: i for i, a in enumerate(AA)}   # self-consistent 20-aa vocab
        self.bos_offset = 1
        self.n_forward = 0
    def _lazy_load(self): pass
    def _forward(self, seq):
        V = len(AA); L = len(seq); self.n_forward += 1
        logits = np.full((L + self.bos_offset + 1, V), -2.0)
        for i, c in enumerate(seq, 1):
            row = i - 1 + self.bos_offset
            if c in self.aa_to_idx:
                logits[row, self.aa_to_idx[c]] = 2.0          # WT favoured -> clears _assert_offset
            for j in range(V):                                 # per-(pos,aa): NOT shift-invariant
                h = int(hashlib.md5(f"{seq[:12]}|{i}|{j}".encode()).hexdigest(), 16)
                logits[row, j] += ((h % 1000) / 1000.0 - 0.5) * 1.5
        return logits - np.log(np.exp(logits).sum(axis=1, keepdims=True))

class SyntheticSequenceProvider:
    """Offline sequence layer: per-gene sequence carrying each variant's ref residue
    at its position. Drop-in for SequenceProvider in standin mode only."""
    def __init__(self, by_gene, filler="A"):
        self.cache = {}
        for g, rs in by_gene.items():
            seq = [filler] * max(r['pos'] for r in rs)
            for r in rs:
                seq[r['pos'] - 1] = r['ref']
            self.cache[g] = "".join(seq)
    def get(self, gene, np_acc=None, prefer=None):
        return self.cache.get(gene)

class OracleScorer:
    """LEAKS the label (Grantham + 0.6*y + noise). Funnel-sensitivity demo only."""
    name = "ORACLE(funnel-sensitivity)"; needs_sequence = False
    def __init__(self, truth): self.truth = truth
    def score_keys(self, keys, recs_by_key, seqs=None):
        rng = np.random.default_rng(RNG_SEED); out = {}
        for k in keys:
            r = recs_by_key[k]
            out[k] = grantham(r['ref'], r['alt'])/215.0 + 0.6*self.truth.get(k,0) + 0.12*rng.standard_normal()
        return out

# ----------------------------------------------------------------------------- 4. SEQUENCE PROVIDER (resumable)
class SequenceProvider:
    """Per-gene canonical protein sequence, cached to disk and verified against the
    variants (seq[pos-1] == ref). Primary source: UniProt (human, reviewed). A
    RefSeq-by-NP_accession path is provided for callers who need ClinVar's exact
    isoform numbering. Fetching is cached + resumable; offline -> no fetch."""
    def __init__(self, path=F_SEQS, allow_network=True):
        self.path = path; self.allow_network = allow_network
        self.cache = {}
        src = _find(path)
        if src and os.path.exists(src):
            try: self.cache = json.load(open(src))
            except Exception: self.cache = {}

    def _save(self):
        tmp = self.path + ".tmp"
        json.dump(self.cache, open(tmp, "w"))
        os.replace(tmp, self.path)

    def _fetch_uniprot(self, gene):
        url = ("https://rest.uniprot.org/uniprotkb/search?"
               "query=gene_exact:{g}+AND+organism_id:9606+AND+reviewed:true"
               "&fields=accession,sequence&format=json&size=1").format(g=gene)
        with urllib.request.urlopen(url, timeout=30) as r:
            js = json.load(r)
        res = js.get("results") or []
        if not res: return None, None
        acc = res[0].get("primaryAccession")
        seq = (res[0].get("sequence") or {}).get("value")
        return seq, acc

    def _fetch_refseq(self, np_acc):
        url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
               "db=protein&id={a}&rettype=fasta&retmode=text").format(a=np_acc)
        with urllib.request.urlopen(url, timeout=30) as r:
            txt = r.read().decode()
        seq = "".join(l.strip() for l in txt.splitlines() if not l.startswith(">"))
        return seq or None

    def get(self, gene, np_acc=None, prefer="uniprot"):
        if gene in self.cache:
            return self.cache[gene]["seq"]
        if not self.allow_network:
            return None
        seq = acc = None; source = None
        try:
            if prefer == "refseq" and np_acc:
                seq = self._fetch_refseq(np_acc); acc = np_acc; source = "refseq"
            if not seq:
                seq, acc = self._fetch_uniprot(gene); source = "uniprot"
            if not seq and np_acc:
                seq = self._fetch_refseq(np_acc); acc = np_acc; source = "refseq"
        except Exception as e:
            print(f"  [seq fetch failed] {gene}: {e}")
            return None
        if seq:
            self.cache[gene] = {"seq": seq, "source": source, "acc": acc}
            self._save()
        return seq

    @staticmethod
    def verify(seq, recs):
        """Return (kept_recs, match_rate). Drops variants whose stated ref aa does
        not match the fetched sequence at that position -- the isoform/numbering gate."""
        if not seq: return [], 0.0
        kept, ok = [], 0
        for r in recs:
            p = r['pos']
            if 1 <= p <= len(seq) and seq[p-1] == r['ref']:
                kept.append(r); ok += 1
        return kept, (ok / max(1, len(recs)))

# ----------------------------------------------------------------------------- 5. RESUMABLE SCORING PASS
def _cache_tag(scorer):
    """Filesystem-safe tag from the scorer's name so each scorer gets its OWN
    cache files: a standin pass and an ESM-C pass can NEVER contaminate each other
    (no manual rm, no 'flip the GPU on and get the mock's numbers back')."""
    name = getattr(scorer, "name", scorer.__class__.__name__).lower()
    return re.sub(r"_+", "_", re.sub(r"[^0-9a-z]+", "_", name)).strip("_")[:48] or "scorer"

def _load_done():
    src = _find(F_DONE)
    if src and os.path.exists(src):
        try: return set(json.load(open(src)))
        except Exception: return set()
    return set()

def _load_scores():
    src = _find(F_SCORES)
    out = {}
    if src and os.path.exists(src):
        with open(src) as fh:
            for line in fh:
                line = line.strip()
                if not line: continue
                try:
                    d = json.loads(line); out[d["key"]] = d["score"]
                except Exception: pass
    return out

def _append_scores(records):
    with open(F_SCORES, "a") as fh:
        for d in records:
            fh.write(json.dumps(d) + "\n")

def _mark_done(genes_done):
    tmp = F_DONE + ".tmp"
    json.dump(sorted(genes_done), open(tmp, "w"))
    os.replace(tmp, F_DONE)

def run_scoring_pass(scorer, eligible_genes, recs_by_gene, seq_provider,
                     resume=True, fresh=False, mock_truth=None):
    """Score the eligible genes with `scorer`, gene-by-gene, RESUMABLY.
    A reset loses at most the in-flight gene. Returns {key: score}."""
    # Namespace the cache PER SCORER so standin and ESM-C can never collide: each
    # scorer reads/writes only esm_{scores,done}.<scorer-tag>.* . So flipping the
    # GPU on after a standin run starts the ESM-C cache EMPTY (it scores every
    # gene), instead of seeing the mock's genes "already done" and silently
    # handing back chance-level numbers relabelled "ESM-C".
    global F_SCORES, F_DONE
    _tag     = _cache_tag(scorer)
    F_SCORES = os.path.join(WORKDIR, f"esm_scores.{_tag}.jsonl")
    F_DONE   = os.path.join(WORKDIR, f"esm_done.{_tag}.json")
    if fresh:
        for f in (F_SCORES, F_DONE):
            if os.path.exists(f): os.remove(f)
    done   = _load_done() if resume else set()
    scores = _load_scores() if resume else {}
    # (standin mock carries no label signal; OracleScorer takes truth at construction)
    todo = [g for g in eligible_genes if g not in done]
    print(f"\nscoring pass [{scorer.name}]: {len(eligible_genes)} eligible genes "
          f"| {len(done)} cached, {len(todo)} to score")
    t0 = time.time()
    for gi, gene in enumerate(todo, 1):
        recs = recs_by_gene[gene]
        # sequence layer (only needed for sequence-based scorers)
        seq = None
        if scorer.needs_sequence:
            np_acc = next((r.get("acc") for r in recs if r.get("acc")), None)
            seq = seq_provider.get(gene, np_acc=np_acc)
            recs_use, match = SequenceProvider.verify(seq, recs) if seq else (recs, 1.0)
            if seq and match < 0.80:
                print(f"  [skip] {gene}: only {match:.0%} of variants match the fetched "
                      f"sequence (isoform/numbering mismatch)")
                done.add(gene); _mark_done(done); continue
        else:
            recs_use = recs
        keys = sorted({variant_key(r) for r in recs_use})
        recs_by_key = {variant_key(r): r for r in recs_use}
        kscores = scorer.score_keys(keys, recs_by_key, {gene: seq} if seq else {})
        _append_scores([{"gene":gene,"key":k,"score":float(v)} for k,v in kscores.items()])
        scores.update(kscores)
        done.add(gene); _mark_done(done)
        if gi % 10 == 0 or gi == len(todo):
            print(f"  ...{gi}/{len(todo)} genes scored ({time.time()-t0:.1f}s)")
    print(f"scoring pass complete: {len(scores):,} variant-keys scored")
    return scores

# ----------------------------------------------------------------------------- 6. FUNNEL + COMPARISON
def eligible_and_trust(recs):
    """LABEL-ONLY gates (predictor-independent). Returns (eligible, trust, by_gene).
    TCC_BENCH_MODE: "trust_nested" -> eligible is TRUST-nested, so the Grantham
    control reproduces {USH2A, BRCA1} and ESM-added genes stay attributable;
    "all_eligible" -> the wider >=20/class base."""
    by_gene = {}
    for r in recs:
        by_gene.setdefault(r['gene'], []).append(r)
    eligible, trust = [], []
    for g, rs in by_gene.items():
        npos = sum(1 for r in rs if r['y'] == 1)
        nneg = sum(1 for r in rs if r['y'] == 0)
        clean = npos + nneg
        is_elig = npos >= ELIGIBLE_PER_CLASS and nneg >= ELIGIBLE_PER_CLASS
        is_trust = (clean >= TRUST_MIN_CLEAN and npos > 0 and nneg > 0 and
                    sum(1 for r in rs if r['stars'] >= 2) / clean >= TRUST_REVIEW_BAR)
        if is_trust:
            trust.append(g)
        if is_elig and (TCC_BENCH_MODE == "all_eligible" or is_trust):
            eligible.append(g)
    return sorted(eligible), sorted(trust), by_gene

def per_gene_table(eligible_genes, recs_by_gene, scores):
    """For each eligible gene: AUROC + bootstrap CI under the given scores."""
    rows = []
    for g in eligible_genes:
        rs = recs_by_gene[g]
        sc, yy = [], []
        for r in rs:
            k = variant_key(r)
            if k in scores:
                sc.append(scores[k]); yy.append(r['y'])
        if sum(yy) < ELIGIBLE_PER_CLASS or (len(yy)-sum(yy)) < ELIGIBLE_PER_CLASS:
            continue  # lost too many to sequence/ref mismatch
        a = auroc(sc, yy)
        lo, hi, w = bootstrap_ci(sc, yy)
        rows.append({"gene":g,"n":len(yy),"npos":int(sum(yy)),"auroc":a,
                     "lo":lo,"hi":hi,"width":w,"measurable":w <= CI_MAX})
    rows.sort(key=lambda d: d["width"])
    return rows

def compare(recs, grantham_scores, esm_scores, esm_label):
    elig, trust, by_gene = eligible_and_trust(recs)
    gt = per_gene_table(elig, by_gene, grantham_scores)
    et = per_gene_table(elig, by_gene, esm_scores)
    g_meas = [r for r in gt if r["measurable"]]
    e_meas = [r for r in et if r["measurable"]]
    g_by = {r["gene"]: r for r in gt}; e_by = {r["gene"]: r for r in et}

    print("\n" + "="*78)
    print("DOES THE MODEL MOVE THE WALL?  per-gene measurability, fixed sample size")
    print("="*78)
    print(f"eligible genes (>= {ELIGIBLE_PER_CLASS}/class, label-only): {len(elig)}")
    print(f"  GRANTHAM   measurable (CI<= {CI_MAX}): {len(g_meas):>3}  "
          f"-> {', '.join(r['gene'] for r in g_meas) or '(none)'}")
    print(f"  {esm_label:<10} measurable (CI<= {CI_MAX}): {len(e_meas):>3}  "
          f"-> {', '.join(r['gene'] for r in e_meas) or '(none)'}")
    gained = [r['gene'] for r in e_meas if r['gene'] not in {x['gene'] for x in g_meas}]
    lost   = [r['gene'] for r in g_meas if r['gene'] not in {x['gene'] for x in e_meas}]
    print(f"  newly measurable under {esm_label}: {gained or '(none)'}")
    if lost: print(f"  no longer measurable: {lost}")

    common = sorted(set(g_by) & set(e_by), key=lambda g: g_by[g]['width'])
    wg = np.array([g_by[g]['width'] for g in common])
    we = np.array([e_by[g]['width'] for g in common])
    if len(common):
        print(f"\nCI width at fixed n (n={len(common)} genes scored by both):")
        print(f"  median width  GRANTHAM {np.median(wg):.3f}  vs  {esm_label} {np.median(we):.3f}")
        print(f"  genes tighter under {esm_label}: {int((we < wg).sum())}/{len(common)}")
        ag = np.array([g_by[g]['auroc'] for g in common])
        ae = np.array([e_by[g]['auroc'] for g in common])
        print(f"  median per-gene AUROC  GRANTHAM {np.median(ag):.3f}  vs  {esm_label} {np.median(ae):.3f}")

    print("\nper-gene detail (eligible genes, sorted by tightest CI under " + esm_label + "):")
    print(f"  {'gene':<12}{'n':>5}{'npos':>6}   {'GRANTHAM AUROC[CI]':<26}{esm_label+' AUROC[CI]':<26}")
    for g in sorted(set(g_by)|set(e_by),
                    key=lambda g: e_by.get(g,{}).get('width', 9)):
        gr = g_by.get(g); es = e_by.get(g)
        def fmt(r):
            if not r: return "-"
            star = "*" if r["measurable"] else " "
            return f"{r['auroc']:.3f}[{r['lo']:.2f},{r['hi']:.2f}]w{r['width']:.2f}{star}"
        n = (es or gr)['n']; npos = (es or gr)['npos']
        print(f"  {g:<12}{n:>5}{npos:>6}   {fmt(gr):<26}{fmt(es):<26}")
    print(f"  (* = measurable: bootstrap 95% CI width <= {CI_MAX:.2f})")

    print("\nINTERPRETATION (the falsifiable fork this sequel was built to settle):")
    if len(e_meas) > len(g_meas):
        print(f"  -> the wall MOVED: {esm_label} makes {len(e_meas)-len(g_meas)} more gene(s)")
        print(f"     measurable at the SAME labels, by reaching the higher-AUROC /")
        print(f"     lower-variance regime. Better models don't just score better --")
        print(f"     they make more genes measurable.")
    elif len(e_meas) == len(g_meas):
        print(f"  -> the wall held: same measurable count. If CIs didn't tighten, the")
        print(f"     bottleneck is sample size, not model skill -- immovable until")
        print(f"     ClinVar grows more two-sided genes.")
    else:
        print(f"  -> {esm_label} measurable count went DOWN -- inspect ties/sequence loss.")
    return gt, et, common, (g_by, e_by)

# ----------------------------------------------------------------------------- 7. FIGURE
def make_figure(common, g_by, e_by, esm_label):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if not common:
        print("  [figure skipped] no genes scored by both predictors"); return
    wg = [g_by[g]['width'] for g in common]; we = [e_by[g]['width'] for g in common]
    ag = [g_by[g]['auroc'] for g in common]; ae = [e_by[g]['auroc'] for g in common]
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 5.2))
    # Panel A: CI width Grantham vs model (below diagonal => model tighter)
    lim = max(max(wg), max(we)) * 1.05
    ax[0].plot([0, lim], [0, lim], "--", color="#888", lw=1, label="equal width")
    ax[0].axvline(CI_MAX, color="#2a7", lw=1); ax[0].axhline(CI_MAX, color="#2a7", lw=1)
    ax[0].scatter(wg, we, s=42, c="#3b6", edgecolor="white", linewidth=0.6, zorder=3)
    for g in common:
        if g_by[g]['measurable'] or e_by[g]['measurable']:
            ax[0].annotate(g, (g_by[g]['width'], e_by[g]['width']),
                           fontsize=7, xytext=(3,3), textcoords="offset points")
    ax[0].fill_between([0, CI_MAX], 0, CI_MAX, color="#3b6", alpha=0.08)
    ax[0].set_xlim(0, lim); ax[0].set_ylim(0, lim)
    ax[0].set_xlabel("Grantham bootstrap 95% CI width")
    ax[0].set_ylabel(f"{esm_label} bootstrap 95% CI width")
    ax[0].set_title("Does the model tighten the CI?\n(points below the diagonal = yes)")
    ax[0].legend(loc="upper left", fontsize=8, frameon=False)
    # Panel B: per-gene AUROC, Grantham vs model
    lo = min(min(ag), min(ae)) - 0.02; hi = max(max(ag), max(ae)) + 0.02
    ax[1].plot([lo, hi], [lo, hi], "--", color="#888", lw=1)
    ax[1].scatter(ag, ae, s=42, c="#46c", edgecolor="white", linewidth=0.6, zorder=3)
    ax[1].set_xlim(lo, hi); ax[1].set_ylim(lo, hi)
    ax[1].set_xlabel("Grantham per-gene AUROC")
    ax[1].set_ylabel(f"{esm_label} per-gene AUROC")
    ax[1].set_title("Per-gene accuracy, same genes\n(above diagonal = model separates better)")
    fig.suptitle("Sequel test: does a real model move the measurability wall?",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(F_FIG, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote figure -> {F_FIG}")

# ----------------------------------------------------------------------------- 8. MODE + DRIVER
def detect_mode():
    try:
        import torch  # noqa
        if torch.cuda.is_available():
            import esm  # noqa
            return "esm"
    except Exception:
        pass
    return "standin"

def main():
    random.seed(RNG_SEED); np.random.seed(RNG_SEED)
    mode = os.environ.get("TCC_MODE") or detect_mode()
    recs = load_clinvar()
    elig, trust, by_gene = eligible_and_trust(recs)
    print(f"\nlabel-only gates: {len(trust)} TRUST | {len(elig)} eligible (>= "
          f"{ELIGIBLE_PER_CLASS}/class)")

    # baseline arm -- always Grantham, always reproduces the known funnel
    g_scorer = GranthamScorer()
    g_keys = sorted({variant_key(r) for g in elig for r in by_gene[g]})
    g_recs_by_key = {variant_key(r): r for g in elig for r in by_gene[g]}
    g_scores = g_scorer.score_keys(g_keys, g_recs_by_key)
    g_meas = {r['gene'] for r in per_gene_table(elig, by_gene, g_scores) if r['measurable']}
    print(f"Grantham control measurable: {sorted(g_meas)}")

    if mode == "esm":
        assert {"USH2A", "BRCA1"} <= g_meas, \
            f"control regressed: {sorted(g_meas)} -- refusing to spend GPU until the base reproduces {{USH2A, BRCA1}}"
        print("\n>>> MODE = esm  (torch + CUDA + esm present -> REAL ESM-C scoring)")
        scorer = ESMCScorer()
        seqp = SequenceProvider(allow_network=True)
        esm_scores = run_scoring_pass(scorer, elig, by_gene, seqp, resume=True)
        esm_label = "ESM-C"
    else:
        print("\n>>> MODE = standin  (no GPU/esm -> read-off-true MOCK, NO label signal; PLUMBING PROOF ONLY)")
        scorer = MockForwardESMC_ReadoffTrue()
        seqp = SyntheticSequenceProvider(by_gene)
        esm_scores = run_scoring_pass(scorer, elig, by_gene, seqp, resume=True)
        esm_label = "STANDIN"

    gt, et, common, (g_by, e_by) = compare(recs, g_scores, esm_scores, esm_label)
    make_figure(common, g_by, e_by, esm_label)
    if mode != "esm":
        print("\n[!] standin run: the MOCK numbers are synthetic and prove only that the")
        print("    swap point, resumable cache, sequence gate, and funnel comparison")
        print("    are wired correctly. Real ESM-C numbers come from a GPU run.")

    # bundle the shippable artifacts for download + commit (figure + scores + index)
    try:
        import shutil
        shipdir = os.path.join(WORKDIR, "tcc_ship")
        os.makedirs(shipdir, exist_ok=True)
        for f in (F_FIG, F_SCORES, F_IDX):
            if os.path.exists(f): shutil.copy(f, shipdir)
        shutil.make_archive(shipdir, "zip", shipdir)
        print(f"\nbundled shippable artifacts -> {shipdir}.zip  (download from the output panel)")
    except Exception as e:
        print(f"[bundle skipped] {e}")
    return gt, et

if __name__ == "__main__":
    main()
