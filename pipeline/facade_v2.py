"""facade_v2.py -- Paper 8 production implementation of the type-aware
correlate() facade (numpy + pandas only, no SciPy / no scikit-learn).

Upgrades over the shipped calibration facade, all exercised by the probe
suite in this folder:
  * exact two-sided p-values (regularised incomplete beta for t and F,
    incomplete gamma for chi-square) instead of the normal-tail
    approximation described in the original manuscript;
  * Bergsma (2013) bias-corrected Cramer's V* alongside raw V;
  * Hays omega^2 reported beside eta^2 for the numeric-categorical branch;
  * Benjamini-Hochberg FDR as a pipeline-level component (bh_reject);
  * the persistence test implemented as code (subsample windows).

Special functions follow the standard continued-fraction / series
expansions (Numerical Recipes, betacf / gcf / gser), cross-checked in the
probe suite against closed-form cases.
"""
import hashlib
import numpy as np
import pandas as pd

# ----------------------------------------------------------------------
# special functions
# ----------------------------------------------------------------------
def _betacf(a, b, x, maxit=200, eps=3e-14, fpmin=1e-300):
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < fpmin: d = fpmin
    d = 1.0 / d
    h = d
    for m in range(1, maxit + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin: d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin: c = fpmin
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin: d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin: c = fpmin
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps: break
    return h

def betainc_reg(a, b, x):
    """Regularised incomplete beta I_x(a, b)."""
    if x <= 0.0: return 0.0
    if x >= 1.0: return 1.0
    lbeta = (np.math.lgamma(a + b) - np.math.lgamma(a) - np.math.lgamma(b)
             + a * np.log(x) + b * np.log(1.0 - x))
    bt = np.exp(lbeta)
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    return 1.0 - bt * _betacf(b, a, 1.0 - x) / b

def t_two_sided_p(t, df):
    """Two-sided p for Student t via incomplete beta."""
    x = df / (df + t * t)
    return betainc_reg(df / 2.0, 0.5, x)

def f_sf(f, d1, d2):
    """Survival function of the F distribution."""
    if f <= 0: return 1.0
    x = d2 / (d2 + d1 * f)
    return betainc_reg(d2 / 2.0, d1 / 2.0, x)

def _gser(a, x, eps=3e-14, itmax=500):
    ap = a
    s = 1.0 / a
    delta = s
    for _ in range(itmax):
        ap += 1.0
        delta *= x / ap
        s += delta
        if abs(delta) < abs(s) * eps: break
    return s * np.exp(-x + a * np.log(x) - np.math.lgamma(a))

def _gcf(a, x, eps=3e-14, itmax=500, fpmin=1e-300):
    b = x + 1.0 - a
    c = 1.0 / fpmin
    d = 1.0 / b
    h = d
    for i in range(1, itmax + 1):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < fpmin: d = fpmin
        c = b + an / c
        if abs(c) < fpmin: c = fpmin
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps: break
    return np.exp(-x + a * np.log(x) - np.math.lgamma(a)) * h

def gammaincc(a, x):
    """Regularised upper incomplete gamma Q(a, x) = P(chi2_df > v), df=2a, v=2x."""
    if x < 0 or a <= 0: return float("nan")
    if x == 0: return 1.0
    if x < a + 1.0:
        return 1.0 - _gser(a, x)
    return _gcf(a, x)

def chi2_sf(chi2, df):
    return gammaincc(df / 2.0, chi2 / 2.0)

# ----------------------------------------------------------------------
# analytic branches
# ----------------------------------------------------------------------
def pearson_r(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = ~(np.isnan(x) | np.isnan(y))
    x, y = x[m], y[m]
    n = len(x)
    xc, yc = x - x.mean(), y - y.mean()
    denom = np.sqrt((xc * xc).sum() * (yc * yc).sum())
    r = float((xc * yc).sum() / denom) if denom > 0 else 0.0
    return r, n

def pearson_p(r, n):
    if n < 4 or abs(r) >= 1.0: return 0.0 if abs(r) >= 1.0 else 1.0
    t = r * np.sqrt((n - 2) / max(1.0 - r * r, 1e-300))
    return float(t_two_sided_p(t, n - 2))

def contingency(a, b):
    return pd.crosstab(pd.Series(a).reset_index(drop=True),
                       pd.Series(b).reset_index(drop=True)).values.astype(float)

def cramers_v_from_table(ct):
    n = ct.sum()
    R, C = ct.sum(axis=1, keepdims=True), ct.sum(axis=0, keepdims=True)
    E = R @ C / n
    chi2 = float(((ct - E) ** 2 / np.where(E < 1e-12, 1e-12, E)).sum())
    r_, c_ = ct.shape
    k = min(r_, c_) - 1
    v = float(np.sqrt(chi2 / (n * max(k, 1))))
    df = (r_ - 1) * (c_ - 1)
    return v, chi2, df, n

def cramers_v_star_from_table(ct):
    """Bergsma (2013) bias-corrected Cramer's V."""
    n = ct.sum()
    r_, c_ = ct.shape
    _, chi2, _, _ = cramers_v_from_table(ct)
    phi2 = chi2 / n
    phi2_tilde = max(0.0, phi2 - ((r_ - 1) * (c_ - 1)) / (n - 1.0))
    r_tilde = r_ - ((r_ - 1) ** 2) / (n - 1.0)
    c_tilde = c_ - ((c_ - 1) ** 2) / (n - 1.0)
    denom = min(r_tilde - 1.0, c_tilde - 1.0)
    if denom <= 0: return 0.0
    return float(np.sqrt(phi2_tilde / denom))

def eta_omega(num, cat):
    """Returns eta^2, omega^2 (Hays), F, df1, df2, n, g."""
    s = pd.Series(num).reset_index(drop=True)
    gser = pd.Series(cat).astype(str).reset_index(drop=True)
    m = s.notna()
    s, gser = s[m], gser[m]
    n = len(s)
    grand = s.mean()
    ss_total = float(((s - grand) ** 2).sum())
    ss_between = 0.0
    ss_within = 0.0
    levels = gser.unique()
    g = len(levels)
    for lvl in levels:
        grp = s[gser == lvl]
        ss_between += len(grp) * (grp.mean() - grand) ** 2
        ss_within += float(((grp - grp.mean()) ** 2).sum())
    eta2 = ss_between / max(ss_total, 1e-300)
    ms_within = ss_within / max(n - g, 1)
    omega2 = (ss_between - (g - 1) * ms_within) / max(ss_total + ms_within, 1e-300)
    df1, df2 = g - 1, n - g
    F = (ss_between / max(df1, 1)) / max(ms_within, 1e-300) if df2 > 0 else 0.0
    return float(eta2), float(omega2), float(F), df1, df2, n, g

# ----------------------------------------------------------------------
# text reductions
# ----------------------------------------------------------------------
def hash_embed(text, d=64):
    v = np.zeros(d)
    for tok in str(text).split():
        h = int(hashlib.md5(tok.encode()).hexdigest(), 16) % d
        v[h] += 1.0
    nrm = np.linalg.norm(v)
    return v / nrm if nrm > 0 else v

def tfidf_embed(texts):
    docs = [str(t).lower().split() for t in texts]
    vocab = {}
    for doc in docs:
        for tok in set(doc):
            vocab[tok] = vocab.get(tok, 0) + 1
    terms = sorted(vocab)
    idx = {t: i for i, t in enumerate(terms)}
    N = len(docs)
    idf = np.array([np.log((1 + N) / (1 + vocab[t])) + 1.0 for t in terms])
    X = np.zeros((N, len(terms)))
    for i, doc in enumerate(docs):
        for tok in doc:
            X[i, idx[tok]] += 1.0
        if len(doc): X[i] /= len(doc)
    X = X * idf
    nrm = np.linalg.norm(X, axis=1, keepdims=True)
    return X / np.where(nrm == 0, 1.0, nrm)

def kmeans(X, k=3, n_iter=100, seed=0):
    rng = np.random.RandomState(seed)
    idx = rng.choice(len(X), k, replace=False)
    mu = X[idx].copy()
    z = np.zeros(len(X), dtype=int)
    for _ in range(n_iter):
        d = ((X[:, None, :] - mu[None, :, :]) ** 2).sum(axis=2)
        z_new = d.argmin(axis=1)
        mu_new = np.array([X[z_new == j].mean(axis=0) if (z_new == j).any()
                           else mu[j] for j in range(k)])
        if np.array_equal(z_new, z) and np.allclose(mu, mu_new): 
            z = z_new; break
        z, mu = z_new, mu_new
    return z

# ----------------------------------------------------------------------
# facade
# ----------------------------------------------------------------------
def detect_type(s):
    s = pd.Series(s)
    if pd.api.types.is_numeric_dtype(s): return "numeric"
    vals = s.dropna().astype(str)
    if len(vals) == 0: return "categorical"
    mean_len = vals.str.len().mean()
    uniq_ratio = vals.nunique() / max(len(vals), 1)
    if mean_len > 18 and uniq_ratio > 0.4: return "text"
    return "categorical"

def correlate(a, b, ta=None, tb=None, text_seed=0):
    """Type-aware correlate(). Returns a uniform result record (dict)."""
    a = pd.Series(a).reset_index(drop=True)
    b = pd.Series(b).reset_index(drop=True)
    ta = ta or detect_type(a)
    tb = tb or detect_type(b)
    rec = {"ta": ta, "tb": tb}
    if ta == "text":
        X = np.vstack([hash_embed(t) for t in a])
        a = pd.Series(kmeans(X, k=3, seed=text_seed)); ta = "categorical"
        rec["text_reduced"] = True
    if tb == "text":
        X = np.vstack([hash_embed(t) for t in b])
        b = pd.Series(kmeans(X, k=3, seed=text_seed)); tb = "categorical"
        rec["text_reduced"] = True
    pair = (ta, tb)
    if pair == ("numeric", "numeric"):
        r, n = pearson_r(a, b)
        rec.update(method="pearson_r", strength=r, abs_strength=abs(r),
                   p_value=pearson_p(r, n), n=n)
    elif pair == ("categorical", "categorical"):
        ct = contingency(a, b)
        v, chi2, df, n = cramers_v_from_table(ct)
        vs = cramers_v_star_from_table(ct)
        rec.update(method="cramers_v", strength=v, abs_strength=abs(v),
                   strength_corrected=vs, p_value=float(chi2_sf(chi2, df)),
                   chi2=chi2, df=df, n=int(n))
    elif pair in (("numeric", "categorical"), ("categorical", "numeric")):
        num, cat = (a, b) if ta == "numeric" else (b, a)
        eta2, omega2, F, df1, df2, n, g = eta_omega(num, cat)
        rec.update(method="eta_squared", strength=eta2, abs_strength=abs(eta2),
                   strength_corrected=max(omega2, 0.0),
                   p_value=float(f_sf(F, df1, df2)) if df2 > 0 else 1.0,
                   F=F, df1=df1, df2=df2, n=n, groups=g)
    else:
        raise ValueError(f"unsupported type pair {pair}")
    return rec

# ----------------------------------------------------------------------
# pipeline-level components
# ----------------------------------------------------------------------
def bh_reject(pvals, q=0.05):
    """Benjamini-Hochberg: returns boolean array of rejections at FDR q."""
    p = np.asarray(pvals, float)
    m = len(p)
    order = np.argsort(p)
    ranked = p[order]
    thresh = q * (np.arange(1, m + 1)) / m
    below = ranked <= thresh
    out = np.zeros(m, dtype=bool)
    if below.any():
        kmax = np.max(np.nonzero(below)[0])
        out[order[:kmax + 1]] = True
    return out

def persistence(a, b, ta, tb, full_strength, n_windows=20, frac=0.5,
                ratio=0.5, min_frac=0.60, seed=123):
    """Subsample-window persistence: share of windows whose strength is at
    least `ratio` x the full-sample strength. Pass if share >= min_frac."""
    a = pd.Series(a).reset_index(drop=True)
    b = pd.Series(b).reset_index(drop=True)
    n = len(a)
    rng = np.random.RandomState(seed)
    w = max(int(frac * n), 30)
    hits = 0
    for _ in range(n_windows):
        idx = rng.choice(n, w, replace=False)
        rec = correlate(a.iloc[idx], b.iloc[idx], ta, tb)
        s = rec["strength_corrected"] if "strength_corrected" in rec else rec["strength"]
        base = full_strength
        if abs(s) >= ratio * abs(base): hits += 1
    share = hits / n_windows
    return share, share >= min_frac
