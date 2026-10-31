"""Optimised branch-statistics engine for CPSS resampling (Paper 8).

Mirrors hidden_files/paper08_work/probes/facade_v2.py exactly in the
statistics used for the selection event (protocol section 2):

  numeric-numeric:     |Pearson r|
  categorical-cat.:    Bergsma (2013) bias-corrected Cramer's V*
  numeric-categorical: max(omega^2, 0)  (Hays)

Categorical columns are pre-encoded once (factorize, NaN -> -1) and
the contingency / group statistics are computed with bincount, so a
half-sample evaluation is a few numpy calls instead of pandas
crosstab / groupby. Levels observed in the subsample only enter the
statistics (matching facade crosstab / unique-level behaviour).

The full-sample strengths, p-values, BH and the paper rule's P3 are
NOT computed here; they use facade_v2 directly.
"""
import numpy as np


class Column:
    """A pre-encoded column. kind is 'numeric' or 'categorical'."""

    def __init__(self, values, kind):
        self.kind = kind
        if kind == "numeric":
            self.x = np.asarray(values, dtype=float)
            self.codes = None
            self.n_levels = 0
        elif kind == "categorical":
            # pd.factorize-like encoding, stable order of first appearance
            seen = {}
            codes = np.empty(len(values), dtype=np.int64)
            for i, v in enumerate(values):
                # values arrive as plain objects/strings; NaN handled by caller
                key = v
                try:
                    if key != key:  # NaN float
                        codes[i] = -1
                        continue
                except Exception:
                    pass
                if key not in seen:
                    seen[key] = len(seen)
                codes[i] = seen[key]
            self.codes = codes
            self.n_levels = len(seen)
            self.x = None
        else:
            raise ValueError(kind)


def _pearson_abs(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    xx, yy = x[m], y[m]
    n = len(xx)
    if n == 0:
        return 0.0
    xc = xx - xx.mean()
    yc = yy - yy.mean()
    denom = np.sqrt((xc * xc).sum() * (yc * yc).sum())
    if denom <= 0:
        return 0.0
    return abs(float((xc * yc).sum() / denom))


def _cramers_v_star(a, b, la, lb):
    m = (a >= 0) & (b >= 0)
    a, b = a[m], b[m]
    n = len(a)
    if n <= 1:
        return 0.0
    ct = np.bincount(a * lb + b, minlength=la * lb).reshape(la, lb).astype(float)
    # Drop unobserved levels, exactly as pd.crosstab on the subsample does.
    row_ok = ct.sum(axis=1) > 0
    col_ok = ct.sum(axis=0) > 0
    ct = ct[np.ix_(row_ok, col_ok)]
    if ct.size == 0 or ct.shape[0] < 2 or ct.shape[1] < 2:
        # facade: k = min(r,c)-1 = 0 -> V = 0, and V* denom <= 0 -> 0
        return 0.0
    r_, c_ = ct.shape
    R = ct.sum(axis=1, keepdims=True)
    C = ct.sum(axis=0, keepdims=True)
    E = R @ C / n
    chi2 = float(((ct - E) ** 2 / np.where(E < 1e-12, 1e-12, E)).sum())
    phi2 = chi2 / n
    phi2_tilde = max(0.0, phi2 - ((r_ - 1) * (c_ - 1)) / (n - 1.0))
    r_tilde = r_ - ((r_ - 1) ** 2) / (n - 1.0)
    c_tilde = c_ - ((c_ - 1) ** 2) / (n - 1.0)
    denom = min(r_tilde - 1.0, c_tilde - 1.0)
    if denom <= 0:
        return 0.0
    return float(np.sqrt(phi2_tilde / denom))


def _omega2_corr(num, cat, n_levels):
    m = np.isfinite(num) & (cat >= 0)
    s, g = num[m], cat[m]
    n = len(s)
    if n == 0:
        return 0.0
    # levels observed in this subsample (after numeric-NaN drop)
    present = np.unique(g)
    cnt = np.bincount(g, minlength=n_levels).astype(float)
    sums = np.bincount(g, weights=s, minlength=n_levels)
    sumsq = np.bincount(g, weights=s * s, minlength=n_levels)
    cnt = cnt[present]
    sums = sums[present]
    grand = s.mean()
    ss_total = float(((s - grand) ** 2).sum())
    # SSB from group means; SSW = SST - SSB (algebraically = facade loop)
    means = sums / cnt
    ss_between = float((cnt * (means - grand) ** 2).sum())
    ss_within = ss_total - ss_between
    glev = len(present)
    ms_within = ss_within / max(n - glev, 1)
    omega2 = (ss_between - (glev - 1) * ms_within) / max(ss_total + ms_within, 1e-300)
    return max(float(omega2), 0.0)


def corrected_abs(col_a, col_b, idx):
    """Corrected strength (absolute for numeric-numeric) of a pair on
    row index array `idx`, using the facade branch statistics."""
    if col_a.kind == "numeric" and col_b.kind == "numeric":
        return _pearson_abs(col_a.x[idx], col_b.x[idx])
    if col_a.kind == "categorical" and col_b.kind == "categorical":
        return _cramers_v_star(col_a.codes[idx], col_b.codes[idx],
                               col_a.n_levels, col_b.n_levels)
    if col_a.kind == "numeric" and col_b.kind == "categorical":
        return _omega2_corr(col_a.x[idx], col_b.codes[idx], col_b.n_levels)
    if col_a.kind == "categorical" and col_b.kind == "numeric":
        return _omega2_corr(col_b.x[idx], col_a.codes[idx], col_a.n_levels)
    raise ValueError((col_a.kind, col_b.kind))


def encode(values, kind):
    """Encode a column, handling missing values for the engine.

    Numeric: NaN preserved. Categorical: None / NaN -> -1 (dropped by
    the statistics, matching crosstab behaviour; facade eta converts
    NaN cats to the string 'nan', but no field/synthetic column here
    contains a missing categorical value, so this case does not arise).
    """
    if kind == "numeric":
        return Column(values, "numeric")
    vals = []
    for v in values:
        if v is None:
            vals.append(np.nan)
        else:
            try:
                if isinstance(v, float) and np.isnan(v):
                    vals.append(np.nan)
                    continue
            except Exception:
                pass
            vals.append(v)
    return Column(vals, "categorical")


def bh_reject(pvals, q=0.05):
    """Identical algorithm to facade_v2.bh_reject (kept here so field
    scripts do not need pandas facade calls in the CPSS loop)."""
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
