"""Probe 8 -- family error-control certificate for the conjoint promotion
rule of Paper 8 (Definition 1: P1 BH q=0.05 AND P2 corrected-strength floor
tau=0.20 AND P3 split-half: corrected strength >= tau on each half).

Sections:
  0. Validate the reference statistics against facade_v2.correlate (the
     production definitions), and the batched world engine against the
     reference statistics.
  A. Monte Carlo calibration of per-half exceedance gamma_tau at half
     sample sizes m = n/2 for total n in {60, 100, 200, 500, 1000}, by
     branch cell, with one-sided Clopper-Pearson upper bounds (per-cell
     95% and simultaneous 95% via Bonferroni over all cells).
  B. All-null world validation (probe-2 design: 40 columns, 780 pairs) at
     n in {60, 100, 500}: full conjoint rule vs the certificate.
  C. Planted worlds (probe-2 planted design, 6 plants) at n = 500.
  D. Cluster stress test: clustered rows, random row split vs
     cluster-level split, at n = 200.

Corrected-strength definitions match facade_v2.py exactly:
  numeric-numeric : |r| (Pearson)
  categorical-cat.: Bergsma V* (cramers_v_star_from_table)
  numeric-cat.    : max(omega^2, 0) (Hays)
p-values: numeric-num via the exact t pivot, cat-cat via chi^2,
numeric-cat via F (scipy.stats tails; validated against facade_v2's
own special functions in Section 0).

NOTE on the validation import: facade_v2 uses np.math.lgamma; the
archived venv ships numpy 2.5, which removed the np.math alias. We
restore it as the stdlib math module (np.math was exactly that alias)
before importing facade_v2. facade_v2.py itself is not modified.

Run with the archived venv python:
  hidden_files/paper08_work/phik_venv/bin/python probe8_certificate.py
All outputs are JSON files written next to this script.
"""
import json
import math
import os
import sys
import time

import numpy as np

# --- numpy-2 shim for facade_v2 validation only (see docstring) ---
if not hasattr(np, "math"):
    np.math = math

from scipy import stats as st
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
PROBES = os.path.join(HERE, "..", "..", "paper08_work", "probes")
sys.path.insert(0, PROBES)

TAU = 0.20
Q_BH = 0.05
DELTA = 0.05  # simultaneous miscoverage budget for calibration bounds

T_START = time.time()


def stamp(section, t0):
    print(f"[{section}] elapsed {time.time() - t0:.1f}s "
          f"(total {time.time() - T_START:.1f}s)", flush=True)


# ----------------------------------------------------------------------
# Reference statistics (definitions identical to facade_v2; used for
# validation and as the slow path for degenerate tables)
# ----------------------------------------------------------------------
def pearson_r_vec(x, y):
    xc = x - x.mean()
    yc = y - y.mean()
    denom = math.sqrt(float((xc * xc).sum() * (yc * yc).sum()))
    return float((xc * yc).sum() / denom) if denom > 0 else 0.0


def pearson_p_from_r(r, n):
    if n < 4:
        return 1.0
    if abs(r) >= 1.0:
        return 0.0
    t = r * math.sqrt((n - 2) / max(1.0 - r * r, 1e-300))
    return float(2.0 * st.t.sf(abs(t), n - 2))


def vstar_and_p_from_counts(ct):
    """Bergsma V* and chi^2 p from a contingency table, after dropping
    zero-margin rows/cols (matches facade_v2's pd.crosstab)."""
    ct = np.asarray(ct, dtype=float)
    ct = ct[ct.sum(axis=1) > 0][:, ct.sum(axis=0) > 0]
    n = ct.sum()
    if n == 0 or ct.shape[0] < 2 or ct.shape[1] < 2:
        return 0.0, 1.0
    r_, c_ = ct.shape
    E = ct.sum(axis=1, keepdims=True) @ ct.sum(axis=0, keepdims=True) / n
    chi2 = float(((ct - E) ** 2 / np.where(E < 1e-12, 1e-12, E)).sum())
    df = (r_ - 1) * (c_ - 1)
    p = float(st.chi2.sf(chi2, df)) if df > 0 else 1.0
    phi2_t = max(0.0, chi2 / n - ((r_ - 1) * (c_ - 1)) / (n - 1.0))
    r_t = r_ - ((r_ - 1) ** 2) / (n - 1.0)
    c_t = c_ - ((c_ - 1) ** 2) / (n - 1.0)
    denom = min(r_t - 1.0, c_t - 1.0)
    vs = float(math.sqrt(phi2_t / denom)) if denom > 0 else 0.0
    return vs, p


def cat_pair_stats(a, b, ra, rb):
    comb = a * rb + b
    ct = np.bincount(comb, minlength=ra * rb).reshape(ra, rb)
    return vstar_and_p_from_counts(ct)


def numcat_stats(x, codes, g):
    """Returns (omega^2 unclipped, F p)."""
    n = len(x)
    counts = np.bincount(codes, minlength=g).astype(float)
    sums = np.bincount(codes, weights=x, minlength=g)
    obs = counts > 0
    g_obs = int(obs.sum())
    if g_obs < 2 or n - g_obs <= 0:
        return 0.0, 1.0
    grand = x.mean()
    sst = float(((x - grand) ** 2).sum())
    if sst <= 0:
        return 0.0, 1.0
    ssb = float((sums[obs] ** 2 / counts[obs]).sum() - n * grand * grand)
    ssw = max(sst - ssb, 0.0)
    msw = ssw / (n - g_obs)
    omega2 = (ssb - (g_obs - 1) * msw) / (sst + msw)
    F = (ssb / (g_obs - 1)) / msw if msw > 0 else 0.0
    p = float(st.f.sf(F, g_obs - 1, n - g_obs)) if F > 0 else 1.0
    return float(omega2), p


def bh_reject(pvals, q=Q_BH):
    """Identical to facade_v2.bh_reject."""
    p = np.asarray(pvals, float)
    m = len(p)
    order = np.argsort(p)
    ranked = p[order]
    thresh = q * (np.arange(1, m + 1)) / m
    below = ranked <= thresh
    out = np.zeros(m, dtype=bool)
    if below.any():
        kmax = int(np.max(np.nonzero(below)[0]))
        out[order[:kmax + 1]] = True
    return out


# ----------------------------------------------------------------------
# Batched world engine (strengths identical to the reference statistics;
# validated in Section 0)
# ----------------------------------------------------------------------
def cat_strength_chi2(a, b, ra, rb):
    """(V*, chi2, df) for one cat-cat pair, strength-only arithmetic."""
    comb = a * rb + b
    ct = np.bincount(comb, minlength=ra * rb).reshape(ra, rb).astype(float)
    rs = ct.sum(axis=1); cs = ct.sum(axis=0)
    if (rs == 0).any() or (cs == 0).any():
        ct = ct[rs > 0][:, cs > 0]
        if ct.shape[0] < 2 or ct.shape[1] < 2:
            return 0.0, 0.0, 1
        rs = ct.sum(axis=1); cs = ct.sum(axis=0)
    n = ct.sum()
    r_, c_ = ct.shape
    E = np.outer(rs, cs) / n
    chi2 = float(((ct - E) ** 2 / np.where(E < 1e-12, 1e-12, E)).sum())
    df = (r_ - 1) * (c_ - 1)
    phi2_t = max(0.0, chi2 / n - ((r_ - 1) * (c_ - 1)) / (n - 1.0))
    r_t = r_ - ((r_ - 1) ** 2) / (n - 1.0)
    c_t = c_ - ((c_ - 1) ** 2) / (n - 1.0)
    denom = min(r_t - 1.0, c_t - 1.0)
    vs = float(math.sqrt(phi2_t / denom)) if denom > 0 else 0.0
    return vs, chi2, df


def numcat_strength_F(x, codes, g):
    """(max(omega^2,0), F, df1, df2) for one numeric-categorical pair."""
    n = len(x)
    counts = np.bincount(codes, minlength=g).astype(float)
    sums = np.bincount(codes, weights=x, minlength=g)
    obs = counts > 0
    g_obs = int(obs.sum())
    if g_obs < 2 or n - g_obs <= 0:
        return 0.0, 0.0, 1, max(n - 2, 1)
    grand = float(x.mean())
    sst = float(((x - grand) ** 2).sum())
    if sst <= 0:
        return 0.0, 0.0, g_obs - 1, n - g_obs
    ssb = float((sums[obs] ** 2 / counts[obs]).sum() - n * grand * grand)
    ssw = max(sst - ssb, 0.0)
    msw = ssw / (n - g_obs)
    omega2 = (ssb - (g_obs - 1) * msw) / (sst + msw)
    F = (ssb / (g_obs - 1)) / msw if msw > 0 else 0.0
    return max(float(omega2), 0.0), float(F), g_obs - 1, n - g_obs


def pair_index(cols):
    """cols: list of {name, kind, levels, data}. Returns pair list with
    branch cell keys and within-kind positions."""
    npos = {}
    cnt = {"num": 0, "cat": 0}
    for idx, c in enumerate(cols):
        npos[idx] = cnt[c["kind"]]
        cnt[c["kind"]] += 1
    pairs = []
    for i in range(len(cols)):
        for j in range(i + 1, len(cols)):
            a, b = cols[i], cols[j]
            if a["kind"] == "num" and b["kind"] == "num":
                cell = "num_num"
            elif a["kind"] == "cat" and b["kind"] == "cat":
                cell = {(3, 3): "cat3x3", (2, 2): "cat2x2",
                        (3, 2): "cat3x2", (2, 3): "cat3x2"}[
                            (a["levels"], b["levels"])]
            else:
                g = a["levels"] if a["kind"] == "cat" else b["levels"]
                cell = f"numcat_g{g}"
            pairs.append({"i": i, "j": j, "cell": cell,
                          "names": (a["name"], b["name"]),
                          "pi": npos[i], "pj": npos[j]})
    return pairs


def eval_pairs(cols, pairs, rows=None, need_p=False):
    """Corrected strengths for all pairs on `rows` (default all rows);
    p-values too if need_p. Returns (S, P or None)."""
    def sub(c):
        return c["data"] if rows is None else c["data"][rows]

    num_cols = [c for c in cols if c["kind"] == "num"]
    cat_cols = [c for c in cols if c["kind"] == "cat"]
    Xn = (np.column_stack([sub(c) for c in num_cols])
          if num_cols else np.zeros((0, 0)))
    m = Xn.shape[0] if num_cols else len(sub(cat_cols[0]))
    cat_data = [sub(c) for c in cat_cols]
    cat_lev = [c["levels"] for c in cat_cols]
    # numeric correlation matrix on the row subset
    Rmat = None
    if num_cols:
        Xc = Xn - Xn.mean(axis=0)
        norms = np.sqrt((Xc * Xc).sum(axis=0))
        norms = np.where(norms > 0, norms, 1.0)
        Z = Xc / norms
        Rmat = np.clip(Z.T @ Z, -1.0, 1.0)
    S = np.empty(len(pairs))
    chi2s, dfs, fs, d1s, d2s, rs = [], [], [], [], [], []
    idx_cc, idx_nc, idx_nn = [], [], []
    for t, pr in enumerate(pairs):
        a, b = cols[pr["i"]], cols[pr["j"]]
        if a["kind"] == "num" and b["kind"] == "num":
            r = float(Rmat[pr["pi"], pr["pj"]])
            S[t] = abs(r)
            rs.append(r); idx_nn.append(t)
        elif a["kind"] == "cat" and b["kind"] == "cat":
            vs, chi2, df = cat_strength_chi2(
                cat_data[pr["pi"]], cat_data[pr["pj"]],
                cat_lev[pr["pi"]], cat_lev[pr["pj"]])
            S[t] = vs
            chi2s.append(chi2); dfs.append(df); idx_cc.append(t)
        else:
            if a["kind"] == "num":
                x = Xn[:, pr["pi"]]; codes = cat_data[pr["pj"]]
                g = cat_lev[pr["pj"]]
            else:
                x = Xn[:, pr["pj"]]; codes = cat_data[pr["pi"]]
                g = cat_lev[pr["pi"]]
            s, F, d1, d2 = numcat_strength_F(x, codes, g)
            S[t] = s
            fs.append(F); d1s.append(d1); d2s.append(d2); idx_nc.append(t)
    if not need_p:
        return S, None
    P = np.ones(len(pairs))
    if idx_nn:
        r = np.asarray(rs)
        ok = np.abs(r) < 1.0
        tt = np.zeros_like(r)
        tt[ok] = (np.abs(r[ok])
                  * np.sqrt((m - 2) / np.maximum(1 - r[ok] ** 2, 1e-300)))
        P[idx_nn] = np.where(ok, 2.0 * st.t.sf(tt, m - 2), 0.0)
    if idx_cc:
        P[idx_cc] = st.chi2.sf(np.asarray(chi2s), np.asarray(dfs))
    if idx_nc:
        P[idx_nc] = st.f.sf(np.asarray(fs), np.asarray(d1s),
                            np.asarray(d2s))
    return S, P


def run_world(cols, pairs, split_h1, split_h2):
    """Full conjoint rule on one world. Returns counts and masks."""
    S_full, P_full = eval_pairs(cols, pairs, need_p=True)
    S1, _ = eval_pairs(cols, pairs, split_h1)
    S2, _ = eval_pairs(cols, pairs, split_h2)
    bh = bh_reject(P_full, Q_BH)
    floor = S_full >= TAU
    p3 = (S1 >= TAU) & (S2 >= TAU)
    promoted = bh & floor & p3
    hec = {}
    for t, pr in enumerate(pairs):
        e, v = hec.get(pr["cell"], (0, 0))
        hec[pr["cell"]] = (e + int(S1[t] >= TAU) + int(S2[t] >= TAU),
                           v + 2)
    return {"n_promoted": int(promoted.sum()),
            "n_bh_floor": int((bh & floor).sum()),
            "n_p3_only": int(p3.sum()),
            "promoted_mask": promoted, "p3_mask": p3,
            "half_exceed": int((S1 >= TAU).sum() + (S2 >= TAU).sum()),
            "half_evals": 2 * len(pairs),
            "half_exceed_by_cell": hec}


def certificate_for(pairs, total_n, gt, truth_mask=None):
    """Sum over (null) pairs of gammaU(cell, n) x gammaU(cell, n)
    (equal halves). Returns (total, pair-counts-by-cell)."""
    tot = 0.0
    comp = {}
    for t, pr in enumerate(pairs):
        if truth_mask is not None and truth_mask[t]:
            continue
        g = gt[(pr["cell"], total_n)]
        tot += g * g
        comp[pr["cell"]] = comp.get(pr["cell"], 0) + 1
    return tot, comp


def base40_columns(rng, n):
    """Probe-2 base design: 24 numeric, 8 cat-3, 8 cat-2, mutually
    independent, i.i.d. rows."""
    cols = []
    for i in range(24):
        cols.append({"name": f"num{i}", "kind": "num", "levels": None,
                     "data": rng.normal(0, 1, n)})
    for i in range(8):
        cols.append({"name": f"cat3_{i}", "kind": "cat", "levels": 3,
                     "data": rng.randint(0, 3, n)})
    for i in range(8):
        cols.append({"name": f"cat2_{i}", "kind": "cat", "levels": 2,
                     "data": rng.randint(0, 2, n)})
    return cols


# ----------------------------------------------------------------------
# Section 0 -- validation: facade_v2 vs reference stats; reference
# stats vs batched engine
# ----------------------------------------------------------------------
def section0_validation():
    import pandas as pd
    import facade_v2 as F

    rng = np.random.RandomState(81001)
    records = []
    max_ds, max_dp = 0.0, 0.0
    specs = []
    for i in range(10):
        specs.append(("nn", 40 + 13 * i))
    for i in range(10):
        specs.append((["cc33", "cc22", "cc32"][i % 3], 40 + 17 * i))
    for i in range(10):
        specs.append((["nc3", "nc2"][i % 2], 40 + 11 * i))
    for kind, n in specs:
        if kind == "nn":
            x = rng.normal(0, 1, n); y = rng.normal(0, 1, n)
            rec = F.correlate(pd.Series(x), pd.Series(y),
                              "numeric", "numeric")
            r = pearson_r_vec(x, y)
            fast_s, fast_p = abs(r), pearson_p_from_r(r, n)
            fac_s = abs(rec["strength"])
        elif kind.startswith("cc"):
            la, lb = {"cc33": (3, 3), "cc22": (2, 2),
                      "cc32": (3, 2)}[kind]
            a = rng.randint(0, la, n); b = rng.randint(0, lb, n)
            a[:la] = np.arange(la); b[:lb] = np.arange(lb)
            rec = F.correlate(pd.Series(a), pd.Series(b),
                              "categorical", "categorical")
            fast_s, fast_p = cat_pair_stats(a, b, la, lb)
            fac_s = rec["strength_corrected"]
        else:
            g = 3 if kind == "nc3" else 2
            x = rng.normal(0, 1, n); c = rng.randint(0, g, n)
            c[:g] = np.arange(g)
            rec = F.correlate(pd.Series(x), pd.Series(c),
                              "numeric", "categorical")
            om, fast_p = numcat_stats(x, c, g)
            fast_s = max(om, 0.0)
            fac_s = rec["strength_corrected"]
        fac_p = float(rec["p_value"])
        ds = abs(fast_s - fac_s); dp = abs(fast_p - fac_p)
        max_ds = max(max_ds, ds); max_dp = max(max_dp, dp)
        records.append({"kind": kind, "n": int(n),
                        "fast_strength": float(fast_s),
                        "facade_strength": float(fac_s),
                        "abs_diff_strength": float(ds),
                        "fast_p": float(fast_p), "facade_p": float(fac_p),
                        "abs_diff_p": float(dp)})
    # engine check: random world, batched engine vs reference per pair
    rng2 = np.random.RandomState(81002)
    cols = base40_columns(rng2, 120)
    pairs = pair_index(cols)
    S_eng, P_eng = eval_pairs(cols, pairs, need_p=True)
    max_es, max_ep = 0.0, 0.0
    for t, pr in enumerate(pairs):
        a, b = cols[pr["i"]], cols[pr["j"]]
        if a["kind"] == "num" and b["kind"] == "num":
            r = pearson_r_vec(a["data"], b["data"])
            s_ref, p_ref = abs(r), pearson_p_from_r(r, 120)
        elif a["kind"] == "cat" and b["kind"] == "cat":
            s_ref, p_ref = cat_pair_stats(a["data"], b["data"],
                                          a["levels"], b["levels"])
        else:
            if a["kind"] == "num":
                x, cd, g = a["data"], b["data"], b["levels"]
            else:
                x, cd, g = b["data"], a["data"], a["levels"]
            om, p_ref = numcat_stats(x, cd, g)
            s_ref = max(om, 0.0)
        max_es = max(max_es, abs(S_eng[t] - s_ref))
        max_ep = max(max_ep, abs(P_eng[t] - p_ref))
    out = {"facade_pairs": len(records), "seed": 81001,
           "max_abs_diff_strength_vs_facade": float(max_ds),
           "max_abs_diff_p_vs_facade": float(max_dp),
           "engine_world_pairs": len(pairs), "engine_seed": 81002,
           "engine_max_abs_diff_strength": float(max_es),
           "engine_max_abs_diff_p": float(max_ep),
           "tolerance_strength": 1e-8, "tolerance_p": 1e-6,
           "pass": bool(max_ds < 1e-8 and max_dp < 1e-6
                        and max_es < 1e-8 and max_ep < 1e-6),
           "facade_import_note": ("facade_v2 imported with np.math "
                                  "restored to the stdlib math module "
                                  "(numpy 2.5 removed the alias); "
                                  "facade_v2.py itself unmodified."),
           "records": records}
    with open(os.path.join(HERE, "probe8_validation.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"[validation] facade: {len(records)} pairs max|dS|={max_ds:.3e} "
          f"max|dp|={max_dp:.3e}; engine: max|dS|={max_es:.3e} "
          f"max|dp|={max_ep:.3e}; pass={out['pass']}", flush=True)
    assert out["pass"], "statistics disagree (facade / engine)"
    return out


# ----------------------------------------------------------------------
# Section A -- calibration of per-half exceedance gamma (vectorised)
# ----------------------------------------------------------------------
CAL_TOTALS = [60, 100, 200, 500, 1000]  # 200 added for the stress test
CAL_BRANCHES = ["num_num", "cat3x3", "cat2x2", "cat3x2",
                "numcat_g3", "numcat_g2"]
CAL_REPS = 50000
CAL_CHUNK = 4000


def cp_upper(k, n, alpha):
    """One-sided Clopper-Pearson upper bound at level 1-alpha."""
    if k >= n:
        return 1.0
    return float(st.beta.ppf(1.0 - alpha, k + 1, n - k))


def _cal_cat(rng, m, reps, ra, rb):
    k = 0
    done = 0
    rt = ra - ((ra - 1) ** 2) / (m - 1.0)
    ct_ = rb - ((rb - 1) ** 2) / (m - 1.0)
    denom = min(rt - 1.0, ct_ - 1.0)
    while done < reps:
        c = min(CAL_CHUNK, reps - done)
        A = rng.randint(0, ra, (c, m)); B = rng.randint(0, rb, (c, m))
        C = A * rb + B
        T = np.stack([(C == v).sum(axis=1)
                      for v in range(ra * rb)], axis=1).astype(float)
        T = T.reshape(c, ra, rb)
        rs = T.sum(axis=2); cs = T.sum(axis=1)
        zm = (rs == 0).any(axis=1) | (cs == 0).any(axis=1)
        E = rs[:, :, None] * cs[:, None, :] / m
        chi2 = (((T - E) ** 2)
                / np.where(E < 1e-12, 1e-12, E)).sum(axis=(1, 2))
        phi2_t = np.maximum(0.0, chi2 / m
                            - ((ra - 1) * (rb - 1)) / (m - 1.0))
        vs = np.sqrt(phi2_t / denom) if denom > 0 else np.zeros(c)
        k += int((vs[~zm] >= TAU).sum())
        for i in np.nonzero(zm)[0]:  # rare degenerate tables: slow path
            s, _ = cat_pair_stats(A[i], B[i], ra, rb)
            if s >= TAU:
                k += 1
        done += c
    return k


def _cal_numcat(rng, m, reps, g):
    k = 0
    done = 0
    while done < reps:
        c = min(CAL_CHUNK, reps - done)
        X = rng.normal(0, 1, (c, m)); C = rng.randint(0, g, (c, m))
        counts = np.stack([(C == v).sum(axis=1)
                           for v in range(g)], axis=1).astype(float)
        sums = np.stack([(X * (C == v)).sum(axis=1)
                         for v in range(g)], axis=1)
        obs = counts > 0
        gobs = obs.sum(axis=1).astype(float)
        grand = X.mean(axis=1)
        sst = ((X - grand[:, None]) ** 2).sum(axis=1)
        ssb = (np.where(obs, sums ** 2 / np.maximum(counts, 1.0), 0.0)
               ).sum(axis=1) - m * grand ** 2
        ssw = np.maximum(sst - ssb, 0.0)
        msw = ssw / np.maximum(m - gobs, 1.0)
        omega = (ssb - (gobs - 1.0) * msw) / np.maximum(sst + msw, 1e-300)
        strength = np.where(gobs >= 2, np.maximum(omega, 0.0), 0.0)
        k += int((strength >= TAU).sum())
        done += c
    return k


def calibrate_branch(branch, m, reps, seed):
    rng = np.random.RandomState(seed)
    if branch == "num_num":
        k = 0
        done = 0
        while done < reps:
            c = min(CAL_CHUNK, reps - done)
            X = rng.normal(0, 1, (c, m)); Y = rng.normal(0, 1, (c, m))
            Xc = X - X.mean(axis=1, keepdims=True)
            Yc = Y - Y.mean(axis=1, keepdims=True)
            r = (Xc * Yc).sum(axis=1) / np.sqrt(
                (Xc * Xc).sum(axis=1) * (Yc * Yc).sum(axis=1))
            k += int((np.abs(r) >= TAU).sum())
            done += c
        return k
    if branch.startswith("cat"):
        ra, rb = {"cat3x3": (3, 3), "cat2x2": (2, 2),
                  "cat3x2": (3, 2)}[branch]
        return _cal_cat(rng, m, reps, ra, rb)
    g = 3 if branch == "numcat_g3" else 2
    return _cal_numcat(rng, m, reps, g)


def section_a_calibration():
    t0 = time.time()
    n_cells = len(CAL_TOTALS) * len(CAL_BRANCHES)
    alpha_sim = DELTA / n_cells  # Bonferroni across all calibration cells
    cells = []
    idx = 0
    for total_n in CAL_TOTALS:
        m1, m2 = total_n // 2, total_n - total_n // 2
        assert m1 == m2
        for branch in CAL_BRANCHES:
            seed = 810000 + idx * 7919
            k = calibrate_branch(branch, m1, CAL_REPS, seed)
            gh = k / CAL_REPS
            cells.append({
                "branch": branch, "total_n": total_n,
                "half_n1": m1, "half_n2": m2, "reps": CAL_REPS,
                "exceedances": int(k), "gamma_hat": float(gh),
                "cp95_upper": cp_upper(k, CAL_REPS, 0.05),
                "cp_simultaneous_upper": cp_upper(k, CAL_REPS, alpha_sim),
                "alpha_simultaneous_per_cell": alpha_sim,
                "seed": seed})
            idx += 1
            print(f"[cal] {branch} total_n={total_n} half={m1}: "
                  f"k={k}/{CAL_REPS} gamma_hat={gh:.6g} "
                  f"U_sim={cells[-1]['cp_simultaneous_upper']:.6g}",
                  flush=True)
    out = {"tau": TAU, "reps_per_cell": CAL_REPS,
           "cells_total": n_cells, "delta_simultaneous": DELTA,
           "note": ("gamma is the per-half exceedance probability "
                    "Pr[S_h >= tau] under the branch null at the half "
                    "sample size; totals are even so half_n1 == half_n2 "
                    "and one estimate serves both halves. num_num: "
                    "independent Gaussians; cat cells: balanced uniform "
                    "levels; numcat: Gaussian numeric vs balanced "
                    "uniform groups. Certificates use "
                    "cp_simultaneous_upper (Bonferroni over all "
                    "30 cells, simultaneous level 0.95)."),
           "cells": cells}
    with open(os.path.join(HERE, "probe8_calibration.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    stamp("calibration", t0)
    return out


def gamma_tables(cal):
    return {(c["branch"], c["total_n"]): c["cp_simultaneous_upper"]
            for c in cal["cells"]}


def gamma_hat_tables(cal):
    return {(c["branch"], c["total_n"]): c["gamma_hat"]
            for c in cal["cells"]}


# ----------------------------------------------------------------------
# Section B -- all-null world validation
# ----------------------------------------------------------------------
def section_b_null_worlds(gt, worlds=50, ghat=None):
    t0 = time.time()
    res = {}
    for n in (60, 100, 500):
        counts, p3counts, fwers = [], [], []
        cert = None
        cert_pt = None
        comp = None
        for w in range(worlds):
            rng = np.random.RandomState(820000 + n * 100 + w)
            cols = base40_columns(rng, n)
            pairs = pair_index(cols)
            if cert is None:
                cert, comp = certificate_for(pairs, n, gt)
                if ghat is not None:
                    cert_pt, _ = certificate_for(pairs, n, ghat)
            srng = np.random.RandomState(830000 + n * 100 + w)
            perm = srng.permutation(n)
            out = run_world(cols, pairs, perm[:n // 2], perm[n // 2:])
            counts.append(out["n_promoted"])
            p3counts.append(out["n_p3_only"])
            fwers.append(1 if out["n_promoted"] > 0 else 0)
        mean_v = float(np.mean(counts))
        res[str(n)] = {
            "n": n, "worlds": worlds, "pairs_per_world": 780,
            "data_seed_base": 820000 + n * 100,
            "split_seed_base": 830000 + n * 100,
            "per_world_false_discoveries": counts,
            "mean_false_discoveries": mean_v,
            "mean_false_se": float(np.std(counts, ddof=1)
                                   / math.sqrt(worlds)),
            "fwer": float(np.mean(fwers)),
            "fdr_equals_fwer_in_all_null": True,
            "per_world_p3_only_counts": p3counts,
            "mean_p3_only_count": float(np.mean(p3counts)),
            "p3_mean_se": float(np.std(p3counts, ddof=1)
                                / math.sqrt(worlds)),
            "certificate_point_estimate_EV":
                (float(cert_pt) if cert_pt is not None else None),
            "certificate_EV_upper": float(cert),
            "certificate_FWER_FDR_upper": float(min(1.0, cert)),
            "certificate_pair_counts_by_cell": comp,
            "mean_le_certificate": bool(mean_v <= cert),
            "p3_mean_le_certificate":
                bool(float(np.mean(p3counts)) <= cert),
            "conservativeness_ratio_cert_over_mean":
                (float(cert / mean_v) if mean_v > 0 else None),
            "excess_certificate_minus_mean": float(cert - mean_v)}
        print(f"[null n={n}] meanV={mean_v:.4f} "
              f"FWER={np.mean(fwers):.3f} "
              f"P3mean={np.mean(p3counts):.3f} cert={cert:.4f}",
              flush=True)
    with open(os.path.join(HERE, "probe8_worlds_null.json"), "w") as fh:
        json.dump(res, fh, indent=2)
    stamp("null worlds", t0)
    return res


# ----------------------------------------------------------------------
# Section C -- planted worlds (probe-2 planted design)
# ----------------------------------------------------------------------
def _popv(p):
    off = (1 - p) / 2
    Pc = np.array([[p if i == j else off for j in range(3)]
                   for i in range(3)])
    J = Pc / 3.0
    E = np.full((3, 3), 1 / 9)
    return float(np.sqrt(float(((J - E) ** 2 / E).sum()) / 2))


P45 = brentq(lambda q: _popv(q) - 0.45, 0.34, 0.999)
P25 = brentq(lambda q: _popv(q) - 0.25, 0.34, 0.999)


def _plant_v_codes(rng, n, p):
    a = rng.randint(0, 3, n)
    off = (1 - p) / 2
    u = rng.rand(n)
    b = np.where(u < p, a, np.where(u < p + off, (a + 1) % 3, (a + 2) % 3))
    return a, b


def planted_columns(rng, n):
    """Probe-2 planted design, transcribed from probe2_fdr.py
    planted_columns(): the 40 base columns plus 12 plant columns = 52
    columns, C(52,2) = 1326 pairs, 6 planted pairs. (probe2's docstring
    and the manuscript say 46 columns / 1035 pairs; the code as written
    builds 52 / 1326. This probe follows the code and flags the gap.)"""
    cols = base40_columns(rng, n)
    x = rng.normal(0, 1, n)
    cols.append({"name": "pnum_a", "kind": "num", "levels": None,
                 "data": x})
    cols.append({"name": "pnum_b", "kind": "num", "levels": None,
                 "data": 0.5 * x + math.sqrt(0.75) * rng.normal(0, 1, n)})
    x2 = rng.normal(0, 1, n)
    cols.append({"name": "pnum_c", "kind": "num", "levels": None,
                 "data": x2})
    cols.append({"name": "pnum_d", "kind": "num", "levels": None,
                 "data": 0.35 * x2 + math.sqrt(1 - 0.35 ** 2)
                 * rng.normal(0, 1, n)})
    g = rng.randint(0, 3, n)
    m25 = np.array([-math.sqrt(0.5), 0.0, math.sqrt(0.5)])
    cols.append({"name": "pcat_a", "kind": "cat", "levels": 3, "data": g})
    cols.append({"name": "pnum_e", "kind": "num", "levels": None,
                 "data": m25[g] + rng.normal(0, 1, n)})
    g2 = rng.randint(0, 3, n)
    m10 = np.array([-math.sqrt(1 / 6), 0.0, math.sqrt(1 / 6)])
    cols.append({"name": "pcat_b", "kind": "cat", "levels": 3, "data": g2})
    cols.append({"name": "pnum_f", "kind": "num", "levels": None,
                 "data": m10[g2] + rng.normal(0, 1, n)})
    a1, b1 = _plant_v_codes(rng, n, P45)
    cols.append({"name": "pcat_c", "kind": "cat", "levels": 3, "data": a1})
    cols.append({"name": "pcat_d", "kind": "cat", "levels": 3, "data": b1})
    a2, b2 = _plant_v_codes(rng, n, P25)
    cols.append({"name": "pcat_e", "kind": "cat", "levels": 3, "data": a2})
    cols.append({"name": "pcat_f", "kind": "cat", "levels": 3, "data": b2})
    truth = {frozenset(("pnum_a", "pnum_b")),
             frozenset(("pnum_c", "pnum_d")),
             frozenset(("pcat_a", "pnum_e")),
             frozenset(("pcat_b", "pnum_f")),
             frozenset(("pcat_c", "pcat_d")),
             frozenset(("pcat_e", "pcat_f"))}
    return cols, truth


def section_c_planted(gt, worlds=20):
    t0 = time.time()
    n = 500
    per_world = []
    cert = None
    null_cell_counts = None
    plant_hits = {}
    for w in range(worlds):
        rng = np.random.RandomState(840000 + w)
        cols, truth = planted_columns(rng, n)
        pairs = pair_index(cols)
        truth_mask = np.array([frozenset(pr["names"]) in truth
                               for pr in pairs])
        if cert is None:
            cert, null_cell_counts = certificate_for(pairs, n, gt,
                                                     truth_mask)
        srng = np.random.RandomState(850000 + w)
        perm = srng.permutation(n)
        out = run_world(cols, pairs, perm[:n // 2], perm[n // 2:])
        prom = out["promoted_mask"]
        n_false = int((prom & ~truth_mask).sum())
        n_found = int((prom & truth_mask).sum())
        for t in np.nonzero(prom & truth_mask)[0]:
            key = " & ".join(sorted(pairs[t]["names"]))
            plant_hits[key] = plant_hits.get(key, 0) + 1
        per_world.append({"world": w, "promoted": int(prom.sum()),
                          "false": n_false, "plants_found": n_found,
                          "fdp": float(n_false / max(int(prom.sum()), 1)),
                          "power": float(n_found / 6.0)})
    res = {"n": n, "worlds": worlds, "columns": 52, "pairs": 1326,
           "null_pairs": 1320, "plants": 6,
           "design_note": ("probe2_fdr.py planted_columns() as written "
                           "builds 52 columns / 1326 pairs, not the 46 / "
                           "1035 of its docstring and the manuscript "
                           "table; this probe follows the code."),
           "plant_conditional_probs": {"p_for_V0.45": float(P45),
                                       "p_for_V0.25": float(P25)},
           "data_seed_base": 840000, "split_seed_base": 850000,
           "per_world": per_world,
           "mean_fdp": float(np.mean([r["fdp"] for r in per_world])),
           "mean_power": float(np.mean([r["power"] for r in per_world])),
           "mean_false": float(np.mean([r["false"] for r in per_world])),
           "plant_detection_counts": plant_hits,
           "certificate_EV_upper_null_pairs_only": float(cert),
           "certificate_null_pair_counts_by_cell": null_cell_counts,
           "mean_false_le_certificate":
               bool(np.mean([r["false"] for r in per_world]) <= cert)}
    with open(os.path.join(HERE, "probe8_worlds_planted.json"), "w") as fh:
        json.dump(res, fh, indent=2)
    print(f"[planted] meanFDP={res['mean_fdp']:.4f} "
          f"power={res['mean_power']:.3f} "
          f"meanFalse={res['mean_false']:.3f} cert_null={cert:.4f}",
          flush=True)
    stamp("planted worlds", t0)
    return res


# ----------------------------------------------------------------------
# Section D -- cluster stress test
# ----------------------------------------------------------------------
def clustered_columns(rng, n, k_clusters, icc_var=1.0, dir_alpha=3.0):
    """40 columns, mutually independent across columns, but rows within
    a column share a cluster random effect (clusters of n/k_clusters).
    Numeric: X = a_cluster + e, a ~ N(0, icc_var), e ~ N(0, 1)
    (ICC = icc_var / (icc_var + 1) = 0.5 at the default).
    Categorical: per-cluster probability vector ~ Dirichlet(alpha * p0)
    with p0 balanced; rows drawn from their cluster's vector."""
    size = n // k_clusters
    assert size * k_clusters == n
    cl = np.repeat(np.arange(k_clusters), size)
    cols = []
    for i in range(24):
        a = rng.normal(0, math.sqrt(icc_var), k_clusters)
        cols.append({"name": f"num{i}", "kind": "num", "levels": None,
                     "data": a[cl] + rng.normal(0, 1, n)})
    specs = [(f"cat3_{i}", 3) for i in range(8)] + \
            [(f"cat2_{i}", 2) for i in range(8)]
    for name, lev in specs:
        p0 = np.full(lev, 1.0 / lev)
        probs = np.array([rng.dirichlet(dir_alpha * p0)
                          for _ in range(k_clusters)])
        u = rng.rand(n)
        codes = np.empty(n, dtype=int)
        for kk in range(k_clusters):
            idx = cl == kk
            codes[idx] = np.searchsorted(np.cumsum(probs[kk]), u[idx])
            codes[idx] = np.minimum(codes[idx], lev - 1)
        cols.append({"name": name, "kind": "cat", "levels": lev,
                     "data": codes})
    return cols, cl


def section_d_cluster(gt, worlds=50):
    t0 = time.time()
    n, k_clusters = 200, 20
    cert, comp = None, None
    agg = {"random_row_split": {"counts": [], "p3": [], "hexceed": 0,
                                "hevals": 0, "hec": {}},
           "cluster_level_split": {"counts": [], "p3": [], "hexceed": 0,
                                   "hevals": 0, "hec": {}}}
    for w in range(worlds):
        rng = np.random.RandomState(860000 + w)
        cols, cl = clustered_columns(rng, n, k_clusters)
        pairs = pair_index(cols)
        if cert is None:
            cert, comp = certificate_for(pairs, n, gt)
        srng = np.random.RandomState(870000 + w)
        perm = srng.permutation(n)
        out_r = run_world(cols, pairs, perm[:n // 2], perm[n // 2:])
        crng = np.random.RandomState(880000 + w)
        cperm = crng.permutation(k_clusters)
        h1 = np.isin(cl, cperm[:k_clusters // 2])
        out_c = run_world(cols, pairs, np.nonzero(h1)[0],
                          np.nonzero(~h1)[0])
        for key, out in (("random_row_split", out_r),
                         ("cluster_level_split", out_c)):
            a = agg[key]
            a["counts"].append(out["n_promoted"])
            a["p3"].append(out["n_p3_only"])
            a["hexceed"] += out["half_exceed"]
            a["hevals"] += out["half_evals"]
            for cell, (e, v) in out["half_exceed_by_cell"].items():
                pe, pv = a["hec"].get(cell, (0, 0))
                a["hec"][cell] = (pe + e, pv + v)
    res = {"n": n, "clusters": k_clusters,
           "cluster_size": n // k_clusters,
           "worlds": worlds, "pairs_per_world": 780,
           "numeric_cluster_model":
               "X = a_c + e, a~N(0,1), e~N(0,1), ICC=0.5",
           "categorical_cluster_model":
               "per-cluster p ~ Dirichlet(3 * balanced), rows ~ p_cluster",
           "data_seed_base": 860000, "rowsplit_seed_base": 870000,
           "clustersplit_seed_base": 880000,
           "certificate_EV_upper_iid": float(cert),
           "certificate_pair_counts_by_cell": comp}
    pair_counts = comp
    for key, a in agg.items():
        mean_v = float(np.mean(a["counts"]))
        hec_out = {}
        indep_pred = 0.0
        for cell, (e, v) in sorted(a["hec"].items()):
            rate = e / v
            hec_out[cell] = {"exceedances": int(e), "evals": int(v),
                             "rate": float(rate)}
            indep_pred += pair_counts[cell] * rate * rate
        res[key] = {
            "per_world_false": a["counts"],
            "mean_false": mean_v,
            "mean_false_se": float(np.std(a["counts"], ddof=1)
                                   / math.sqrt(worlds)),
            "fwer": float(np.mean([1 if c > 0 else 0
                                    for c in a["counts"]])),
            "per_world_p3_only": a["p3"],
            "mean_p3_only": float(np.mean(a["p3"])),
            "p3_mean_se": float(np.std(a["p3"], ddof=1)
                                / math.sqrt(worlds)),
            "empirical_per_half_exceedance":
                float(a["hexceed"] / a["hevals"]),
            "half_exceedance_by_cell": hec_out,
            "independence_prediction_P3_from_clustered_marginals":
                float(indep_pred),
            "dependence_factor_p3_over_independence_prediction":
                float(np.mean(a["p3"]) / indep_pred),
            "mean_le_certificate": bool(mean_v <= cert),
            "p3_mean_le_certificate":
                bool(float(np.mean(a["p3"])) <= cert),
            "certificate_violated_by_mean": bool(mean_v > cert),
            "certificate_violated_by_p3_mean":
                bool(float(np.mean(a["p3"])) > cert)}
        print(f"[cluster {key}] meanV={mean_v:.3f} "
              f"P3mean={np.mean(a['p3']):.3f} "
              f"halfExc={a['hexceed'] / a['hevals']:.4f} "
              f"cert={cert:.4f}", flush=True)
    with open(os.path.join(HERE, "probe8_cluster.json"), "w") as fh:
        json.dump(res, fh, indent=2)
    stamp("cluster stress", t0)
    return res


def main():
    t_all = time.time()
    section0_validation()
    cal = section_a_calibration()
    gt = gamma_tables(cal)
    section_b_null_worlds(gt, worlds=50)
    section_c_planted(gt, worlds=20)
    section_d_cluster(gt, worlds=50)
    print(f"ALL DONE in {time.time() - t_all:.1f}s", flush=True)


if __name__ == "__main__":
    main()
