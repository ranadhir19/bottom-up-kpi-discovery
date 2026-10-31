"""Probe 10 -- cluster-level family error certificate for Paper 8, cycle 2.

Design, seeds and decision rules are pre-registered in
cluster_preregistration.md (written before this script was run).

Sections:
  0. Statistic validation: cycle-1 scalar references vs facade_v2 on 30
     clustered null pairs; probe-10 batched calibration kernels vs the
     scalar references on 2,000 clustered draws per branch.
  A. Calibration of per-half exceedance gamma^cl under the clustered
     sampling, per branch, for two designs:
       P (primary):   halves of 10 clusters x 10 rows (from 20 x 10 worlds)
       S (sensitivity): halves of 20 clusters x 5 rows (from 40 x 5 worlds)
     Numeric ICC 0.5 (a ~ N(0,1), e ~ N(0,1)); categorical per-cluster
     probabilities ~ Dirichlet(3 * balanced) -- exactly the cycle-1
     clustered design (probe8_certificate.clustered_columns).
     50,000 independent null pairs per cell; Clopper-Pearson bounds,
     simultaneous over each design's 6 cells (delta = 0.05, Bonferroni).
  B. Fresh all-null clustered worlds, design P (60 worlds): conjoint rule
     and P3-only counts under the cluster-level split (certified) and,
     descriptively, under a random row split (not covered).
  C. Fresh all-null clustered worlds, design S (50 worlds), cluster split.

Worlds are generated with probe8_certificate.clustered_columns imported
unmodified; the world engine is probe8_certificate.eval_pairs/run_world,
imported unmodified. Run with the archived venv python:
  hidden_files/paper08_work/phik_venv/bin/python probe10_cluster_certificate.py
All outputs are JSON files written next to this script.
"""
import json
import math
import os
import sys
import time

import numpy as np

if not hasattr(np, "math"):
    np.math = math

HERE = os.path.dirname(os.path.abspath(__file__))
THEOREM_DIR = os.path.join(HERE, "..", "theorem")
sys.path.insert(0, os.path.normpath(THEOREM_DIR))

import probe8_certificate as P8  # noqa: E402  (also puts probes/ on sys.path)

TAU = P8.TAU           # 0.20
DELTA = 0.05
BRANCHES = ["num_num", "cat3x3", "cat2x2", "cat3x2",
            "numcat_g3", "numcat_g2"]
CAT_LEVELS = {"cat3x3": (3, 3), "cat2x2": (2, 2), "cat3x2": (3, 2)}
CAL_REPS = 50000
CAL_CHUNK = 2500
ICC_VAR = 1.0          # numeric cluster-effect variance -> ICC 0.5
DIR_ALPHA = 3.0        # Dirichlet total concentration (cycle-1 design)

DESIGNS = {
    "P_primary_20x10": {"K_total": 20, "size": 10, "K_half": 10,
                        "cal_seed_base": 910000},
    "S_sens_40x5": {"K_total": 40, "size": 5, "K_half": 20,
                    "cal_seed_base": 920000},
}

T_START = time.time()


def stamp(section, t0):
    print(f"[{section}] elapsed {time.time() - t0:.1f}s "
          f"(total {time.time() - T_START:.1f}s)", flush=True)


# ----------------------------------------------------------------------
# Clustered column generators (batched; distributional design identical
# to P8.clustered_columns: numeric a_c + e with a ~ N(0, ICC_VAR);
# categorical per-cluster p ~ Dirichlet(DIR_ALPHA * balanced))
# ----------------------------------------------------------------------
def gen_num_clustered(rng, reps, K, s):
    a = rng.normal(0.0, math.sqrt(ICC_VAR), (reps, K, 1))
    e = rng.normal(0.0, 1.0, (reps, K, s))
    return (a + e).reshape(reps, K * s)


def gen_cat_clustered(rng, reps, K, s, lev):
    alpha = DIR_ALPHA * np.full(lev, 1.0 / lev)
    probs = rng.dirichlet(alpha, size=(reps, K))          # (reps, K, lev)
    u = rng.rand(reps, K, s)
    cum_excl = np.cumsum(probs, axis=2)[:, :, :-1]        # (reps, K, lev-1)
    codes = (u[:, :, :, None] > cum_excl[:, :, None, :]).sum(axis=3)
    return codes.reshape(reps, K * s).astype(int)


def gen_pair_batch(rng, branch, reps, K, s):
    """Two mutually independent clustered columns of the branch type."""
    if branch == "num_num":
        return (gen_num_clustered(rng, reps, K, s),
                gen_num_clustered(rng, reps, K, s))
    if branch in CAT_LEVELS:
        ra, rb = CAT_LEVELS[branch]
        return (gen_cat_clustered(rng, reps, K, s, ra),
                gen_cat_clustered(rng, reps, K, s, rb))
    g = 3 if branch == "numcat_g3" else 2
    return (gen_num_clustered(rng, reps, K, s),
            gen_cat_clustered(rng, reps, K, s, g))


# ----------------------------------------------------------------------
# Batched strength kernels (arithmetic copied from cycle 1:
# P8._cal_cat / P8._cal_numcat computation blocks; Pearson as in
# P8's num_num calibration)
# ----------------------------------------------------------------------
def strengths_num_num(X, Y):
    Xc = X - X.mean(axis=1, keepdims=True)
    Yc = Y - Y.mean(axis=1, keepdims=True)
    r = (Xc * Yc).sum(axis=1) / np.sqrt(
        (Xc * Xc).sum(axis=1) * (Yc * Yc).sum(axis=1))
    return np.abs(r)


def strengths_cat(A, B, ra, rb):
    """Bergsma V* batch; degenerate tables fall back to P8.cat_pair_stats
    (the cycle-1 slow path), exactly as P8._cal_cat does."""
    c, m = A.shape
    rt = ra - ((ra - 1) ** 2) / (m - 1.0)
    ct_ = rb - ((rb - 1) ** 2) / (m - 1.0)
    denom = min(rt - 1.0, ct_ - 1.0)
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
    vs = np.where(zm, 0.0, vs)
    for i in np.nonzero(zm)[0]:
        s, _ = P8.cat_pair_stats(A[i], B[i], ra, rb)
        vs[i] = s
    return vs


def strengths_numcat(X, C, g):
    c, m = X.shape
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
    return np.where(gobs >= 2, np.maximum(omega, 0.0), 0.0)


def kernel_strengths(branch, A, B):
    if branch == "num_num":
        return strengths_num_num(A, B)
    if branch in CAT_LEVELS:
        ra, rb = CAT_LEVELS[branch]
        return strengths_cat(A, B, ra, rb)
    g = 3 if branch == "numcat_g3" else 2
    return strengths_numcat(A, B, g)


def scalar_strength(branch, a, b):
    if branch == "num_num":
        return abs(P8.pearson_r_vec(a, b))
    if branch in CAT_LEVELS:
        ra, rb = CAT_LEVELS[branch]
        return P8.cat_pair_stats(a, b, ra, rb)[0]
    g = 3 if branch == "numcat_g3" else 2
    return max(P8.numcat_stats(a, b, g)[0], 0.0)


# ----------------------------------------------------------------------
# Section 0 -- validation
# ----------------------------------------------------------------------
def _single_clustered_column(rng, kind, levels, K, s):
    """One clustered column via the scalar path (same design)."""
    cl = np.repeat(np.arange(K), s)
    if kind == "num":
        a = rng.normal(0.0, math.sqrt(ICC_VAR), K)
        return a[cl] + rng.normal(0.0, 1.0, K * s)
    p0 = np.full(levels, 1.0 / levels)
    probs = np.array([rng.dirichlet(DIR_ALPHA * p0) for _ in range(K)])
    u = rng.rand(K * s)
    codes = np.empty(K * s, dtype=int)
    for kk in range(K):
        idx = cl == kk
        codes[idx] = np.searchsorted(np.cumsum(probs[kk]), u[idx])
        codes[idx] = np.minimum(codes[idx], levels - 1)
    # validation arrays only: guarantee every level appears so the
    # facade's crosstab path is exercised non-degenerately
    present = set(np.unique(codes).tolist())
    for lev in range(levels):
        if lev not in present:
            codes[lev * s] = lev
    return codes


def section0_validation():
    import pandas as pd
    import facade_v2 as F

    t0 = time.time()
    K, s = 10, 10                      # primary half size: 100 rows
    rng = np.random.RandomState(91001)
    records = []
    max_ds, max_dp = 0.0, 0.0
    specs = ([("nn", None)] * 10
             + [("cc33", (3, 3)), ("cc22", (2, 2)), ("cc32", (3, 2))] * 3
             + [("cc33", (3, 3))]
             + [("nc3", 3), ("nc2", 2)] * 5)
    assert len(specs) == 30
    for kind, par in specs:
        if kind == "nn":
            x = _single_clustered_column(rng, "num", None, K, s)
            y = _single_clustered_column(rng, "num", None, K, s)
            rec = F.correlate(pd.Series(x), pd.Series(y),
                              "numeric", "numeric")
            r = P8.pearson_r_vec(x, y)
            fast_s, fast_p = abs(r), P8.pearson_p_from_r(r, K * s)
            fac_s = abs(rec["strength"])
        elif kind.startswith("cc"):
            la, lb = par
            a = _single_clustered_column(rng, "cat", la, K, s)
            b = _single_clustered_column(rng, "cat", lb, K, s)
            rec = F.correlate(pd.Series(a), pd.Series(b),
                              "categorical", "categorical")
            fast_s, fast_p = P8.cat_pair_stats(a, b, la, lb)
            fac_s = rec["strength_corrected"]
        else:
            g = par
            x = _single_clustered_column(rng, "num", None, K, s)
            c = _single_clustered_column(rng, "cat", g, K, s)
            rec = F.correlate(pd.Series(x), pd.Series(c),
                              "numeric", "categorical")
            om, fast_p = P8.numcat_stats(x, c, g)
            fast_s = max(om, 0.0)
            fac_s = rec["strength_corrected"]
        fac_p = float(rec["p_value"])
        ds = abs(fast_s - fac_s); dp = abs(fast_p - fac_p)
        max_ds = max(max_ds, ds); max_dp = max(max_dp, dp)
        records.append({"kind": kind, "n": K * s,
                        "fast_strength": float(fast_s),
                        "facade_strength": float(fac_s),
                        "abs_diff_strength": float(ds),
                        "abs_diff_p": float(dp)})
    # (b) batched kernels vs scalar references, 2,000 clustered draws/cell
    kernel_checks = {}
    max_kd = 0.0
    all_counts_agree = True
    for bi, branch in enumerate(BRANCHES):
        rngk = np.random.RandomState(91002 + bi * 131)
        A, B = gen_pair_batch(rngk, branch, 2000, K, s)
        ks = kernel_strengths(branch, A, B)
        ss = np.array([scalar_strength(branch, A[i], B[i])
                       for i in range(2000)])
        md = float(np.max(np.abs(ks - ss)))
        agree = bool(int((ks >= TAU).sum()) == int((ss >= TAU).sum()))
        max_kd = max(max_kd, md)
        all_counts_agree = all_counts_agree and agree
        kernel_checks[branch] = {
            "draws": 2000, "max_abs_diff_strength": md,
            "kernel_exceedances": int((ks >= TAU).sum()),
            "scalar_exceedances": int((ss >= TAU).sum()),
            "exceedance_counts_agree": agree}
    out = {"facade_pairs": len(records), "seed": 91001,
           "design": "clustered half samples, 10 clusters x 10 rows, "
                     "ICC 0.5 numeric / Dirichlet(3*balanced) categorical",
           "max_abs_diff_strength_vs_facade": float(max_ds),
           "max_abs_diff_p_vs_facade": float(max_dp),
           "kernel_vs_scalar": kernel_checks,
           "kernel_max_abs_diff_strength": float(max_kd),
           "kernel_exceedance_counts_all_agree": bool(all_counts_agree),
           "tolerance_strength": 1e-8, "tolerance_p": 1e-6,
           "pass": bool(max_ds < 1e-8 and max_dp < 1e-6
                        and max_kd < 1e-8 and all_counts_agree),
           "records": records}
    with open(os.path.join(HERE, "probe10_validation.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"[validation] facade 30 clustered pairs: max|dS|={max_ds:.3e} "
          f"max|dp|={max_dp:.3e}; kernels vs scalar: max|dS|={max_kd:.3e} "
          f"counts_agree={all_counts_agree}; pass={out['pass']}",
          flush=True)
    assert out["pass"], "statistic validation failed"
    stamp("validation", t0)
    return out


# ----------------------------------------------------------------------
# Section A -- calibration of gamma^cl
# ----------------------------------------------------------------------
def calibrate_cell(branch, K, s, reps, seed):
    rng = np.random.RandomState(seed)
    k = 0
    done = 0
    while done < reps:
        c = min(CAL_CHUNK, reps - done)
        A, B = gen_pair_batch(rng, branch, c, K, s)
        k += int((kernel_strengths(branch, A, B) >= TAU).sum())
        done += c
    return k


def pair_counts():
    rng = np.random.RandomState(91003)
    cols = P8.base40_columns(rng, 200)
    pairs = P8.pair_index(cols)
    counts = {}
    for pr in pairs:
        counts[pr["cell"]] = counts.get(pr["cell"], 0) + 1
    return counts


def certificate_from_cells(cells, counts, key):
    return float(sum(counts[c["branch"]] * c[key] ** 2 for c in cells))


def delta_se_point(cells, counts, reps):
    var = 0.0
    for c in cells:
        g = c["gamma_hat_cl"]
        var += (counts[c["branch"]] ** 2) * (2.0 * g) ** 2 \
            * g * (1.0 - g) / reps
    return float(math.sqrt(var))


def section_a_calibration(counts):
    t0 = time.time()
    out = {"tau": TAU, "reps_per_cell": CAL_REPS, "delta": DELTA,
           "numeric_model": "X = a_c + e, a~N(0,1), e~N(0,1), ICC=0.5",
           "categorical_model": "per-cluster p ~ Dirichlet(3 * balanced), "
                                "rows ~ p_cluster",
           "note": ("gamma^cl is the per-half exceedance probability "
                    "Pr[S_h >= tau] under the clustered null at the "
                    "half's cluster count. Simultaneous bounds are "
                    "Clopper-Pearson at per-cell level 1 - delta/6 "
                    "(Bonferroni over the design's 6 cells); each design "
                    "is a separate 0.95-simultaneous certificate."),
           "pair_counts_by_cell": counts, "designs": {}}
    alpha_sim = DELTA / len(BRANCHES)
    for dname, d in DESIGNS.items():
        cells = []
        for bi, branch in enumerate(BRANCHES):
            seed = d["cal_seed_base"] + bi * 7919
            k = calibrate_cell(branch, d["K_half"], d["size"],
                               CAL_REPS, seed)
            gh = k / CAL_REPS
            cell = {"branch": branch, "K_half": d["K_half"],
                    "cluster_size": d["size"],
                    "half_rows": d["K_half"] * d["size"],
                    "reps": CAL_REPS, "exceedances": int(k),
                    "gamma_hat_cl": float(gh),
                    "cp95_upper": P8.cp_upper(k, CAL_REPS, 0.05),
                    "cp_simultaneous_upper":
                        P8.cp_upper(k, CAL_REPS, alpha_sim),
                    "alpha_simultaneous_per_cell": alpha_sim,
                    "seed": seed}
            cells.append(cell)
            print(f"[cal {dname}] {branch} K_half={d['K_half']}: "
                  f"k={k}/{CAL_REPS} gamma_cl={gh:.6g} "
                  f"U_sim={cell['cp_simultaneous_upper']:.6g}", flush=True)
        out["designs"][dname] = {
            "K_total": d["K_total"], "cluster_size": d["size"],
            "K_half": d["K_half"], "cells": cells,
            "certificate_point": certificate_from_cells(
                cells, counts, "gamma_hat_cl"),
            "certificate_simultaneous_upper": certificate_from_cells(
                cells, counts, "cp_simultaneous_upper"),
            "certificate_point_delta_se": delta_se_point(
                cells, counts, CAL_REPS)}
    with open(os.path.join(HERE, "probe10_calibration_cluster.json"),
              "w") as fh:
        json.dump(out, fh, indent=2)
    stamp("calibration", t0)
    return out


def iid_contrast(counts):
    """Cycle-1 i.i.d. certificate at total n = 200, recomputed from
    probe8_calibration.json with the same pair counts."""
    path = os.path.join(HERE, "..", "theorem", "probe8_calibration.json")
    with open(os.path.normpath(path)) as fh:
        cal8 = json.load(fh)
    pt = 0.0
    up = 0.0
    detail = {}
    for c in cal8["cells"]:
        if c["total_n"] != 200:
            continue
        n = counts[c["branch"]]
        pt += n * c["gamma_hat"] ** 2
        up += n * c["cp_simultaneous_upper"] ** 2
        detail[c["branch"]] = {"gamma_hat_iid": c["gamma_hat"],
                               "gammaU_sim_iid":
                                   c["cp_simultaneous_upper"]}
    return {"certificate_point_iid_n200": float(pt),
            "certificate_simultaneous_upper_iid_n200": float(up),
            "per_branch": detail,
            "source": "probe8_calibration.json (cycle 1), total_n=200 "
                      "cells, 50,000 i.i.d. null pairs per cell"}


# ----------------------------------------------------------------------
# Sections B/C -- fresh clustered worlds
# ----------------------------------------------------------------------
def _agg_half_exceed(out_list):
    hec = {}
    for out in out_list:
        for cell, (e, v) in out["half_exceed_by_cell"].items():
            pe, pv = hec.get(cell, (0, 0))
            hec[cell] = (pe + e, pv + v)
    return hec


def _summarise(out_list, counts, worlds):
    false_c = [o["n_promoted"] for o in out_list]
    p3_c = [o["n_p3_only"] for o in out_list]
    hec = _agg_half_exceed(out_list)
    hec_out = {}
    indep_pred = 0.0
    for cell, (e, v) in sorted(hec.items()):
        rate = e / v
        hec_out[cell] = {"exceedances": int(e), "evals": int(v),
                         "rate": float(rate)}
        indep_pred += counts[cell] * rate * rate
    mean_v = float(np.mean(false_c))
    mean_p3 = float(np.mean(p3_c))
    fwer = float(np.mean([1 if c > 0 else 0 for c in false_c]))
    return {
        "per_world_false": [int(c) for c in false_c],
        "mean_false": mean_v,
        "mean_false_se": float(np.std(false_c, ddof=1)
                               / math.sqrt(worlds)),
        "fwer": fwer,
        "fwer_se": float(math.sqrt(fwer * (1.0 - fwer) / worlds)),
        "fdr_equals_fwer_in_all_null": True,
        "per_world_p3_only": [int(c) for c in p3_c],
        "mean_p3_only": mean_p3,
        "p3_mean_se": float(np.std(p3_c, ddof=1) / math.sqrt(worlds)),
        "half_exceedance_by_cell": hec_out,
        "independence_prediction_P3_from_world_marginals":
            float(indep_pred),
        "dependence_factor_p3_over_independence_prediction":
            float(mean_p3 / indep_pred) if indep_pred > 0 else None}


def run_worlds(K_total, worlds, data_base, csplit_base,
               rowsplit_base=None):
    size = 200 // K_total
    cl_out, rw_out = [], []
    for w in range(worlds):
        rng = np.random.RandomState(data_base + w)
        cols, cl = P8.clustered_columns(rng, 200, K_total)
        pairs = P8.pair_index(cols)
        crng = np.random.RandomState(csplit_base + w)
        cperm = crng.permutation(K_total)
        h1 = np.isin(cl, cperm[:K_total // 2])
        cl_out.append(P8.run_world(cols, pairs, np.nonzero(h1)[0],
                                   np.nonzero(~h1)[0]))
        if rowsplit_base is not None:
            srng = np.random.RandomState(rowsplit_base + w)
            perm = srng.permutation(200)
            rw_out.append(P8.run_world(cols, pairs, perm[:100],
                                       perm[100:]))
    return cl_out, rw_out


def compare_block(summary, cal_design, counts):
    """Observed P3/conjoint vs point prediction, certificate, and the
    cycle-1 i.i.d. certificate (filled in by caller)."""
    c_pt = cal_design["certificate_point"]
    c_u = cal_design["certificate_simultaneous_upper"]
    se_cal = cal_design["certificate_point_delta_se"]
    se_comb = math.sqrt(summary["p3_mean_se"] ** 2 + se_cal ** 2)
    return {
        "point_prediction_sum_gamma_cl_sq": float(c_pt),
        "point_prediction_delta_se_calibration": float(se_cal),
        "certificate_CclU_simultaneous": float(c_u),
        "certificate_FWER_FDR_upper": float(min(1.0, c_u)),
        "p3_minus_point_in_combined_se":
            float((summary["mean_p3_only"] - c_pt) / se_comb),
        "p3_le_certificate": bool(summary["mean_p3_only"] <= c_u),
        "p3_excess_over_certificate_in_world_se":
            float((summary["mean_p3_only"] - c_u)
                  / summary["p3_mean_se"]),
        "conjoint_le_certificate":
            bool(summary["mean_false"] <= c_u),
        "conjoint_excess_over_certificate_in_world_se":
            float((summary["mean_false"] - c_u)
                  / summary["mean_false_se"])}


def section_b_primary(cal, counts):
    t0 = time.time()
    worlds = 60
    cl_out, rw_out = run_worlds(20, worlds, 930000, 940000,
                                rowsplit_base=950000)
    summ_c = _summarise(cl_out, counts, worlds)
    summ_r = _summarise(rw_out, counts, worlds)
    cal_p = cal["designs"]["P_primary_20x10"]
    res = {"design": "P: 20 clusters x 10 rows, n=200, ICC 0.5 numeric / "
                     "Dirichlet(3*balanced) categorical (cycle-1 design)",
           "worlds": worlds, "pairs_per_world": 780,
           "data_seed_base": 930000, "clustersplit_seed_base": 940000,
           "rowsplit_seed_base": 950000,
           "cluster_level_split": summ_c,
           "cluster_split_vs_cluster_certificate":
               compare_block(summ_c, cal_p, counts),
           "random_row_split_descriptive_NOT_certified": summ_r,
           "iid_certificate_contrast": iid_contrast(counts)}
    with open(os.path.join(HERE, "probe10_worlds_cluster.json"),
              "w") as fh:
        json.dump(res, fh, indent=2)
    print(f"[worlds P cluster split] meanV={summ_c['mean_false']:.3f} "
          f"P3mean={summ_c['mean_p3_only']:.3f} "
          f"Cpt={cal_p['certificate_point']:.3f} "
          f"CU={cal_p['certificate_simultaneous_upper']:.3f}", flush=True)
    print(f"[worlds P row split, descriptive] "
          f"meanV={summ_r['mean_false']:.3f} "
          f"P3mean={summ_r['mean_p3_only']:.3f}", flush=True)
    stamp("primary worlds", t0)
    return res


def section_c_sensitivity(cal, counts):
    t0 = time.time()
    worlds = 50
    cl_out, _ = run_worlds(40, worlds, 960000, 970000)
    summ = _summarise(cl_out, counts, worlds)
    cal_s = cal["designs"]["S_sens_40x5"]
    res = {"design": "S: 40 clusters x 5 rows, n=200, ICC 0.5 numeric / "
                     "Dirichlet(3*balanced) categorical; halves of "
                     "20 clusters x 5 rows (100 rows, as design P)",
           "worlds": worlds, "pairs_per_world": 780,
           "data_seed_base": 960000, "clustersplit_seed_base": 970000,
           "cluster_level_split": summ,
           "cluster_split_vs_cluster_certificate":
               compare_block(summ, cal_s, counts),
           "iid_certificate_contrast": iid_contrast(counts)}
    with open(os.path.join(HERE, "probe10_worlds_sensitivity.json"),
              "w") as fh:
        json.dump(res, fh, indent=2)
    print(f"[worlds S cluster split] meanV={summ['mean_false']:.3f} "
          f"P3mean={summ['mean_p3_only']:.3f} "
          f"Cpt={cal_s['certificate_point']:.3f} "
          f"CU={cal_s['certificate_simultaneous_upper']:.3f}", flush=True)
    stamp("sensitivity worlds", t0)
    return res


def main():
    t_all = time.time()
    section0_validation()
    counts = pair_counts()
    cal = section_a_calibration(counts)
    section_b_primary(cal, counts)
    section_c_sensitivity(cal, counts)
    print(f"ALL DONE in {time.time() - t_all:.1f}s", flush=True)


if __name__ == "__main__":
    main()
