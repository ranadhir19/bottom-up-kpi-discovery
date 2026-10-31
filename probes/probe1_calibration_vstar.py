"""Probe 1 -- multi-seed calibration audit + Bergsma V* (Paper 8).

Generators copied line-for-line from the shipped calibration facade
(journal repo 08, experiments/kpi_discovery_facade.py, seed stream
replicated per seed). Population targets are computed exactly from the
generators as coded -- NOT taken from the manuscript's stated ground
truths, two of which do not match the coded generators (see output).
100 seeds, n = 500 throughout.
"""
import json
import numpy as np
import pandas as pd
import facade_v2 as F

N = 500
SEEDS = range(100)

def pop_v(p_diag_cond):
    """Population Cramer's V for a 3x3 symmetric generator whose
    conditional diagonal is p_diag_cond and whose off-diagonal mass is
    split evenly; marginals uniform."""
    off = (1.0 - p_diag_cond) / 2.0
    Pc = np.array([[p_diag_cond if i == j else off for j in range(3)]
                   for i in range(3)])
    J = Pc / 3.0
    E = np.full((3, 3), 1.0 / 9.0)
    chi2n = float(((J - E) ** 2 / E).sum())
    return float(np.sqrt(chi2n / 2.0))

# population targets, generators as coded
POP = {
    "B1_rho": 0.5 / np.sqrt(0.25 + 0.85 ** 2),            # y = .5x + N(0,.85)
    "B2_eta2": 1.5 / (1.5 + 1.0),                          # means 0,1.5,3; sd 1
    # B3 as coded: copy w.p. .82 else uniform over all 3 labels
    "B3_V_coded": pop_v(0.82 + 0.18 / 3.0),
    # B3 variant: redraw excludes the matching label -> diag cond .82
    "B3_V_variant": pop_v(0.82),
}

def run_branch_set(seed):
    rng = np.random.RandomState(seed)
    out = {}
    x = rng.normal(0, 1, N); y = 0.5 * x + rng.normal(0, 0.85, N)
    out["r"], _ = F.pearson_r(x, y)
    cat = rng.choice(["A", "B", "C"], N)
    means = {"A": 0.0, "B": 1.5, "C": 3.0}
    num = np.array([rng.normal(means[c], 1.0) for c in cat])
    out["eta2"], out["omega2"], *_ = F.eta_omega(num, cat)
    a = rng.choice(["X", "Y", "Z"], N)
    flip = rng.rand(N) < 0.18
    b = np.where(flip, rng.choice(["P", "Q", "R"], N),
                 np.where(a == "X", "P", np.where(a == "Y", "Q", "R")))
    ct = F.contingency(a, b)
    out["V"], *_ = F.cramers_v_from_table(ct)
    out["Vstar"] = F.cramers_v_star_from_table(ct)
    # variant generator: redraw among the two non-matching labels
    flip2 = rng.rand(N) < 0.18
    alt = np.where(a == "X", rng.choice(["Q", "R"], N),
                   np.where(a == "Y", rng.choice(["P", "R"], N),
                            rng.choice(["P", "Q"], N)))
    b2 = np.where(flip2, alt, np.where(a == "X", "P",
                                       np.where(a == "Y", "Q", "R")))
    ct2 = F.contingency(a, b2)
    out["V_variant"], *_ = F.cramers_v_from_table(ct2)
    out["Vstar_variant"] = F.cramers_v_star_from_table(ct2)
    g1, g2 = rng.normal(0, 1, N), rng.normal(0, 1, N)
    out["null_r"], _ = F.pearson_r(g1, g2)
    pos = ["love this product", "great support team", "fast checkout",
           "very satisfied"]
    neg = ["terrible bug", "slow as molasses", "support never replied",
           "lost my data"]
    texts, label = [], []
    for _ in range(N):
        if rng.rand() < 0.5:
            texts.append(rng.choice(pos)); label.append(1)
        else:
            texts.append(rng.choice(neg)); label.append(0)
    cl = F.kmeans(np.vstack([F.hash_embed(t) for t in texts]), k=3, seed=seed)
    out["V_text"] = F.cramers_v_from_table(
        F.contingency(pd.Series(cl), pd.Series(label)))[0]
    return out

rows = [run_branch_set(s) for s in SEEDS]
df = pd.DataFrame(rows)

def summ(col, target=None):
    s = df[col]
    d = {"mean": float(s.mean()), "sd": float(s.std()),
         "q05": float(s.quantile(0.05)), "q95": float(s.quantile(0.95))}
    if target is not None:
        d["target"] = float(target)
        d["bias"] = float(s.mean() - target)
        d["rmse"] = float(np.sqrt(((s - target) ** 2).mean()))
    return d

# extended null audit: 1000 fresh pairs at seed stream 5000+
null_abs = []
for s in range(1000):
    rng = np.random.RandomState(5000 + s)
    r, _ = F.pearson_r(rng.normal(0, 1, N), rng.normal(0, 1, N))
    null_abs.append(abs(r))
null_abs = np.array(null_abs)

res = {
    "n": N, "seeds": len(list(SEEDS)),
    "population_targets": POP,
    "B1_pearson_r": summ("r", POP["B1_rho"]),
    "B2_eta2": summ("eta2", POP["B2_eta2"]),
    "B2_omega2": summ("omega2", POP["B2_eta2"]),
    "B3_V_raw_vs_coded_target": summ("V", POP["B3_V_coded"]),
    "B3_Vstar_vs_coded_target": summ("Vstar", POP["B3_V_coded"]),
    "B3_V_raw_vs_variant_target": summ("V_variant", POP["B3_V_variant"]),
    "B3_Vstar_vs_variant_target": summ("Vstar_variant", POP["B3_V_variant"]),
    "B3_V_raw_bias_if_target_were_0.733": float(df["V"].mean() - 0.733),
    "B6_text_V_hash64": summ("V_text"),
    "null_single_pair_abs_r": {
        "mean": float(np.abs(df["null_r"]).mean()),
        "max": float(np.abs(df["null_r"]).max())},
    "null_1000_pairs": {
        "mean_abs_r": float(null_abs.mean()),
        "max_abs_r": float(null_abs.max()),
        "n_exceeding_0.20": int((null_abs > 0.20).sum())},
    "manuscript_stated_truths": {"B1": 0.489, "B2": 0.60, "B3": 0.733},
}
with open("results/probe1_calibration_vstar.json", "w") as fh:
    json.dump(res, fh, indent=2)
print(json.dumps(res, indent=2))
