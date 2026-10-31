"""Probe 3 -- the persistence test as a tested component, and an audit of
the addendum's operating-region claim (alpha_eff = 0.4 * 0.05 = 0.02).

Persistence (manuscript Sec. 3.8): 20 windows, each a 50% subsample;
a window passes if its strength >= 0.5 x full-sample strength; the pair
persists if >= 60% of windows pass. We measure, per pair class:
  pi  = mean per-window pass probability;
  pass rate of the 60% rule;
for null pairs (the false-pass rate is the per-edge Type-I error of the
persistence gate alone) and planted pairs (power cost). n in {100, 500}.
"""
import json
import numpy as np
import pandas as pd
import facade_v2 as F

def persist_stats(a, b, ta, tb, reps_seed):
    rec = F.correlate(pd.Series(a), pd.Series(b), ta, tb)
    base = rec.get("strength_corrected", rec["strength"])
    base = abs(base)
    rng = np.random.RandomState(reps_seed)
    n = len(a)
    w = max(int(0.5 * n), 30)
    hits = 0
    for _ in range(20):
        idx = rng.choice(n, w, replace=False)
        r2 = F.correlate(pd.Series(a).iloc[idx], pd.Series(b).iloc[idx], ta, tb)
        s = abs(r2.get("strength_corrected", r2["strength"]))
        if s >= 0.5 * base: hits += 1
    return hits / 20.0, base

def run(kind, n, reps=200):
    shares, bases = [], []
    for s in range(reps):
        rng = np.random.RandomState(60000 + s)
        if kind == "null_num":
            a = rng.normal(0, 1, n); b = rng.normal(0, 1, n); ta = tb = "numeric"
        elif kind == "plant_rho03":
            a = rng.normal(0, 1, n); b = 0.3 * a + np.sqrt(0.91) * rng.normal(0, 1, n); ta = tb = "numeric"
        elif kind == "plant_rho05":
            a = rng.normal(0, 1, n); b = 0.5 * a + np.sqrt(0.75) * rng.normal(0, 1, n); ta = tb = "numeric"
        elif kind == "null_cat":
            a = rng.choice(["a", "b", "c"], n); b = rng.choice(["p", "q", "r"], n)
            ta = tb = "categorical"
        share, base = persist_stats(pd.Series(a), pd.Series(b), ta, tb, 900 + s)
        shares.append(share); bases.append(base)
    shares = np.array(shares); bases = np.array(bases)
    out = {"pi_mean_window_pass": float(shares.mean()),
           "pass_rate_60pct_rule": float((shares >= 0.60).mean()),
           "mean_full_strength": float(bases.mean())}
    if kind == "null_num":
        m = bases >= 0.20
        out["n_pairs_with_full_strength>=0.20"] = int(m.sum())
        if m.sum() > 0:
            out["pass_rate_given_strength>=0.20"] = float((shares[m] >= 0.60).mean())
    return out

res = {}
for kind in ["null_num", "plant_rho03", "plant_rho05", "null_cat"]:
    for n in ([100, 500] if kind != "null_cat" else [500]):
        res[f"{kind}_n={n}"] = run(kind, n)
res["addendum_claim"] = {"alpha_eff_claimed": 0.02,
                         "basis": "0.4 * 0.05 under a 60% window-survival rule (binomial)"}
with open("results/probe3_persistence.json", "w") as fh:
    json.dump(res, fh, indent=2)
print(json.dumps(res, indent=2))
