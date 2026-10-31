"""Probe 3b -- disjoint split-half reproducibility, the replacement gate.

Rule: split rows into two disjoint halves; the pair passes iff the
(corrected) strength in EACH half meets the reportable floor 0.20
(for numeric pairs, |r| >= 0.20 in both halves). Halves are disjoint,
so under the null the two events are independent -- unlike the
overlapping subsample windows of the original persistence rule.
200 pairs per cell.
"""
import json
import numpy as np
import pandas as pd
import facade_v2 as F

def splithalf_pass(a, b, ta, tb, floor=0.20, seed=5):
    a = pd.Series(a).reset_index(drop=True); b = pd.Series(b).reset_index(drop=True)
    n = len(a)
    rng = np.random.RandomState(seed)
    idx = rng.permutation(n)
    h1, h2 = idx[:n // 2], idx[n // 2:]
    ok = True
    for h in (h1, h2):
        r = F.correlate(a.iloc[h], b.iloc[h], ta, tb)
        s = abs(r.get("strength_corrected", r["strength"]))
        if s < floor: ok = False
    return ok

def run(kind, n, reps=200):
    hits = 0
    for s in range(reps):
        rng = np.random.RandomState(70000 + s)
        if kind == "null_num":
            a = rng.normal(0, 1, n); b = rng.normal(0, 1, n); ta = tb = "numeric"
        elif kind == "plant_rho03":
            a = rng.normal(0, 1, n); b = 0.3 * a + np.sqrt(0.91) * rng.normal(0, 1, n); ta = tb = "numeric"
        elif kind == "plant_rho05":
            a = rng.normal(0, 1, n); b = 0.5 * a + np.sqrt(0.75) * rng.normal(0, 1, n); ta = tb = "numeric"
        elif kind == "null_cat":
            a = rng.choice(["a", "b", "c"], n); b = rng.choice(["p", "q", "r"], n); ta = tb = "categorical"
        elif kind == "plant_V030":
            a = rng.choice(["a", "b", "c"], n)
            idx = np.array([["a", "b", "c"].index(v) for v in a])
            p = 0.5333; off = (1 - p) / 2
            u = rng.rand(n)
            bidx = np.where(u < p, idx, np.where(u < p + off, (idx + 1) % 3, (idx + 2) % 3))
            b = np.array([["a", "b", "c"][i] for i in bidx]); ta = tb = "categorical"
        if splithalf_pass(pd.Series(a), pd.Series(b), ta, tb, seed=300 + s):
            hits += 1
    return hits / reps

res = {}
for kind in ["null_num", "plant_rho03", "plant_rho05", "null_cat", "plant_V030"]:
    for n in [100, 500]:
        res[f"{kind}_n={n}"] = {"pass_rate": run(kind, n)}
with open("results/probe3b_splithalf.json", "w") as fh:
    json.dump(res, fh, indent=2)
print(json.dumps(res, indent=2))
