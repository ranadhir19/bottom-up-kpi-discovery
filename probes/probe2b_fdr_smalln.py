"""Probe 2b -- (a) null world at n = 60, categorical-heavy: where the
manuscript's strength-only rule collapses and BH + V* carries the load.
(b) per-plant detection record for rule R4 at n = 500 (5 worlds).
"""
import json
import itertools
import numpy as np
import pandas as pd
import facade_v2 as F
from probe2_fdr import base_columns, planted_columns, mine

# (a) small-n null world
res = {}
for n, worlds in [(60, 10)]:
    F_N = n
    agg = {k: [] for k in ["R0", "R1", "R3", "R4"]}
    for w in range(worlds):
        rng = np.random.RandomState(41000 + w)
        cols = {}
        import probe2_fdr as P2
        # rebuild base columns at this n (probe2 uses global N=500; do it inline)
        for i in range(16):
            cols[f"num{i}"] = (rng.normal(0, 1, n), "numeric")
        for i in range(16):
            cols[f"cat3_{i}"] = (rng.choice(["a", "b", "c"], n), "categorical")
        for i in range(8):
            cols[f"cat2_{i}"] = (rng.choice(["u", "v"], n), "categorical")
        names = list(cols)
        recs = []
        for u, v in itertools.combinations(names, 2):
            (x, tx), (y, ty) = cols[u], cols[v]
            r = F.correlate(pd.Series(x), pd.Series(y), tx, ty)
            recs.append({"pair": (u, v), "raw": abs(r["strength"]),
                         "corr": abs(r.get("strength_corrected", r["strength"])),
                         "p": r["p_value"]})
        p = np.array([r["p"] for r in recs])
        bh = F.bh_reject(p, 0.05)
        agg["R0"].append(sum(r["raw"] >= 0.20 for r in recs))
        agg["R1"].append(sum(r["corr"] >= 0.20 for r in recs))
        agg["R3"].append(int(bh.sum()))
        agg["R4"].append(sum(1 for r, b in zip(recs, bh) if b and r["corr"] >= 0.20))
    res[f"null_n={n}"] = {k: {"mean": float(np.mean(v)), "max": int(np.max(v))}
                          for k, v in agg.items()}

# (b) per-plant detection under R4 at n = 500
det = {}
for w in range(5):
    rng = np.random.RandomState(31000 + w)
    cols, truth = planted_columns(rng)
    recs = mine(cols)
    p = np.array([r["p"] for r in recs])
    bh = F.bh_reject(p, 0.05)
    found = {r["pair"] for r, b in zip(recs, bh) if b and r["corr"] >= 0.20}
    for t in truth:
        key = " & ".join(sorted(t))
        det.setdefault(key, []).append(t in found)
res["plant_detection_R4"] = {k: {"detected": f"{sum(v)}/{len(v)}"} for k, v in det.items()}
with open("results/probe2b_fdr_smalln.json", "w") as fh:
    json.dump(res, fh, indent=2)
print(json.dumps(res, indent=2))
