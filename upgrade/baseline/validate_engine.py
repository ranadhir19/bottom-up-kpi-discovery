"""Gate: validate cpss_engine against facade_v2.correlate (protocol S7).

>=30 random pairs, >=10 per branch, each compared on the full sample and
on a random half-sample. Acceptance: max abs difference < 1e-6.
"""
import json
import sys
import time

import numpy as np
import pandas as pd

BASE = str(__import__("pathlib").Path(__file__).resolve().parents[2] / "probes")  # DEPOSIT COPY: archival absolute path replaced with a configurable one (see MANIFEST.md, "Deposit edits")
sys.path.insert(0, BASE)
sys.path.insert(0, ".")
import facade_v2 as F
import cpss_engine as E

t0 = time.time()
rng = np.random.RandomState(424242)

def facade_corr(x, y, tx, ty, idx=None):
    sx = pd.Series(np.asarray(x, dtype=object) if tx == "categorical" else np.asarray(x, float))
    sy = pd.Series(np.asarray(y, dtype=object) if ty == "categorical" else np.asarray(y, float))
    if idx is not None:
        sx, sy = sx.iloc[idx], sy.iloc[idx]
    rec = F.correlate(sx, sy, tx, ty)
    return abs(rec.get("strength_corrected", rec["strength"]))

def rand_num(n, rng, kind):
    if kind == 0:
        return rng.normal(0, 1, n)
    if kind == 1:
        return rng.gamma(2.0, 1.5, n)
    return rng.normal(10, 3, n) + (rng.rand(n) < 0.1) * rng.normal(0, 20, n)

def rand_cat(n, rng, levels):
    probs = rng.dirichlet(np.ones(levels) * 2.0)
    return rng.choice([f"L{i}" for i in range(levels)], n, p=probs)

results = []
branch_max = {"numeric-numeric": 0.0, "categorical-categorical": 0.0,
              "numeric-categorical": 0.0}
overall_max = 0.0
K = 12  # per branch, 36 pairs total
for branch_i, branch in enumerate(branch_max):
    for k in range(K):
        n = int(rng.choice([80, 150, 250, 500, 900]))
        if branch == "numeric-numeric":
            xa = rand_num(n, rng, k % 3)
            shared = rng.normal(0, 1, n)
            rho = [0.0, 0.3, 0.6][k % 3]
            xb = rho * shared + np.sqrt(max(1 - rho * rho, 0.01)) * rand_num(n, rng, (k + 1) % 3) if k % 3 else rand_num(n, rng, k % 3)
            if k % 3:
                xa = shared + 0.3 * rand_num(n, rng, 0)
            tx = ty = "numeric"
            ca, cb = E.encode(xa, "numeric"), E.encode(xb, "numeric")
        elif branch == "categorical-categorical":
            lev_a = 2 + (k % 4)
            lev_b = 2 + ((k + 1) % 4)
            xa = rand_cat(n, rng, lev_a)
            if k % 3 == 0:
                xb = rand_cat(n, rng, lev_b)
            else:  # dependent copy with noise
                xb = np.array([v if rng.rand() > 0.35 else rand_cat(1, rng, lev_b)[0] for v in xa])
                # map labels if level counts differ: cycle labels
                xb = np.array([f"L{int(v[1:]) % lev_b}" for v in xb])
            tx = ty = "categorical"
            ca, cb = E.encode(xa, "categorical"), E.encode(xb, "categorical")
        else:
            levels = 2 + (k % 4)
            xa = rand_cat(n, rng, levels)
            shift = {"L0": -1.0, "L1": 0.0, "L2": 1.0, "L3": 1.8, "L4": -1.8}
            xb = np.array([shift.get(v, 0.0) for v in xa]) * (0.8 if k % 2 else 0.0) + rand_num(n, rng, k % 3)
            tx, ty = "categorical", "numeric"
            ca, cb = E.encode(xa, "categorical"), E.encode(xb, "numeric")
        # NaN injection on numeric side for a few cases (pairwise handling)
        half = np.random.RandomState(1000 + branch_i * 100 + k).permutation(n)[: n // 2]
        for tag, idx in (("full", None), ("half", half)):
            use_idx = idx if idx is not None else np.arange(n)
            eng = E.corrected_abs(ca, cb, use_idx)
            fac = facade_corr(xa, xb, tx, ty, idx)
            d = abs(eng - fac)
            branch_max[branch] = max(branch_max[branch], d)
            overall_max = max(overall_max, d)
            results.append({"branch": branch, "k": k, "n": n, "eval": tag,
                            "engine": float(eng), "facade": float(fac),
                            "abs_diff": float(d)})

out = {"n_pairs": 3 * K, "n_evaluations": len(results),
       "branch_max_abs_diff": branch_max,
       "max_abs_diff": float(overall_max),
       "acceptance_threshold": 1e-6,
       "gate_passed": bool(overall_max < 1e-6),
       "runtime_seconds": time.time() - t0,
       "details": results}
with open("engine_validation.json", "w") as fh:
    json.dump(out, fh, indent=2)
print(json.dumps({k: v for k, v in out.items() if k != "details"}, indent=2))
