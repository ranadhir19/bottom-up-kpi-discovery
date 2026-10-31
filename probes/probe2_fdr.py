"""Probe 2 -- FDR as a TESTED pipeline component (Paper 8).

Mining worlds of mixed-type columns, n = 500 rows. Two worlds:
  NULL:    40 independent columns (24 numeric, 8 cat-3, 8 cat-2);
           all C(40,2) = 780 pairs are null.
  PLANTED: 46 columns; 6 planted pairs (rho .50/.35, eta2 .25/.10,
           V .45/.25), the other 1029 pairs null.
Rules compared per pair:
  R0  manuscript rule: raw strength >= 0.20
  R1  corrected strength >= 0.20 (V* / omega^2 / r)
  R2  uncorrected p < 0.05
  R3  BH q = 0.05 alone
  R4  BH q = 0.05 AND corrected strength >= 0.20   (v2 pipeline rule)
  R5  R4 AND persistence pass (60% of 20 half-sample windows at
      half-strength), evaluated on R4 survivors only, as in the pipeline
20 worlds per condition. Reports discovery counts, FWER (null world),
power and realised FDP (planted world).
"""
import json
import itertools
import numpy as np
import pandas as pd
import facade_v2 as F

N = 500
WORLDS = 20

def base_columns(rng):
    cols = {}
    for i in range(24):
        cols[f"num{i}"] = (rng.normal(0, 1, N), "numeric")
    for i in range(8):
        cols[f"cat3_{i}"] = (rng.choice(["a", "b", "c"], N), "categorical")
    for i in range(8):
        cols[f"cat2_{i}"] = (rng.choice(["u", "v"], N), "categorical")
    return cols

def planted_columns(rng):
    cols = base_columns(rng)
    # num-num plants
    x = rng.normal(0, 1, N)
    cols["pnum_a"] = (x, "numeric")
    cols["pnum_b"] = (0.5 * x + np.sqrt(1 - 0.25) * rng.normal(0, 1, N), "numeric")
    x2 = rng.normal(0, 1, N)
    cols["pnum_c"] = (x2, "numeric")
    cols["pnum_d"] = (0.35 * x2 + np.sqrt(1 - 0.35 ** 2) * rng.normal(0, 1, N), "numeric")
    # num-cat plants: eta2 .25 -> between var = 1/3 (within 1); eta2 .10 -> 1/9
    g = rng.choice(["a", "b", "c"], N)
    m25 = {"a": -np.sqrt(0.5), "b": 0.0, "c": np.sqrt(0.5)}  # var of means = 1/3
    cols["pcat_a"] = (g, "categorical")
    cols["pnum_e"] = (np.array([m25[v] for v in g]) + rng.normal(0, 1, N), "numeric")
    g2 = rng.choice(["a", "b", "c"], N)
    m10 = {"a": -np.sqrt(1 / 6), "b": 0.0, "c": np.sqrt(1 / 6)}  # var = 1/9
    cols["pcat_b"] = (g2, "categorical")
    cols["pnum_f"] = (np.array([m10[v] for v in g2]) + rng.normal(0, 1, N), "numeric")
    # cat-cat plants: V .45 / .25 via conditional diagonals solved on the fly
    def plant_v(vpop):
        from scipy.optimize import brentq
        def popv(p):
            off = (1 - p) / 2
            Pc = np.array([[p if i == j else off for j in range(3)] for i in range(3)])
            J = Pc / 3.0; E = np.full((3, 3), 1 / 9)
            return float(np.sqrt(float(((J - E) ** 2 / E).sum()) / 2))
        p = brentq(lambda q: popv(q) - vpop, 0.34, 0.999)
        a = rng.choice(["a", "b", "c"], N)
        idx = np.array([["a", "b", "c"].index(v) for v in a])
        u = rng.rand(N)
        off = (1 - p) / 2
        bidx = np.where(u < p, idx, np.where(u < p + off, (idx + 1) % 3, (idx + 2) % 3))
        return a, np.array([["a", "b", "c"][i] for i in bidx])
    a1, b1 = plant_v(0.45); cols["pcat_c"] = (a1, "categorical"); cols["pcat_d"] = (b1, "categorical")
    a2, b2 = plant_v(0.25); cols["pcat_e"] = (a2, "categorical"); cols["pcat_f"] = (b2, "categorical")
    truth = {frozenset(("pnum_a", "pnum_b")), frozenset(("pnum_c", "pnum_d")),
             frozenset(("pcat_a", "pnum_e")), frozenset(("pcat_b", "pnum_f")),
             frozenset(("pcat_c", "pcat_d")), frozenset(("pcat_e", "pcat_f"))}
    return cols, truth

def mine(cols):
    names = list(cols)
    recs = []
    for u, v in itertools.combinations(names, 2):
        (x, tx), (y, ty) = cols[u], cols[v]
        r = F.correlate(pd.Series(x), pd.Series(y), tx, ty)
        corr = r.get("strength_corrected", r["strength"])
        recs.append({"pair": frozenset((u, v)), "raw": abs(r["strength"]),
                     "corr": abs(corr), "p": r["p_value"],
                     "x": pd.Series(x), "y": pd.Series(y), "tx": tx, "ty": ty,
                     "corr_signed": corr})
    return recs

def apply_rules(recs, seed):
    p = np.array([r["p"] for r in recs])
    bh = F.bh_reject(p, 0.05)
    out = {}
    out["R0"] = {r["pair"] for r in recs if r["raw"] >= 0.20}
    out["R1"] = {r["pair"] for r in recs if r["corr"] >= 0.20}
    out["R2"] = {r["pair"] for r, b in zip(recs, p) if b < 0.05}
    out["R3"] = {r["pair"] for r, b in zip(recs, bh) if b}
    r4 = [r for r, b in zip(recs, bh) if b and r["corr"] >= 0.20]
    out["R4"] = {r["pair"] for r in r4}
    r5 = set()
    for r in r4:
        share, ok = F.persistence(r["x"], r["y"], r["tx"], r["ty"],
                                  r["corr_signed"], seed=seed)
        if ok: r5.add(r["pair"])
    out["R5"] = r5
    return out

res = {"null": {}, "planted": {}}
for cond in ["null", "planted"]:
    agg = {k: [] for k in ["R0", "R1", "R2", "R3", "R4", "R5"]}
    fwer = {k: 0 for k in agg}
    powers, fdps = [], []
    for w in range(WORLDS):
        rng = np.random.RandomState(31000 + w)
        if cond == "null":
            cols, truth = base_columns(rng), set()
        else:
            cols, truth = planted_columns(rng)
        recs = mine(cols)
        sel = apply_rules(recs, seed=777 + w)
        for k, s in sel.items():
            agg[k].append(len(s))
            if cond == "null" and len(s) > 0: fwer[k] += 1
        if cond == "planted":
            s4 = sel["R4"]
            powers.append(len(s4 & truth) / len(truth))
            fdps.append(len(s4 - truth) / max(len(s4), 1))
    res[cond] = {"mean_discoveries": {k: float(np.mean(v)) for k, v in agg.items()},
                 "max_discoveries": {k: int(np.max(v)) for k, v in agg.items()}}
    if cond == "null":
        res[cond]["FWER"] = {k: fwer[k] / WORLDS for k in fwer}
    else:
        res[cond]["R4_power_mean"] = float(np.mean(powers))
        res[cond]["R4_realised_FDP_mean"] = float(np.mean(fdps))
with open("results/probe2_fdr.json", "w") as fh:
    json.dump(res, fh, indent=2)
print(json.dumps(res, indent=2))
