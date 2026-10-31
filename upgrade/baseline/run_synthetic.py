"""Synthetic CPSS head-to-head (protocol S5).

World generators are copied VERBATIM from
hidden_files/paper08_work/probes/probe2_fdr.py (base_columns /
planted_columns). probe2_fdr.py runs its own experiment at import time,
so it cannot be imported as a library; the two generator functions are
reproduced character-for-character in logic and constants.
Full-sample strengths/p-values and the paper rule's split-half use
facade_v2 exactly as the probes do; CPSS resampling uses cpss_engine
(validated in engine_validation.json).
"""
import itertools
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

N = 500
WORLDS = 10          # per condition (protocol: >= 8)
B_PAIRS = 50         # complementary pairs per world (protocol S3 reduction)
TAU = 0.20
RULES = ["paper_conjoint", "cpss_pi090", "cpss_pi080", "hybrid_pi090_bh"]

# ----------------------------------------------------------------------
# probe2_fdr.py generators, verbatim
# ----------------------------------------------------------------------
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
            J = Pc / 3.0; E_ = np.full((3, 3), 1 / 9)
            return float(np.sqrt(float(((J - E_) ** 2 / E_).sum()) / 2))
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

# ----------------------------------------------------------------------
def corrected(rec):
    return rec.get("strength_corrected", rec["strength"])

def splithalf_ok(cols, u, v, n):
    """Paper rule P3, probes 5/7 form: RandomState(7) permutation,
    corrected strength >= tau on each disjoint half."""
    (x, tx), (y, ty) = cols[u], cols[v]
    idx = np.random.RandomState(7).permutation(n)
    for h in (idx[:n // 2], idx[n // 2:]):
        r = F.correlate(pd.Series(np.asarray(x)[h]), pd.Series(np.asarray(y)[h]), tx, ty)
        if abs(corrected(r)) < TAU:
            return False
    return True

def run_world(cond, w):
    t0 = time.time()
    rng = np.random.RandomState(31000 + w)
    if cond == "null":
        cols, truth = base_columns(rng), set()
    else:
        cols, truth = planted_columns(rng)
    names = list(cols)
    pair_list = list(itertools.combinations(names, 2))
    # ---- full sample via facade (probe2 mine) ----
    pvals = np.empty(len(pair_list))
    corr_abs = np.empty(len(pair_list))
    for i, (u, v) in enumerate(pair_list):
        (x, tx), (y, ty) = cols[u], cols[v]
        r = F.correlate(pd.Series(x), pd.Series(y), tx, ty)
        pvals[i] = r["p_value"]
        corr_abs[i] = abs(corrected(r))
    bh = F.bh_reject(pvals, 0.05)
    # ---- CPSS frequencies via engine ----
    enc = {nm: E.encode(cols[nm][0], cols[nm][1]) for nm in names}
    counts = np.zeros(len(pair_list), dtype=int)
    rng_c = np.random.RandomState(91000 + w)
    for _ in range(B_PAIRS):
        perm = rng_c.permutation(N)
        for half in (perm[:N // 2], perm[N // 2:]):
            for i, (u, v) in enumerate(pair_list):
                if E.corrected_abs(enc[u], enc[v], half) >= TAU:
                    counts[i] += 1
    freq = counts / (2 * B_PAIRS)
    # ---- rules ----
    sel = {}
    paper = set()
    for i, (u, v) in enumerate(pair_list):
        if bh[i] and corr_abs[i] >= TAU and splithalf_ok(cols, u, v, N):
            paper.add(frozenset((u, v)))
    sel["paper_conjoint"] = paper
    sel["cpss_pi090"] = {frozenset(pair_list[i]) for i in range(len(pair_list)) if freq[i] >= 0.90}
    sel["cpss_pi080"] = {frozenset(pair_list[i]) for i in range(len(pair_list)) if freq[i] >= 0.80}
    sel["hybrid_pi090_bh"] = {frozenset(pair_list[i]) for i in range(len(pair_list))
                              if freq[i] >= 0.90 and bh[i]}
    per_rule = {}
    for rule in RULES:
        s = sel[rule]
        false_disc = len(s - truth)
        per_rule[rule] = {
            "discoveries": int(len(s)),
            "false_discoveries": int(false_disc),
            "true_discoveries": int(len(s & truth)),
            "power": (float(len(s & truth) / len(truth)) if truth else None),
            "realised_fdp": float(false_disc / max(len(s), 1)),
        }
    return {"condition": cond, "world": w, "n_columns": len(names),
            "n_pairs": len(pair_list),
            "pair_names": [f"{u}||{v}" for u, v in pair_list],
            "stability_frequencies": [float(x) for x in freq],
            "full_corr_abs": [float(x) for x in corr_abs],
            "per_rule": per_rule,
            "runtime_seconds": time.time() - t0}

def aggregate(worlds, cond):
    agg = {}
    for rule in RULES:
        discs = [x["per_rule"][rule]["discoveries"] for x in worlds]
        falses = [x["per_rule"][rule]["false_discoveries"] for x in worlds]
        entry = {"mean_discoveries": float(np.mean(discs)),
                 "mean_false_discoveries": float(np.mean(falses)),
                 "max_discoveries": int(np.max(discs))}
        if cond == "null":
            entry["FWER"] = float(np.mean([d > 0 for d in discs]))
        else:
            entry["mean_power"] = float(np.mean([x["per_rule"][rule]["power"] for x in worlds]))
            entry["mean_realised_FDP"] = float(np.mean([x["per_rule"][rule]["realised_fdp"] for x in worlds]))
        agg[rule] = entry
    return agg

def main():
    t_all = time.time()
    out = {"design": {"n": N, "worlds_per_condition": WORLDS,
                      "B_complementary_pairs_synthetic": B_PAIRS,
                      "tau": TAU, "pi_primary": 0.90, "pi_sensitivity": 0.80,
                      "bh_q": 0.05,
                      "note_probe2_docstring": "probe2_fdr.py docstring says planted = 46 cols / 1029 other pairs; the verbatim code produces the counts recorded per world in this file."},
           "conditions": {}}
    for cond in ["null", "planted"]:
        worlds = []
        for w in range(WORLDS):
            res = run_world(cond, w)
            worlds.append(res)
            print(cond, w, {r: res["per_rule"][r]["discoveries"] for r in RULES},
                  f"{res['runtime_seconds']:.1f}s", flush=True)
        out["conditions"][cond] = {"n_columns": worlds[0]["n_columns"],
                                   "n_pairs": worlds[0]["n_pairs"],
                                   "aggregate": aggregate(worlds, cond),
                                   "worlds": worlds}
    out["total_runtime_seconds"] = time.time() - t_all
    with open("synthetic_results.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print(json.dumps({c: out["conditions"][c]["aggregate"] for c in out["conditions"]}, indent=2))
    print("TOTAL", out["total_runtime_seconds"])

if __name__ == "__main__":
    main()
