"""Probe 1b -- where the V / eta^2 bias corrections actually bite: small n.

Cramer's V: symmetric 3x3 generators with population V in {0.73, 0.30},
n in {30, 60, 100, 200, 500}; 300 seeds each. Reports mean raw V, mean
Bergsma V*, and bias of each. eta^2: g = 6 groups, population eta^2 =
0.30, n in {60, 120, 600}; eta^2 vs omega^2.
"""
import json
import numpy as np
import facade_v2 as F

def pop_v_from_cond(p_diag):
    off = (1.0 - p_diag) / 2.0
    Pc = np.array([[p_diag if i == j else off for j in range(3)]
                   for i in range(3)])
    J = Pc / 3.0
    E = np.full((3, 3), 1.0 / 9.0)
    return float(np.sqrt(float(((J - E) ** 2 / E).sum()) / 2.0))

# conditional diag giving population V = 0.30: solve numerically
from scipy.optimize import brentq
p30 = brentq(lambda p: pop_v_from_cond(p) - 0.30, 0.34, 0.99)
print("cond diag for V=0.30:", round(p30, 4), " check V:", round(pop_v_from_cond(p30), 4))
print("cond diag 0.82 -> V:", round(pop_v_from_cond(0.82), 4))

def draw_cat(rng, n, p_diag, labels=("X", "Y", "Z")):
    a = rng.choice(labels, n)
    off = (1.0 - p_diag) / 2.0
    u = rng.rand(n)
    mapped = np.array([labels.index(v) for v in a])
    b_idx = np.where(u < p_diag, mapped,
                     np.where(u < p_diag + off, (mapped + 1) % 3,
                              (mapped + 2) % 3))
    return a, np.array([labels[i] for i in b_idx])

out = {"cramers": {}, "eta": {}}
for vpop, p_diag in [(0.73, 0.82), (0.30, p30)]:
    for n in [30, 60, 100, 200, 500]:
        vs, vss = [], []
        for s in range(300):
            rng = np.random.RandomState(10000 + s)
            a, b = draw_cat(rng, n, p_diag)
            ct = F.contingency(a, b)
            vs.append(F.cramers_v_from_table(ct)[0])
            vss.append(F.cramers_v_star_from_table(ct))
        out["cramers"][f"Vpop={vpop}_n={n}"] = {
            "V_raw_mean": float(np.mean(vs)), "V_raw_bias": float(np.mean(vs) - vpop),
            "Vstar_mean": float(np.mean(vss)), "Vstar_bias": float(np.mean(vss) - vpop),
            "V_raw_sd": float(np.std(vs))}
# eta^2 small-sample: g=6 means equally spaced to give pop eta2=0.30
# means m_j = d*(j-2.5); between var = d^2 * var(0..5)= d^2*35/12... choose d s.t. eta2=0.3
# between = (35/12) d^2 ; eta2 = between/(between+1) = 0.3 -> between = 3/7
d = np.sqrt((3.0 / 7.0) / (35.0 / 12.0))
means6 = d * (np.arange(6) - 2.5)
for n in [60, 120, 600]:
    es, os_ = [], []
    for s in range(300):
        rng = np.random.RandomState(20000 + s)
        cat = rng.choice(6, n)
        num = means6[cat] + rng.normal(0, 1, n)
        e, o, *_ = F.eta_omega(num, cat)
        es.append(e); os_.append(o)
    out["eta"][f"g=6_eta2pop=0.30_n={n}"] = {
        "eta2_mean": float(np.mean(es)), "eta2_bias": float(np.mean(es) - 0.30),
        "omega2_mean": float(np.mean(os_)), "omega2_bias": float(np.mean(os_) - 0.30)}
with open("results/probe1b_bias_smalln.json", "w") as fh:
    json.dump(out, fh, indent=2)
print(json.dumps(out, indent=2))
