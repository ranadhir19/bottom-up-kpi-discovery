"""Probe 7 -- strengthening replication run on IBM Telco (second dataset).
Protocol frozen in ../telco_preregistration.md before this script ran.
Responds to Hypothesis Tester defect D8 (single field dataset).
Mirrors probe 5's machinery exactly: facade_v2 correlate, BH q=0.05,
corrected-strength floor 0.20, split-half gate (seed 7).
"""
import json
import itertools
import numpy as np
import pandas as pd
import facade_v2 as F

import os
FIELD = os.environ.get("BUKPI_FIELD_DATA", "data/field")  # DEPOSIT COPY: archival absolute path replaced with a configurable one (see MANIFEST.md, "Deposit edits")

df = pd.read_csv(f"{FIELD}/telco_churn.csv")
df["TotalCharges"] = pd.to_numeric(df["TotalCharges"], errors="coerce")
n_rows = len(df)

rng = np.random.RandomState(20261027)
df["rand_num"] = rng.normal(0, 1, len(df))
perm = np.random.RandomState(1)
df["monthly_perm"] = perm.permutation(df["MonthlyCharges"].values)
df = df.reset_index(drop=True)

NUM = ["tenure", "MonthlyCharges", "TotalCharges", "rand_num", "monthly_perm"]
CAT = ["gender", "SeniorCitizen", "Partner", "Dependents", "PhoneService",
       "MultipleLines", "InternetService", "OnlineSecurity", "OnlineBackup",
       "DeviceProtection", "TechSupport", "StreamingTV", "StreamingMovies",
       "Contract", "PaperlessBilling", "PaymentMethod", "Churn"]
types = {c: "numeric" for c in NUM}
types.update({c: "categorical" for c in CAT})

def corrected(rec):
    return rec.get("strength_corrected", rec["strength"])

rows = []
cols = NUM + CAT
for u, v in itertools.combinations(cols, 2):
    r = F.correlate(df[u], df[v], types[u], types[v])
    rows.append({"a": u, "b": v, "method": r["method"],
                 "strength_raw": abs(r["strength"]),
                 "strength_corr": abs(corrected(r)), "p": r["p_value"],
                 "n": r["n"]})
fam = pd.DataFrame(rows)
assert len(fam) == 231, len(fam)
fam["bh"] = F.bh_reject(fam["p"].values, 0.05)
fam["floor"] = fam["strength_corr"] >= 0.20

def splithalf_ok(a, b):
    idx = np.random.RandomState(7).permutation(len(df))
    for h in (idx[:len(df) // 2], idx[len(df) // 2:]):
        r = F.correlate(df[a].iloc[h], df[b].iloc[h], types[a], types[b])
        if abs(corrected(r)) < 0.20:
            return False
    return True

ADDONS = {"OnlineSecurity", "OnlineBackup", "DeviceProtection",
          "TechSupport", "StreamingTV", "StreamingMovies"}

def classify(a, b):
    pair = {a, b}
    if "TotalCharges" in pair and pair & {"tenure", "MonthlyCharges"}:
        return "a_accumulated_billing"
    if ("InternetService" in pair and pair & ADDONS) or (pair <= ADDONS):
        return "b_structural_service_level"
    return "d_substantive"

promoted = []
for _, r in fam[fam["bh"] & fam["floor"]].iterrows():
    ok = splithalf_ok(r["a"], r["b"])
    promoted.append({**r.to_dict(), "splithalf": bool(ok),
                     "classification": classify(r["a"], r["b"])})
promoted_final = [p for p in promoted if p["splithalf"]]
controls = [p for p in promoted_final
            if "rand_num" in (p["a"], p["b"]) or "monthly_perm" in (p["a"], p["b"])]
top15 = fam.reindex(fam["strength_corr"].sort_values(ascending=False).index).head(15)

out = {
    "dataset": "IBM Telco customer churn (telco_churn.csv)",
    "analysis_rows": int(n_rows),
    "family_pairs": int(len(fam)),
    "bh_rejections": int(fam["bh"].sum()),
    "bh_and_floor": int((fam["bh"] & fam["floor"]).sum()),
    "promoted_final": [{k: (float(v) if isinstance(v, (int, float, np.floating)) else v)
                        for k, v in p.items()} for p in promoted_final],
    "promoted_but_splithalf_failed": [
        {"a": p["a"], "b": p["b"], "strength_corr": p["strength_corr"]}
        for p in promoted if not p["splithalf"]],
    "control_pairs_promoted": len(controls),
    "top15_by_corrected_strength": [
        {"a": r["a"], "b": r["b"], "method": r["method"],
         "strength_corr": float(r["strength_corr"]), "p": float(r["p"]),
         "bh": bool(r["bh"])} for _, r in top15.iterrows()],
}
fam.to_csv("results/probe7_full_pair_table.csv", index=False)
with open("results/probe7_telco_run.json", "w") as fh:
    json.dump(out, fh, indent=2)
print(json.dumps(out, indent=2))
