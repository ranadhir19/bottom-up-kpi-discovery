"""Probe 9 -- third-dataset replication run on UCI Bank Marketing.
Protocol frozen in bank_preregistration.md (written 2026-10-27 15:04 EDT,
before any data file was downloaded) before this script was written or run.
Mirrors probes 5/7 machinery exactly: facade_v2 correlate, BH q=0.05,
corrected-strength floor 0.20, split-half gate P3 (seed 7).
No sampling, no row filtering, no recoding; y is an ordinary column.
"""
import json
import itertools
import sys
import hashlib
import datetime
import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[2] / "probes"))  # DEPOSIT COPY: archival absolute path replaced with a configurable one (see MANIFEST.md, "Deposit edits")
import facade_v2 as F

THIRD = str(__import__("pathlib").Path(__file__).resolve().parent)  # DEPOSIT COPY: archival absolute path replaced with a configurable one (see MANIFEST.md, "Deposit edits")
CSV = f"{THIRD}/bank-additional/bank-additional-full.csv"
RUN_START = datetime.datetime.now(datetime.timezone.utc).isoformat()

df = pd.read_csv(CSV, sep=";")

EXPECTED_COLS = ["age", "job", "marital", "education", "default", "housing",
                 "loan", "contact", "month", "day_of_week", "duration",
                 "campaign", "pdays", "previous", "poutcome",
                 "emp.var.rate", "cons.price.idx", "cons.conf.idx",
                 "euribor3m", "nr.employed", "y"]
assert df.shape == (41188, 21), df.shape
assert list(df.columns) == EXPECTED_COLS, list(df.columns)
n_rows = len(df)

rng = np.random.RandomState(20261027)
df["rand_num"] = rng.normal(0, 1, len(df))
perm = np.random.RandomState(1)
df["age_perm"] = perm.permutation(df["age"].values)
df = df.reset_index(drop=True)

NUM = ["age", "duration", "campaign", "pdays", "previous",
       "emp.var.rate", "cons.price.idx", "cons.conf.idx", "euribor3m",
       "nr.employed", "rand_num", "age_perm"]
CAT = ["job", "marital", "education", "default", "housing", "loan",
       "contact", "month", "day_of_week", "poutcome", "y"]
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
assert len(fam) == 253, len(fam)
fam["bh"] = F.bh_reject(fam["p"].values, 0.05)
fam["floor"] = fam["strength_corr"] >= 0.20

def splithalf_ok(a, b):
    idx = np.random.RandomState(7).permutation(len(df))
    for h in (idx[:len(df) // 2], idx[len(df) // 2:]):
        r = F.correlate(df[a].iloc[h], df[b].iloc[h], types[a], types[b])
        if abs(corrected(r)) < 0.20:
            return False
    return True

HISTORY = {"pdays", "previous", "poutcome"}
MACRO = {"emp.var.rate", "cons.price.idx", "cons.conf.idx", "euribor3m",
         "nr.employed"}
CONTROLS = {"rand_num", "age_perm"}

def classify(a, b):
    pair = {a, b}
    if pair & CONTROLS:
        return "0_control"
    if "duration" in pair:
        return "a_leakage_duration"
    if pair <= HISTORY:
        return "b_campaign_history"
    if pair & MACRO:
        return "c_shared_context_macro"
    return "e_substantive"

promoted = []
for _, r in fam[fam["bh"] & fam["floor"]].iterrows():
    ok = splithalf_ok(r["a"], r["b"])
    promoted.append({**r.to_dict(), "splithalf": bool(ok),
                     "classification": classify(r["a"], r["b"])})
promoted_final = [p for p in promoted if p["splithalf"]]
controls_promoted = [p for p in promoted_final
                     if p["classification"] == "0_control"]
control_rows = fam[fam["a"].isin(CONTROLS) | fam["b"].isin(CONTROLS)]
from collections import Counter
class_counts = dict(Counter(p["classification"] for p in promoted_final))
top15 = fam.reindex(fam["strength_corr"].sort_values(ascending=False).index).head(15)

def clean(p):
    return {k: (float(v) if isinstance(v, (int, float, np.floating)) else v)
            for k, v in p.items()}

def file_hash(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

out = {
    "dataset": "UCI Bank Marketing, bank-additional-full.csv (official UCI archive, dataset 222)",
    "source_archive_url": "https://archive.ics.uci.edu/static/public/222/bank+marketing.zip",
    "source_dataset_page": "https://archive.ics.uci.edu/dataset/222/bank+marketing",
    "archive_sha256": file_hash(f"{THIRD}/bank+marketing.zip"),
    "inner_archive_sha256": file_hash(f"{THIRD}/bank-additional.zip"),
    "csv_sha256": file_hash(CSV),
    "preregistration_file": "bank_preregistration.md",
    "preregistration_written_edt": "2026-10-27 15:04 EDT (19:04 UTC), before any data download",
    "run_started_utc": RUN_START,
    "analysis_rows": int(n_rows),
    "analysis_columns_before_controls": 21,
    "analysis_columns_after_controls": 23,
    "columns": cols,
    "family_pairs": int(len(fam)),
    "bh_rejections": int(fam["bh"].sum()),
    "bh_and_floor": int((fam["bh"] & fam["floor"]).sum()),
    "promoted_final": [clean(p) for p in promoted_final],
    "promoted_final_count": len(promoted_final),
    "promoted_but_splithalf_failed": [
        {"a": p["a"], "b": p["b"], "strength_corr": float(p["strength_corr"]),
         "classification": p["classification"]}
        for p in promoted if not p["splithalf"]],
    "classification_counts": class_counts,
    "control_pairs_promoted": len(controls_promoted),
    "control_pairs_total": int(len(control_rows)),
    "control_pairs_bh_rejected": int(control_rows["bh"].sum()),
    "control_pairs_bh_and_floor": int((control_rows["bh"] & control_rows["floor"]).sum()),
    "control_pair_outcomes": [
        {"a": r["a"], "b": r["b"], "strength_corr": float(r["strength_corr"]),
         "p": float(r["p"]), "bh": bool(r["bh"]), "floor": bool(r["floor"])}
        for _, r in control_rows.iterrows()],
    "top15_by_corrected_strength": [
        {"a": r["a"], "b": r["b"], "method": r["method"],
         "strength_corr": float(r["strength_corr"]), "p": float(r["p"]),
         "bh": bool(r["bh"])} for _, r in top15.iterrows()],
    "deviations_from_preregistration": [],
}
out["run_finished_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
fam.to_csv(f"{THIRD}/probe9_full_pair_table.csv", index=False)
with open(f"{THIRD}/probe9_bank_run.json", "w") as fh:
    json.dump(out, fh, indent=2)
print(json.dumps(out, indent=2))
