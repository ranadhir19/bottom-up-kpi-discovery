"""Field CPSS head-to-head (protocol S6): Olist + Telco.

Analysis tables rebuilt EXACTLY as in probe5_olist_run.py /
probe7_telco_run.py (same construction, sample, seeds, controls).
Full-sample strengths/p-values, BH and paper-rule P3 use facade_v2;
CPSS resampling (B=100 complementary pairs) uses cpss_engine.
Paper-rule sets are checked against the stored probe 5 / probe 7
promoted sets (13 Olist / 88 Telco).
"""
import itertools
import json
import sys
import time

import numpy as np
import pandas as pd

PROBES = str(__import__("pathlib").Path(__file__).resolve().parents[2] / "probes")  # DEPOSIT COPY: archival absolute path replaced with a configurable one (see MANIFEST.md, "Deposit edits")
import os
FIELD = os.environ.get("BUKPI_FIELD_DATA", "data/field")  # DEPOSIT COPY: archival absolute path replaced with a configurable one (see MANIFEST.md, "Deposit edits")
sys.path.insert(0, PROBES)
sys.path.insert(0, ".")
import facade_v2 as F
import cpss_engine as E

TAU = 0.20
B = 100
RULES = ["paper_conjoint", "cpss_pi090", "cpss_pi080", "hybrid_pi090_bh"]

def corrected(rec):
    return rec.get("strength_corrected", rec["strength"])

# ----------------------------------------------------------------------
# table builders: probe 5 / probe 7 construction code, verbatim logic
# ----------------------------------------------------------------------
def build_olist():
    orders = pd.read_csv(f"{FIELD}/olist_orders_dataset.csv")
    items = pd.read_csv(f"{FIELD}/olist_order_items_dataset.csv")
    pays = pd.read_csv(f"{FIELD}/olist_order_payments_dataset.csv")
    revs = pd.read_csv(f"{FIELD}/olist_order_reviews_dataset.csv")
    cust = pd.read_csv(f"{FIELD}/olist_customers_dataset.csv")
    orders = orders[orders["order_status"] == "delivered"].copy()
    for c in ["order_purchase_timestamp", "order_delivered_customer_date",
              "order_estimated_delivery_date"]:
        orders[c] = pd.to_datetime(orders[c])
    it = items.groupby("order_id").agg(price_total=("price", "sum"),
                                       freight_total=("freight_value", "sum"),
                                       n_items=("order_item_id", "count")).reset_index()
    pay = pays.sort_values("payment_value").groupby("order_id").tail(1)[
        ["order_id", "payment_type", "payment_installments"]]
    payv = pays.groupby("order_id")["payment_value"].sum().rename("payment_value").reset_index()
    pay = pay.merge(payv, on="order_id").rename(columns={"payment_installments": "installments_max"})
    rv = revs[["order_id", "review_score", "review_comment_message"]].drop_duplicates("order_id")
    cu = cust[["customer_id", "customer_state"]]
    df = (orders.merge(it, on="order_id").merge(pay, on="order_id")
          .merge(rv, on="order_id").merge(cu, on="customer_id", how="left"))
    df["delivery_days"] = (df["order_delivered_customer_date"] - df["order_purchase_timestamp"]).dt.total_seconds() / 86400
    df["delay_days"] = (df["order_delivered_customer_date"] - df["order_estimated_delivery_date"]).dt.total_seconds() / 86400
    df["comment_len"] = df["review_comment_message"].fillna("").astype(str).str.len()
    df["has_comment"] = np.where(df["comment_len"] >= 10, "yes", "no")
    df["purchase_dow"] = df["order_purchase_timestamp"].dt.day_name()
    df = df.dropna(subset=["delivery_days", "delay_days", "customer_state"])
    n_full = len(df)
    if n_full > 25000:
        df = df.sample(n=25000, random_state=20261027)
    rng = np.random.RandomState(20261027)
    df["rand_num"] = rng.normal(0, 1, len(df))
    perm = np.random.RandomState(1)
    df["review_score_perm"] = perm.permutation(df["review_score"].values)
    df = df.reset_index(drop=True)
    NUM = ["price_total", "freight_total", "payment_value", "n_items",
           "installments_max", "review_score", "delivery_days", "delay_days",
           "comment_len", "rand_num", "review_score_perm"]
    CAT = ["payment_type", "customer_state", "purchase_dow", "has_comment"]
    types = {c: "numeric" for c in NUM}; types.update({c: "categorical" for c in CAT})
    # commented subset + text reductions, probe 5 verbatim
    sub = df[df["comment_len"] >= 10].copy()
    sub["review_band"] = np.where(sub["review_score"] <= 2, "low",
                                   np.where(sub["review_score"] == 3, "mid", "high"))
    Xt = F.tfidf_embed(sub["review_comment_message"].astype(str).tolist())
    sub["text_tfidf"] = F.kmeans(Xt, k=3, seed=42)
    sub = sub.reset_index(drop=True)
    return df, sub, NUM, CAT, types, n_full

def build_telco():
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
    types = {c: "numeric" for c in NUM}; types.update({c: "categorical" for c in CAT})
    return df, NUM, CAT, types, n_rows

# ----------------------------------------------------------------------
# generic machinery over "frames": {frame_name: (df, types_dict)}
# family entry: (frame, a, b)
# ----------------------------------------------------------------------
def full_sample(frame_map, family):
    rows = []
    for fr, a, b in family:
        df, types = frame_map[fr]
        r = F.correlate(df[a], df[b], types[a], types[b])
        rows.append({"frame": fr, "a": a, "b": b, "method": r["method"],
                     "strength_corr": float(abs(corrected(r))),
                     "p": float(r["p_value"]), "n": int(r["n"])})
    return rows

def splithalf_ok(frame_map, fr, a, b):
    df, types = frame_map[fr]
    idx = np.random.RandomState(7).permutation(len(df))
    for h in (idx[:len(df) // 2], idx[len(df) // 2:]):
        r = F.correlate(df[a].iloc[h], df[b].iloc[h], types[a], types[b])
        if abs(corrected(r)) < TAU:
            return False
    return True

def cpss_freq(frame_map, family, seeds):
    """One frequency per family entry; permutations per frame from the
    frame's declared seed (protocol S3/S6)."""
    by_frame = {}
    for i, (fr, a, b) in enumerate(family):
        by_frame.setdefault(fr, []).append(i)
    freq = np.zeros(len(family))
    for fr, idxs in by_frame.items():
        df, types = frame_map[fr]
        n = len(df)
        enc = {}
        for c in df.columns:
            if c in types:
                enc[c] = E.encode(df[c].tolist() if types[c] == "categorical" else df[c].values,
                                  types[c])
        counts = np.zeros(len(idxs), dtype=int)
        rng = np.random.RandomState(seeds[fr])
        for _ in range(B):
            perm = rng.permutation(n)
            for half in (perm[:n // 2], perm[n // 2:]):
                for j, i in enumerate(idxs):
                    _, a, b = family[i]
                    if E.corrected_abs(enc[a], enc[b], half) >= TAU:
                        counts[j] += 1
        for j, i in enumerate(idxs):
            freq[i] = counts[j] / (2 * B)
    return freq

def evaluate(dataset, frame_map, family, seeds, control_cols, reference_pairs):
    t0 = time.time()
    rows = full_sample(frame_map, family)
    pvals = np.array([r["p"] for r in rows])
    bh = F.bh_reject(pvals, 0.05)
    t_cpss0 = time.time()
    freq = cpss_freq(frame_map, family, seeds)
    cpss_seconds = time.time() - t_cpss0
    for r, b_, f_ in zip(rows, bh, freq):
        r["bh"] = bool(b_)
        r["stability_freq"] = float(f_)
    # rules
    sets = {rule: set() for rule in RULES}
    for i, r in enumerate(rows):
        key = frozenset((r["a"], r["b"]))
        if r["bh"] and r["strength_corr"] >= TAU and splithalf_ok(frame_map, r["frame"], r["a"], r["b"]):
            sets["paper_conjoint"].add(key)
        if r["stability_freq"] >= 0.90:
            sets["cpss_pi090"].add(key)
        if r["stability_freq"] >= 0.80:
            sets["cpss_pi080"].add(key)
        if r["stability_freq"] >= 0.90 and r["bh"]:
            sets["hybrid_pi090_bh"].add(key)
    lookup = {frozenset((r["a"], r["b"])): r for r in rows}
    def describe(key):
        r = lookup[key]
        return {"a": r["a"], "b": r["b"], "frame": r["frame"],
                "strength_corr": r["strength_corr"],
                "stability_freq": r["stability_freq"], "p": r["p"]}
    paper = sets["paper_conjoint"]
    ref = {frozenset(p) for p in reference_pairs}
    control_keys = [frozenset((r["a"], r["b"])) for r in rows
                    if r["a"] in control_cols or r["b"] in control_cols]
    per_rule = {}
    for rule in RULES:
        s = sets[rule]
        per_rule[rule] = {
            "promoted_count": len(s),
            "overlap_with_paper": len(s & paper),
            "promoted_pairs": sorted([sorted(k) for k in s]),
            "extra_vs_paper": [describe(k) for k in sorted(s - paper, key=lambda k: sorted(k))],
            "missed_vs_paper": [describe(k) for k in sorted(paper - s, key=lambda k: sorted(k))],
            "control_pairs_promoted": [sorted(k) for k in sorted((s & set(control_keys)), key=lambda k: sorted(k))],
            "control_pairs_promoted_count": len(s & set(control_keys)),
        }
    out = {"dataset": dataset,
           "family_pairs": len(rows),
           "control_pairs_in_family": len(control_keys),
           "max_stability_freq_among_controls": float(max(lookup[k]["stability_freq"] for k in control_keys)),
           "paper_set_matches_probe_reference": bool(paper == ref),
           "paper_reference_count": len(ref), "paper_recomputed_count": len(paper),
           "per_rule": per_rule,
           "pair_table": [{**r, "paper": bool(frozenset((r["a"], r["b"])) in paper),
                           "cpss_pi090": bool(frozenset((r["a"], r["b"])) in sets["cpss_pi090"]),
                           "cpss_pi080": bool(frozenset((r["a"], r["b"])) in sets["cpss_pi080"]),
                           "hybrid": bool(frozenset((r["a"], r["b"])) in sets["hybrid_pi090_bh"])}
                          for r in rows],
           "cpss_runtime_seconds": cpss_seconds,
           "total_runtime_seconds": time.time() - t0}
    return out

def main():
    results = {}
    # ---------- Olist ----------
    df, sub, NUM, CAT, types, n_full = build_olist()
    sub_types = dict(types); sub_types.update({"text_tfidf": "categorical", "review_band": "categorical"})
    frame_map = {"main": (df, types), "text": (sub, sub_types)}
    family = [("main", u, v) for u, v in itertools.combinations(NUM + CAT, 2)]
    family += [("text", "text_tfidf", "review_band"), ("text", "text_tfidf", "payment_type")]
    assert len(family) == 107, len(family)
    ref5 = json.load(open(f"{PROBES}/results/probe5_olist_run.json"))
    ref_pairs = [(p["a"], p["b"]) for p in ref5["promoted_final"]]
    res_o = evaluate("olist", frame_map, family,
                     {"main": 20261101, "text": 20261103},
                     {"rand_num", "review_score_perm"}, ref_pairs)
    res_o["analysis_rows_full"] = int(n_full); res_o["analysis_rows_used"] = int(len(df))
    res_o["commented_subset_rows"] = int(len(sub))
    results["olist"] = res_o
    pd.DataFrame(res_o["pair_table"]).to_csv("field_olist_pair_table.csv", index=False)
    print("OLIST done", res_o["per_rule"]["paper_conjoint"]["promoted_count"],
          res_o["paper_set_matches_probe_reference"], f"{res_o['total_runtime_seconds']:.1f}s", flush=True)
    # ---------- Telco ----------
    tdf, TNUM, TCAT, ttypes, n_rows = build_telco()
    tframe = {"main": (tdf, ttypes)}
    tfamily = [("main", u, v) for u, v in itertools.combinations(TNUM + TCAT, 2)]
    assert len(tfamily) == 231, len(tfamily)
    ref7 = json.load(open(f"{PROBES}/results/probe7_telco_run.json"))
    ref7_pairs = [(p["a"], p["b"]) for p in ref7["promoted_final"]]
    res_t = evaluate("telco", tframe, tfamily, {"main": 20261102},
                     {"rand_num", "monthly_perm"}, ref7_pairs)
    res_t["analysis_rows"] = int(n_rows)
    results["telco"] = res_t
    pd.DataFrame(res_t["pair_table"]).to_csv("field_telco_pair_table.csv", index=False)
    print("TELCO done", res_t["per_rule"]["paper_conjoint"]["promoted_count"],
          res_t["paper_set_matches_probe_reference"], f"{res_t['total_runtime_seconds']:.1f}s", flush=True)
    with open("field_results.json", "w") as fh:
        json.dump(results, fh, indent=2)
    print("FIELD TOTAL", results["olist"]["total_runtime_seconds"] + results["telco"]["total_runtime_seconds"])

if __name__ == "__main__":
    main()
