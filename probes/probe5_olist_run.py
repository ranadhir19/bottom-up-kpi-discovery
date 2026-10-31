"""Probe 5 -- pre-registered end-to-end discovery run on Olist (real data).
Protocol frozen in ../olist_preregistration.md before this script ran.
"""
import json
import itertools
import numpy as np
import pandas as pd
import facade_v2 as F

import os
FIELD = os.environ.get("BUKPI_FIELD_DATA", "data/field")  # DEPOSIT COPY: archival absolute path replaced with a configurable one (see MANIFEST.md, "Deposit edits")

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

# text pairs on commented subset
sub = df[df["comment_len"] >= 10].copy()
sub["review_band"] = np.where(sub["review_score"] <= 2, "low",
                               np.where(sub["review_score"] == 3, "mid", "high"))
Xt = F.tfidf_embed(sub["review_comment_message"].astype(str).tolist())
sub["text_tfidf"] = F.kmeans(Xt, k=3, seed=42)
Xh = np.vstack([F.hash_embed(t, 64) for t in sub["review_comment_message"].astype(str)])
sub["text_hash"] = F.kmeans(Xh, k=3, seed=42)
for red, tgt in [("text_tfidf", "review_band"), ("text_tfidf", "payment_type"),
                 ("text_hash", "review_band"), ("text_hash", "payment_type")]:
    r = F.correlate(sub[red], sub[tgt], "categorical",
                    "numeric" if tgt == "review_score" else "categorical")
    rows.append({"a": red, "b": tgt, "method": r["method"],
                 "strength_raw": abs(r["strength"]),
                 "strength_corr": abs(corrected(r)), "p": r["p_value"],
                 "n": r["n"], "text_pair": True})

tab = pd.DataFrame(rows)
# mining family = the 105 column pairs + 2 TF-IDF text pairs (pre-registered 107)
fam = tab[(~tab.get("text_pair", False).fillna(False)) |
          (tab["a"] == "text_tfidf")].copy()
assert len(fam) == 107, len(fam)
fam["bh"] = F.bh_reject(fam["p"].values, 0.05)
fam["floor"] = fam["strength_corr"] >= 0.20

def splithalf_ok(a, b):
    if a.startswith("text_") or b.startswith("text_"):
        d = sub; acol, bcol = a, b
    else:
        d = df; acol, bcol = a, b
    idx = np.random.RandomState(7).permutation(len(d))
    for h in (idx[:len(d) // 2], idx[len(d) // 2:]):
        r = F.correlate(d[acol].iloc[h], d[bcol].iloc[h],
                        types.get(acol, "categorical"), types.get(bcol, "categorical"))
        if abs(corrected(r)) < 0.20: return False
    return True

promoted = []
for _, r in fam[fam["bh"] & fam["floor"]].iterrows():
    ok = splithalf_ok(r["a"], r["b"])
    # descriptive persistence share for promoted pairs (not a gate)
    d = sub if r["a"].startswith("text_") else df
    share, _ = F.persistence(d[r["a"]], d[r["b"]],
                             types.get(r["a"], "categorical"),
                             types.get(r["b"], "categorical"),
                             r["strength_corr"], seed=99)
    promoted.append({**r.to_dict(), "splithalf": bool(ok),
                     "persistence_share_orig_rule": float(share)})

promoted_final = [p for p in promoted if p["splithalf"]]
controls = [p for p in promoted_final
            if "rand_num" in (p["a"], p["b"]) or "review_score_perm" in (p["a"], p["b"])]
top15 = fam.reindex(fam["strength_corr"].sort_values(ascending=False).index).head(15)
out = {
    "analysis_rows_full": int(n_full), "analysis_rows_used": int(len(df)),
    "commented_subset_rows": int(len(sub)),
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
    "text_pairs_both_reductions": [
        {"a": r["a"], "b": r["b"], "strength_corr": float(r["strength_corr"]),
         "p": float(r["p"])} for _, r in tab[tab.get("text_pair", False).fillna(False)].iterrows()],
}
tab.to_csv("results/probe5_full_pair_table.csv", index=False)
with open("results/probe5_olist_run.json", "w") as fh:
    json.dump(out, fh, indent=2)
print(json.dumps(out, indent=2))
