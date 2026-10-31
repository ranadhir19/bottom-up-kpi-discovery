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

# text pairs on commented subset
sub = df[df["comment_len"] >= 10].copy()
sub["review_band"] = np.where(sub["review_score"] <= 2, "low",
                               np.where(sub["review_score"] == 3, "mid", "high"))
Xt = F.tfidf_embed(sub["review_comment_message"].astype(str).tolist())
sub["text_tfidf"] = F.kmeans(Xt, k=3, seed=42)
Xh = np.vstack([F.hash_embed(t, 64) for t in sub["review_comment_message"].astype(str)])
sub["text_hash"] = F.kmeans(Xh, k=3, seed=42)
for red, tgt in [("text_tfidf", "review_band"), ("text_tfidf", "payment_type"),

# --- probe 5b proper: signed strengths + abs-fixed persistence shares ---
def corr_signed(rec): return rec.get("strength_corrected", rec["strength"])
promoted_pairs = [("price_total","freight_total"),("price_total","payment_value"),
("price_total","installments_max"),("freight_total","payment_value"),
("freight_total","n_items"),("freight_total","installments_max"),
("payment_value","installments_max"),("review_score","delivery_days"),
("review_score","delay_days"),("review_score","comment_len"),
("delivery_days","delay_days"),("comment_len","has_comment"),
("text_tfidf","review_band")]
out=[]
for a,b in promoted_pairs:
    d = sub if a.startswith("text_") else df
    ta = types.get(a,"categorical"); tb = types.get(b,"categorical")
    r = F.correlate(d[a], d[b], ta, tb)
    share,_ = F.persistence(d[a], d[b], ta, tb, corr_signed(r), seed=99)
    out.append({"a":a,"b":b,"signed_strength":float(corr_signed(r)),
                "persistence_share_fixed":float(share)})
print(json.dumps(out, indent=2))
json.dump(out, open("results/probe5b_promoted_signed.json","w"), indent=2)
