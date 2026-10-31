"""Probe 6 -- phi_K cross-check on the Olist field family (D7 remedy).
Rebuilds probe 5's analysis table with the identical construction and
seeds, then computes Baak et al. (2020) phi_K for the 105 non-text
pairs with the phik library (v0.12.5, default settings, noise
correction on) and the library's chi-square significance matrix.
Analogue promotion sets: A = phi_K >= 0.20; B = A and BH(q=0.05)
over chi-square p-values computed with SciPy on the same binned
contingency tables phik uses internally (phik's own
significance_matrix returns z-scores, not p-values, in v0.12.5;
the binning is phik's `bin_data` with default settings).
Compared against the paper's 12 non-text promoted pairs from
probe 5 / probe 5b.
Run with the phik venv python.
"""
import json
import itertools
import numpy as np
import pandas as pd
import phik  # noqa: F401  (registers pandas accessors)
from phik import phik_matrix
from phik.binning import bin_data
from scipy.stats import chi2_contingency

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
cols = NUM + CAT
sub = df[cols].copy()

pm = phik_matrix(sub, interval_cols=NUM)
binned = bin_data(sub, cols=NUM)

def chi2_p(u, v):
    ct = pd.crosstab(binned[u], binned[v])
    return float(chi2_contingency(ct.values, correction=False)[1])

def bh_reject(pvals, q=0.05):
    p = np.asarray(pvals, dtype=float)
    m = len(p)
    order = np.argsort(p)
    ranked = p[order]
    thresh = q * (np.arange(1, m + 1)) / m
    below = ranked <= thresh
    reject = np.zeros(m, dtype=bool)
    if below.any():
        k = np.max(np.nonzero(below)[0])
        reject[order[:k + 1]] = True
    return reject

pairs = []
for u, v in itertools.combinations(cols, 2):
    pairs.append({"a": u, "b": v, "phik": float(pm.loc[u, v]),
                  "p": chi2_p(u, v)})
tab = pd.DataFrame(pairs)
assert len(tab) == 105, len(tab)
tab["bh"] = bh_reject(tab["p"].values, 0.05)
tab["floor"] = tab["phik"] >= 0.20
setA = {tuple(sorted((r.a, r.b))) for r in tab[tab["floor"]].itertuples()}
setB = {tuple(sorted((r.a, r.b))) for r in tab[tab["floor"] & tab["bh"]].itertuples()}

paper13 = [("price_total", "freight_total"), ("price_total", "payment_value"),
           ("price_total", "installments_max"), ("freight_total", "payment_value"),
           ("freight_total", "n_items"), ("freight_total", "installments_max"),
           ("payment_value", "installments_max"), ("review_score", "delivery_days"),
           ("review_score", "delay_days"), ("review_score", "comment_len"),
           ("delivery_days", "delay_days"), ("comment_len", "has_comment")]
paper12 = {tuple(sorted(p)) for p in paper13}
phik_of = {tuple(sorted((r.a, r.b))): r.phik for r in tab.itertuples()}

out = {
    "phik_version": phik.__version__,
    "analysis_rows_full": int(n_full), "analysis_rows_used": int(len(df)),
    "family_pairs": 105,
    "setA_phik_floor_only": sorted(map(list, setA)),
    "setB_phik_floor_and_BH": sorted(map(list, setB)),
    "paper_promoted_nontext": sorted(map(list, paper12)),
    "overlap_A_with_paper": sorted(map(list, setA & paper12)),
    "overlap_B_with_paper": sorted(map(list, setB & paper12)),
    "phik_only_vs_paper": sorted(map(list, setB - paper12)),
    "paper_only_vs_phik": sorted(map(list, paper12 - setB)),
    "phik_values_for_paper_pairs": {f"{a}|{b}": phik_of[tuple(sorted((a, b)))]
                                    for a, b in paper13},
    "controls_in_setB": [list(p) for p in setB
                         if "rand_num" in p or "review_score_perm" in p],
}
with open("results/probe6_phik_crosscheck.json", "w") as fh:
    json.dump(out, fh, indent=2)
print(json.dumps(out, indent=2))
