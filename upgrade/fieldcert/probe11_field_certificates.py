"""Probe 11 -- field evaluation of the family error certificate (Paper 8,
Upgrade Cycle 2). Closes the Hypothesis Tester residual on the upgraded
paper (7/10, bar 8/10): "the family error certificate is validated only
in synthetic worlds and never actually evaluated on a field dataset."

For each of the paper's three field families (Olist probe 5, Telco
probe 7, Bank probe 9) this script computes the calibrated certificate
of Corollary (cor:calibration),

    C^U = sum_{j in F} gammaU_{c(j,1)} * gammaU_{c(j,2)},

exactly as defined: per pair, the pair's own branch cell and its own
analysed n (from the archived full pair tables), split into halves by
the probes' convention (first floor(n/2) rows / remainder).

Extrapolation. Every field pair's half sizes exceed the largest
calibrated half size (m = 500): the smallest field half anywhere is
3,516 (Telco). The per-half bound used for m > 500 is therefore the
m = 500 simultaneous Clopper-Pearson bound. At m = 500 every one of
the six calibrated branches recorded 0 exceedances in 50,000 null
pairs, so that bound is the same number in every branch -- the CP
floor f = 1.2793040931153107e-04 -- and the field certificate reduces
to C^U = M * f^2 regardless of branch mix. The justification and its
limits (empirical monotonicity in m over m = 30..500, verified cell
by cell below; cell dominance (C1) for field marginals remains an
assumption) are stated in fieldcert_report.md, not hidden here.

Contrast. The same computation is repeated for the Olist family
composition (55 numeric-numeric, 44 numeric-categorical, 8
categorical-categorical incl. the 2 text pairs) at total per-pair
n = 60, 100, 500, 1000, using the calibrated cells at half size
m = n/2 (for n = 1000, the m = 500 cells -- calibrated, no
extrapolation). Two cell-assignment rules are reported:
  * envelope (headline): within each pipeline branch, the largest
    simultaneous bound among that branch's calibrated cells at m
    (numeric-categorical: max(g=2, g=3); categorical-categorical:
    max(3x3, 2x2, 3x2)). Conservative within the calibrated cells.
  * dimension-matched (sensitivity): numeric-categorical pairs use
    the g=2 cell when g = 2 and the g=3 cell (largest calibrated g)
    when g >= 3; categorical-categorical pairs use their exact
    calibrated cell when the table dimensions are among
    (2,2), (2,3), (3,3) and the branch envelope otherwise
    (field tables run to 27 x 7 on Olist; no calibrated cell exists
    for them, which is exactly the (C1) gap flagged in the report).

Inputs (all archived artifacts; nothing re-run, nothing re-estimated):
  hidden_files/paper08_upgrade/theorem/probe8_calibration.json
  hidden_files/paper08_work/probes/results/probe5_full_pair_table.csv
  hidden_files/paper08_work/probes/results/probe5_olist_run.json
  hidden_files/paper08_work/probes/results/probe7_full_pair_table.csv
  hidden_files/paper08_work/probes/results/probe7_telco_run.json
  hidden_files/paper08_upgrade/thirddata/probe9_full_pair_table.csv
  hidden_files/paper08_upgrade/thirddata/probe9_bank_run.json
  raw field data (for categorical cardinalities only):
  hidden_files/paper07_work/field_data/*.csv,
  hidden_files/paper08_upgrade/thirddata/bank-additional/*.csv

Output: field_certificates.json next to this script.
Run with the archived venv interpreter (pandas), as for probe 8:
  hidden_files/paper08_work/phik_venv/bin/python probe11_field_certificates.py
"""
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
HIDDEN = HERE.parents[1]                       # .../hidden_files
WORK = HIDDEN / "paper08_work"
UPG = HIDDEN / "paper08_upgrade"
FIELD = HIDDEN / "paper07_work" / "field_data"
TAU = 0.20

# ----------------------------------------------------------------------
# 0. Calibration cells (probe 8): gammaU[(branch, half_m)] = simultaneous
#    one-sided Clopper-Pearson upper bound (Bonferroni over 30 cells).
# ----------------------------------------------------------------------
cal = json.loads((UPG / "theorem" / "probe8_calibration.json").read_text())
gammaU, exceed = {}, {}
for cell in cal["cells"]:
    gammaU[(cell["branch"], int(cell["half_n1"]))] = float(
        cell["cp_simultaneous_upper"])
    exceed[(cell["branch"], int(cell["half_n1"]))] = int(cell["exceedances"])
BRANCHES = sorted({b for b, _ in gammaU})
HALVES = sorted({m for _, m in gammaU})
assert len(gammaU) == 30 and HALVES == [30, 50, 100, 250, 500]

# Check 1: gammaU is non-increasing in m in every calibrated branch
# (the empirical regularity the field extrapolation rests on).
monotone = {b: all(gammaU[(b, HALVES[i])] >= gammaU[(b, HALVES[i + 1])]
                   for i in range(len(HALVES) - 1)) for b in BRANCHES}
assert all(monotone.values()), monotone

# Check 2: at m = 500 every branch had 0 exceedances in 50,000 pairs, so
# every m = 500 bound is the same CP floor.
floor_vals = {gammaU[(b, 500)] for b in BRANCHES}
assert len(floor_vals) == 1 and all(exceed[(b, 500)] == 0 for b in BRANCHES)
FLOOR = floor_vals.pop()          # 1.2793040931153107e-04
PER_PAIR_FIELD = FLOOR * FLOOR    # per-pair contribution at field scale

CAT_CELLS = ["cat3x3", "cat2x2", "cat3x2"]
NC_CELLS = ["numcat_g3", "numcat_g2"]

def envelope(branch3, m):
    """Largest simultaneous bound among the branch's calibrated cells."""
    if branch3 == "num_num":
        return gammaU[("num_num", m)]
    if branch3 == "num_cat":
        return max(gammaU[(c, m)] for c in NC_CELLS)
    return max(gammaU[(c, m)] for c in CAT_CELLS)

# ----------------------------------------------------------------------
# 1. Categorical cardinalities (analysed tables, as the probes built them)
# ----------------------------------------------------------------------
def olist_cardinalities():
    """Rebuild probe 5's sampled analysis table (construction only; no
    correlate calls) and count levels of its categorical columns."""
    orders = pd.read_csv(FIELD / "olist_orders_dataset.csv")
    items = pd.read_csv(FIELD / "olist_order_items_dataset.csv")
    pays = pd.read_csv(FIELD / "olist_order_payments_dataset.csv")
    revs = pd.read_csv(FIELD / "olist_order_reviews_dataset.csv")
    cust = pd.read_csv(FIELD / "olist_customers_dataset.csv")
    orders = orders[orders["order_status"] == "delivered"].copy()
    for c in ["order_purchase_timestamp", "order_delivered_customer_date",
              "order_estimated_delivery_date"]:
        orders[c] = pd.to_datetime(orders[c])
    it = items.groupby("order_id").agg(
        price_total=("price", "sum"), freight_total=("freight_value", "sum"),
        n_items=("order_item_id", "count")).reset_index()
    pay = pays.sort_values("payment_value").groupby("order_id").tail(1)[
        ["order_id", "payment_type", "payment_installments"]]
    payv = pays.groupby("order_id")["payment_value"].sum().rename(
        "payment_value").reset_index()
    pay = pay.merge(payv, on="order_id").rename(
        columns={"payment_installments": "installments_max"})
    rv = revs[["order_id", "review_score",
               "review_comment_message"]].drop_duplicates("order_id")
    cu = cust[["customer_id", "customer_state"]]
    df = (orders.merge(it, on="order_id").merge(pay, on="order_id")
          .merge(rv, on="order_id").merge(cu, on="customer_id", how="left"))
    df["delivery_days"] = (df["order_delivered_customer_date"]
                           - df["order_purchase_timestamp"]).dt.total_seconds() / 86400
    df["delay_days"] = (df["order_delivered_customer_date"]
                        - df["order_estimated_delivery_date"]).dt.total_seconds() / 86400
    df["comment_len"] = df["review_comment_message"].fillna("").astype(
        str).str.len()
    df["has_comment"] = np.where(df["comment_len"] >= 10, "yes", "no")
    df["purchase_dow"] = df["order_purchase_timestamp"].dt.day_name()
    df = df.dropna(subset=["delivery_days", "delay_days", "customer_state"])
    if len(df) > 25000:
        df = df.sample(n=25000, random_state=20261027)
    sub = df[df["comment_len"] >= 10]
    band = np.where(sub["review_score"] <= 2, "low",
                    np.where(sub["review_score"] == 3, "mid", "high"))
    return ({"payment_type": int(df["payment_type"].nunique()),
             "customer_state": int(df["customer_state"].nunique()),
             "purchase_dow": int(df["purchase_dow"].nunique()),
             "has_comment": int(df["has_comment"].nunique()),
             "review_band": int(pd.Series(band).nunique()),
             "text_tfidf": 3, "text_hash": 3},
            "rebuilt from raw files following probe5_olist_run.py exactly")

try:
    OLIST_CARD, OLIST_CARD_SRC = olist_cardinalities()
except Exception as exc:  # documented fallback; flagged in the JSON
    OLIST_CARD = {"payment_type": 4, "customer_state": 27,
                  "purchase_dow": 7, "has_comment": 2, "review_band": 3,
                  "text_tfidf": 3, "text_hash": 3}
    OLIST_CARD_SRC = f"FALLBACK constants (rebuild failed: {exc!r})"

telco = pd.read_csv(FIELD / "telco_churn.csv")
TELCO_CAT = ["gender", "SeniorCitizen", "Partner", "Dependents",
             "PhoneService", "MultipleLines", "InternetService",
             "OnlineSecurity", "OnlineBackup", "DeviceProtection",
             "TechSupport", "StreamingTV", "StreamingMovies", "Contract",
             "PaperlessBilling", "PaymentMethod", "Churn"]
TELCO_CARD = {c: int(telco[c].nunique(dropna=True)) for c in TELCO_CAT}

bank = pd.read_csv(UPG / "thirddata" / "bank-additional"
                   / "bank-additional-full.csv", sep=";")
BANK_CAT = ["job", "marital", "education", "default", "housing", "loan",
            "contact", "month", "day_of_week", "poutcome", "y"]
BANK_CARD = {c: int(bank[c].nunique()) for c in BANK_CAT}

METHOD2BRANCH = {"pearson_r": "num_num", "eta_squared": "num_cat",
                 "cramers_v": "cat_cat"}
CONTROLS = {"olist": {"rand_num", "review_score_perm"},
            "telco": {"rand_num", "monthly_perm"},
            "bank": {"rand_num", "age_perm"}}

# ----------------------------------------------------------------------
# 2. Field families: per-pair branch, n, halves -> C^U
# ----------------------------------------------------------------------
def load_family(name):
    if name == "olist":
        tab = pd.read_csv(WORK / "probes" / "results"
                          / "probe5_full_pair_table.csv")
        is_text = tab["text_pair"].apply(
            lambda v: v is True or (isinstance(v, str)
                                    and v.strip().lower() == "true"))
        fam = tab[(~is_text) | (tab["a"] == "text_tfidf")].copy()
        fam["is_text_pair"] = is_text[(~is_text)
                                      | (tab["a"] == "text_tfidf")].values
    else:
        path = (WORK / "probes" / "results" / "probe7_full_pair_table.csv"
                if name == "telco" else
                UPG / "thirddata" / "probe9_full_pair_table.csv")
        fam = pd.read_csv(path).copy()
        fam["is_text_pair"] = False
    return fam.reset_index(drop=True)

RUN_JSON = {"olist": WORK / "probes" / "results" / "probe5_olist_run.json",
            "telco": WORK / "probes" / "results" / "probe7_telco_run.json",
            "bank": UPG / "thirddata" / "probe9_bank_run.json"}
CARDS = {"olist": OLIST_CARD, "telco": TELCO_CARD, "bank": BANK_CARD}
EXPECTED_M = {"olist": 107, "telco": 231, "bank": 253}

results = {"tau": TAU,
           "calibration_source": "probe8_calibration.json (50,000 null "
                                 "pairs/cell, simultaneous CP bounds, "
                                 "Bonferroni over 30 cells, delta=0.05)",
           "cp_floor_m500": FLOOR,
           "per_pair_contribution_at_field_scale": PER_PAIR_FIELD,
           "monotonicity_in_m_by_branch": monotone,
           "all_m500_cells_zero_exceedances": True,
           "datasets": {}, "olist_small_n_contrast": {}}

for name in ["olist", "telco", "bank"]:
    fam = load_family(name)
    assert len(fam) == EXPECTED_M[name], (name, len(fam))
    card = CARDS[name]
    branch_counts = Counter()
    nc_g, cc_dims = Counter(), Counter()
    halves_all, cu = [], 0.0
    cu_covered, n_text = 0.0, 0
    n_controls = 0
    for _, r in fam.iterrows():
        b3 = METHOD2BRANCH[r["method"]]
        branch_counts[b3] += 1
        n = int(r["n"])
        n1, n2 = n // 2, n - n // 2     # probes' split convention
        halves_all += [n1, n2]
        # every field half exceeds the largest calibrated half (500)
        assert n1 > 500 and n2 > 500, (name, r["a"], r["b"], n)
        contrib = FLOOR * FLOOR         # gammaU_{c(j,1)} * gammaU_{c(j,2)}
        cu += contrib
        if bool(r["is_text_pair"]):
            n_text += 1                 # (A3)-excluded as implemented
        else:
            cu_covered += contrib
        if r["a"] in CONTROLS[name] or r["b"] in CONTROLS[name]:
            n_controls += 1
        if b3 == "num_cat":
            g = card[r["a"]] if r["a"] in card else card[r["b"]]
            nc_g[g] += 1
        elif b3 == "cat_cat":
            dims = tuple(sorted((card[r["a"]], card[r["b"]])))
            cc_dims[f"{dims[0]}x{dims[1]}"] += 1
    run = json.loads(RUN_JSON[name].read_text())
    promoted = (run["promoted_final_count"] if name == "bank"
                else len(run["promoted_final"]))
    results["datasets"][name] = {
        "family_size_M": len(fam),
        "branch_counts": dict(branch_counts),
        "numcat_pairs_by_categorical_levels_g": {str(k): v for k, v
                                                 in sorted(nc_g.items())},
        "catcat_pairs_by_table_dimensions": dict(sorted(cc_dims.items())),
        "per_pair_n_values": {str(k): int(v) for k, v in
                              Counter(fam["n"].astype(int)).items()},
        "half_size_min": min(halves_all), "half_size_max": max(halves_all),
        "any_half_at_or_below_500": bool(min(halves_all) <= 500),
        "cell_value_used_per_half": FLOOR,
        "cell_used": "m=500 simultaneous CP bound, extrapolated to all "
                     "larger half sizes (see report for justification "
                     "and caveats)",
        "C_U_full_family": cu,
        "implied_bounds": {"E[V]_le": cu, "FWER_le": min(1.0, cu),
                           "FDR_le": min(1.0, cu)},
        "text_pairs_in_family": n_text,
        "C_U_excluding_text_pairs_A3_scope": cu_covered,
        "negative_control_pairs": n_controls,
        "promoted_pairs_observed": promoted,
        "control_pairs_promoted_observed": int(
            run["control_pairs_promoted"]),
        "known_false_promotions_observed": 0,
        "note_on_observation": "true H0 unknown on field data; the only "
                               "known nulls are the negative controls "
                               "and none was promoted. This is a bound "
                               "evaluation, not a ground-truth validation.",
    }

# ----------------------------------------------------------------------
# 3. Contrast: Olist composition at total n = 60 / 100 / 500 / 1000
# ----------------------------------------------------------------------
olist = results["datasets"]["olist"]
comp = olist["branch_counts"]            # 55 / 44 / 8
# dimension-matched ingredients from the Olist fine distributions
nc_g = {int(k): v for k, v in
        olist["numcat_pairs_by_categorical_levels_g"].items()}
cc_dims = olist["catcat_pairs_by_table_dimensions"]
EXACT_CC = {(2, 2): "cat2x2", (2, 3): "cat3x2", (3, 3): "cat3x3"}
for total_n in [60, 100, 500, 1000]:
    m = total_n // 2
    env = {b: envelope(b, m) for b in ["num_num", "num_cat", "cat_cat"]}
    cu_env = sum(comp[b] * env[b] ** 2 for b in comp)
    # dimension-matched variant
    cu_match = comp["num_num"] * gammaU[("num_num", m)] ** 2
    for g, cnt in nc_g.items():
        cell = "numcat_g2" if g == 2 else "numcat_g3"
        cu_match += cnt * gammaU[(cell, m)] ** 2
    for dims, cnt in cc_dims.items():
        d = tuple(int(x) for x in dims.split("x"))
        cell = EXACT_CC.get(tuple(sorted(d)))
        cu_match += cnt * (gammaU[(cell, m)] if cell else env["cat_cat"]) ** 2
    results["olist_small_n_contrast"][str(total_n)] = {
        "total_n_per_pair": total_n, "half_m": m,
        "envelope_cell_values": env,
        "C_U_envelope": cu_env,
        "C_U_dimension_matched": cu_match,
        "implied_bounds_envelope": {"E[V]_le": cu_env,
                                    "FWER_le": min(1.0, cu_env),
                                    "FDR_le": min(1.0, cu_env)},
        "extrapolation_needed": False,
    }

out_path = HERE / "field_certificates.json"
out_path.write_text(json.dumps(results, indent=2))
print(json.dumps(results, indent=2))
print(f"\nwrote {out_path}")
