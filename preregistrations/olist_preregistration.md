# Pre-registration — Paper 8 real-data discovery run (Olist)

**Written before the run is executed** (2026-10-27, builder, Paper 8
production). The run is the end-to-end evidence for the revised
manuscript: the v2 pipeline mining one real warehouse-style dataset.
Anything not specified here that the run later reports will be labelled
exploratory.

## Data
Olist Brazilian E-Commerce public dataset, CSVs as archived in this
project's Paper 7 field-study folder (`olist_orders_dataset.csv`,
`olist_order_items_dataset.csv`, `olist_order_payments_dataset.csv`,
`olist_order_reviews_dataset.csv`, `olist_customers_dataset.csv`).
The products file is not in the archive, so no product-category column
exists in this run; that absence is recorded, not worked around.

## Analysis table
One row per order with `order_status = delivered`, at least one item,
at least one payment record, and a review score. Columns:
- numeric: `price_total`, `freight_total`, `payment_value`, `n_items`,
  `installments_max`, `review_score`, `delivery_days`
  (delivered − purchase, days), `delay_days` (delivered − estimated,
  days), `comment_len` (characters, 0 when no comment);
- categorical: `payment_type` (type of the largest payment row),
  `customer_state`, `purchase_dow`, `has_comment`;
- text: `review_comment_message` (pairs computed on rows whose comment
  has ≥ 10 characters);
- negative controls: `rand_num` (standard normal, seed 20261027) and
  `review_score_perm` (review_score permuted, seed 1).
If the table exceeds 25,000 rows it is sampled to 25,000 rows with
seed 20261027 before mining.

## Mining protocol
All pairs among the 9 numeric + 4 categorical + 2 control columns
(C(15,2) = 105 pairs) plus 2 text pairs (text × review band derived
from review_score as categorical {1–2, 3, 4–5}; text × payment_type)
= 107 candidate pairs. Type dispatch per the v2 facade; corrected
strengths (Bergsma V*, Hays ω², Pearson r) are the reportable
strengths; exact p-values.

**Text reduction (primary): word TF–IDF + Lloyd k-means (k = 3).**
Chosen because the accompanying text-branch benchmark (probe 4,
including a task on this dataset's comments) showed TF–IDF dominating
the shipped 64-d hash sketch on real text (V = 0.621 vs 0.330).
Hash-64 results for the same two pairs are reported alongside.

## Promotion rule (fixed in advance)
A pair is a *promoted candidate KPI link* iff all of:
1. Benjamini–Hochberg rejection at q = 0.05 across all 107 p-values;
2. corrected strength ≥ 0.20;
3. split-half reproducibility: corrected strength ≥ 0.20 in both
   disjoint halves (random split, seed 7);
4. not flagged high spurious-risk by the facade heuristics
   (|strength| > 0.95 with n < 30, or a hint substring — neither can
   fire at this sample size / with these column names; recorded as
   checked).
The original persistence share (20 half-sample windows, half-strength,
60% rule) is reported descriptively for promoted pairs; it is **not**
a gate (probe 3 showed it passes 73–100% of null pairs).

## Pre-registered expectations
- Neither negative-control column appears in any promoted pair. If one
  does, the run's error control is falsified and that outcome is
  reported as such.
- No prediction is made about how many genuine pairs promote; zero is
  a reportable outcome.

## Reporting
Full promoted list with strengths and p-values; the 15 strongest pairs
by corrected strength regardless of promotion; control outcomes; text
pairs under both reductions; the analysis-table size and the sampling
actually applied.
