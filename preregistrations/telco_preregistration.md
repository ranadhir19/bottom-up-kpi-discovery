# Pre-registration — Paper 8 strengthening replication run (IBM Telco)

**Written before the run is executed** (2026-10-27, builder, Paper 8
strengthening pass responding to Hypothesis Tester defect D8:
single field dataset). Frozen in this project record before
`probe7_telco_run.py` is run; the same in-record limitation as the
Olist pre-registration applies (no external registry timestamp —
recorded openly, per tester defect D5).

## Data
IBM Telco customer churn public sample (`telco_churn.csv`, as
archived in Paper 7's field-study folder): one row per customer,
21 columns including the `customerID` identifier (excluded) and
the `Churn` label.

## Analysis table
All rows of the file. `TotalCharges` is blank for a small number
of rows in the standard release; pairs are computed pairwise-complete
(the facade's per-pair `n_XY`), as in the Olist run.
- numeric: `tenure`, `MonthlyCharges`, `TotalCharges`;
- categorical: `gender`, `SeniorCitizen`, `Partner`, `Dependents`,
  `PhoneService`, `MultipleLines`, `InternetService`,
  `OnlineSecurity`, `OnlineBackup`, `DeviceProtection`,
  `TechSupport`, `StreamingTV`, `StreamingMovies`, `Contract`,
  `PaperlessBilling`, `PaymentMethod`, `Churn` (17 columns);
- negative controls appended before mining: `rand_num` (standard
  normal) and `monthly_perm` (`MonthlyCharges` permuted), seeds
  fixed in the run script.
Family: all pairs among 3 + 17 + 2 = 22 columns → C(22,2) = 231 pairs.
No text columns exist in this dataset; no text pairs.

## Promotion rule (identical to the Olist run)
1. Benjamini–Hochberg rejection at q = 0.05 across all 231 p-values;
2. corrected strength ≥ 0.20 (Bergsma V*, Hays ω², |Pearson r|);
3. split-half reproducibility: corrected strength ≥ 0.20 in both
   disjoint halves (random split, seed 7).

## Classification rules — written before any result is inspected
Applied to the promoted set, in this order:
- **(a) Accumulated-billing coupling (partly mechanical):** pairs
  involving `TotalCharges` together with `tenure` or
  `MonthlyCharges`. In this schema TotalCharges is the accumulated
  billed amount over the customer's tenure at roughly the monthly
  charge, so association is expected by construction.
- **(b) Structural (shared service level):** pairs between
  `InternetService` and any of the six internet add-on columns
  (`OnlineSecurity`, `OnlineBackup`, `DeviceProtection`,
  `TechSupport`, `StreamingTV`, `StreamingMovies`), and pairs
  among those six add-on columns: each add-on column takes the
  level "No internet service" exactly when `InternetService` = No,
  so the pairs share a deterministic level by construction.
- **(c) Definitional / identity:** a pair where one column is a
  deterministic function of the other in the raw schema (none is
  known a priori in this file; if one appears it is classified
  here with the function stated).
- **(d) Substantive:** everything else.

## Pre-registered expectations
- Neither negative-control column appears in any promoted pair;
  if one does, the run's error control is falsified and reported.
- No prediction is made about the number of promoted pairs; zero
  is a reportable outcome.

## Reporting
BH rejection count; promoted list with branch, n_XY, corrected
strengths, split-half outcome; control outcomes; classification
of every promoted pair under the rules above.
