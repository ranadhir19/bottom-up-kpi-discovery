# Pre-registration — Paper 8 third-dataset replication run (UCI Bank Marketing)

**Written before any data file is downloaded or inspected**
(2026-10-27 15:04 EDT / 19:04 UTC, builder, Paper 8 upgrade pass).
This document was written from the public UCI dataset documentation
page only (UCI Machine Learning Repository, "Bank Marketing",
dataset 222, fetched as page text on 2026-10-27 before writing) and
from the fixed pipeline specification of the two earlier Paper 8
field runs. No Bank Marketing data file had been downloaded, opened,
or inspected at the time of writing. Frozen in this project record
before `probe9_bank_run.py` is written or run.

**Honesty standard (same as the earlier runs).** This is an
in-record plan-before-run with **no external registry** — there is
no registry timestamp; the discipline is that this file exists in
the project record before the run, exactly as for the Olist
pre-registration (`hidden_files/paper08_work/olist_preregistration.md`)
and the Telco pre-registration
(`hidden_files/paper08_work/telco_preregistration.md`), and any
deviation discovered later must be disclosed in the run report
rather than silently adapted.

## Data (fixed in advance)

- **Dataset:** UCI Bank Marketing (Moro, Cortez & Rita), direct
  marketing campaigns of a Portuguese banking institution.
- **Version (fixed):** `bank-additional-full.csv` from the official
  UCI Bank Marketing archive — the "bank-additional" release with
  **20 input attributes + the output attribute `y`**, per the UCI
  documentation. The older 17-attribute `bank-full.csv` /
  `bank.csv` files are **not** used.
- **Source (fixed):** the official UCI archive only
  (`archive.ics.uci.edu`, dataset 222). No Kaggle, GitHub, or other
  mirror copy may be substituted. If the official archive cannot be
  obtained, the run stops and that is reported; it does not fall
  back to a mirror.
- **Expected file shape (to be verified after download, before
  mining):** semicolon-delimited, one header row, **41,188 rows**,
  **21 columns** (20 inputs + `y`), with columns in the documented
  order:
  `age, job, marital, education, default, housing, loan, contact,
  month, day_of_week, duration, campaign, pdays, previous,
  poutcome, emp.var.rate, cons.price.idx, cons.conf.idx,
  euribor3m, nr.employed, y`.
  If the downloaded file does not match the documented row/column
  counts and column names, the run stops and the mismatch is
  reported as a deviation; the plan is not quietly re-fit to
  whatever file arrived.

## Documented schema (from the UCI documentation page alone)

Bank client data:
- `age` (numeric)
- `job` (categorical: "admin.", "blue-collar", "entrepreneur",
  "housemaid", "management", "retired", "self-employed",
  "services", "student", "technician", "unemployed", "unknown")
- `marital` (categorical: "divorced", "married", "single",
  "unknown"; "divorced" means divorced or widowed)
- `education` (categorical: "basic.4y", "basic.6y", "basic.9y",
  "high.school", "illiterate", "professional.course",
  "university.degree", "unknown")
- `default` (categorical: "no", "yes", "unknown") — has credit in
  default?
- `housing` (categorical: "no", "yes", "unknown") — has housing
  loan?
- `loan` (categorical: "no", "yes", "unknown") — has personal loan?

Last contact of the current campaign:
- `contact` (categorical: "cellular", "telephone")
- `month` (categorical: "jan", …, "dec")
- `day_of_week` (categorical: "mon", "tue", "wed", "thu", "fri")
- `duration` (numeric, seconds) — last contact duration. The UCI
  documentation's own note, quoted because a classification rule
  below rests on it: this attribute "highly affects the output
  target (e.g., if duration=0 then y='no'). Yet, the duration is
  not known before a call is performed. Also, after the end of the
  call y is obviously known. Thus, this input should only be
  included for benchmark purposes and should be discarded if the
  intention is to have a realistic predictive model."

Other attributes:
- `campaign` (numeric) — number of contacts performed during this
  campaign and for this client (includes last contact)
- `pdays` (numeric) — number of days that passed by after the
  client was last contacted from a previous campaign (999 means
  the client was not previously contacted)
- `previous` (numeric) — number of contacts performed before this
  campaign and for this client
- `poutcome` (categorical: "failure", "nonexistent", "success") —
  outcome of the previous marketing campaign

Social and economic context attributes:
- `emp.var.rate` (numeric) — employment variation rate, quarterly
  indicator
- `cons.price.idx` (numeric) — consumer price index, monthly
  indicator
- `cons.conf.idx` (numeric) — consumer confidence index, monthly
  indicator
- `euribor3m` (numeric) — Euribor 3-month rate, daily indicator
- `nr.employed` (numeric) — number of employees, quarterly
  indicator

Output:
- `y` (categorical: "no", "yes") — has the client subscribed a
  term deposit?

## Analysis table

All rows of `bank-additional-full.csv`, exactly as delivered. No
row filtering, no sampling, no recoding, no dropping of the
documented "unknown" levels (they are levels of the schema, not
missing data), and no transformation of any column (in particular
`pdays` is used as the documented numeric column, sentinel value
999 included, exactly as in the schema). Pairs are computed
pairwise-complete via the facade's per-pair `n`, as in the earlier
runs.

Type assignment (fixed from the documentation above):
- **numeric (10):** `age`, `duration`, `campaign`, `pdays`,
  `previous`, `emp.var.rate`, `cons.price.idx`, `cons.conf.idx`,
  `euribor3m`, `nr.employed`
- **categorical (11):** `job`, `marital`, `education`, `default`,
  `housing`, `loan`, `contact`, `month`, `day_of_week`,
  `poutcome`, `y`
- **`y` is included as an ordinary column** in the pair universe,
  exactly as `Churn` was in the Telco run; it receives no special
  treatment in mining (its known leakage partner, `duration`, is
  handled by the classification rules below, not by exclusion).

## Negative controls (appended before mining, fixed in advance)

Exactly two, in the same style as probes 5 and 7:
- `rand_num` — standard normal, `numpy.random.RandomState(20261027)`
  (the same seed as both earlier runs);
- `age_perm` — a permuted copy of the real column `age`,
  `numpy.random.RandomState(1)` (the same permutation seed as both
  earlier runs).
Both are typed numeric. No other synthetic column is added.

## Pair universe (fixed)

All unordered pairs of the resulting 21 + 2 = **23 columns**:
C(23,2) = **253 candidate pairs**. There are no text columns in
this dataset; there are no text pairs. The run must assert the
family size is exactly 253; any other count stops the run.

## Pipeline (fixed; identical machinery to probes 5 and 7)

Statistics are those of `facade_v2.py` (imported unchanged from
`hidden_files/paper08_work/probes/`): type-aware `correlate()`
with exact p-values; reportable (corrected) strength is Bergsma
V* for categorical×categorical, Hays ω² (floored at 0) for
numeric×categorical, and |Pearson r| for numeric×numeric — i.e.
`abs(strength_corrected)` where the facade provides one, else
`abs(strength)`, exactly as in probes 5/7. A pair is a *promoted
candidate KPI link* iff all of:
1. Benjamini–Hochberg rejection at **q = 0.05** across all 253
   p-values (`facade_v2.bh_reject`);
2. corrected strength ≥ **τ = 0.20**;
3. **split-half gate P3:** corrected strength ≥ 0.20 in both
   disjoint halves of a random split with **seed 7**
   (`numpy.random.RandomState(7).permutation`, first half /
   second half), exactly as in probes 5/7.
No other gate, filter, or post-hoc adjustment may be applied.

## Classification rules — NAME-BASED, written before any result

Applied to every promoted pair, **in this order**; the first rule
that fires determines the class. The rules use column names and
the documented schema only — they were written without seeing any
result of this run.

- **(0) Control:** the pair involves `rand_num` or `age_perm`.
  Any such promoted pair falsifies the run's error control (see
  expectations) and is reported as a control outcome, not as a
  finding.
- **(a) Leakage / mechanical — `duration`:** the pair involves
  `duration`. By the documentation quoted above, `duration` is
  only known after the call ends and strongly affects `y`; it is
  a post-outcome measurement, so any association involving it is
  leakage/mechanical by rule, including `duration`–`y` and
  `duration` paired with any other column.
- **(b) Campaign-history coupling:** **both** columns are in
  {`pdays`, `previous`, `poutcome`}. In the documented schema
  these three describe the same previous-campaign history and are
  coupled by construction: a client with no previous campaign
  (`previous` = 0) necessarily has `poutcome` = "nonexistent" and
  `pdays` = 999, so pairs among them share a deterministic
  "was previously contacted" structure. A pair with exactly one
  campaign-history column does **not** fire this rule and falls
  through.
- **(c) Shared-context (macro indicators):** the pair involves at
  least one of {`emp.var.rate`, `cons.price.idx`, `cons.conf.idx`,
  `euribor3m`, `nr.employed`}. These are national social/economic
  indicators that vary only over calendar time, not across clients
  contacted at the same time; any association involving one is
  attributed to shared temporal context by rule.
- **(d) Definitional / identity:** one column is a deterministic
  function of the other in the raw documented schema (none is
  known a priori in this file; if such a pair is promoted it is
  classified here with the function stated).
- **(e) Substantive:** everything else. In particular, pairs among
  the sociodemographic / client columns {`age`, `job`, `marital`,
  `education`, `default`, `housing`, `loan`} — and any other pair
  no earlier rule claims — are substantive by this rule.

## Pre-registered expectations

- Neither negative-control column appears in any promoted pair.
  If one does, the run's error control is falsified and that
  outcome is reported as such.
- No prediction is made about how many genuine pairs promote;
  zero is a reportable outcome.
- No prediction is made about the classification mix; the counts
  under the rules above are reported whatever they are.

## Reporting (fixed)

Row and column counts as verified against the documentation;
pair count; BH rejection count; BH ∧ floor count; split-half
removals (named pairs); promoted final count and full list with
methods, per-pair `n`, corrected strengths, and p-values; both
control columns' outcomes across the pipeline; classification
counts under the rules above, with every promoted pair's class;
the full 253-row pair table saved as CSV; dataset provenance
(official UCI URL, archive/file hash, verification result).
