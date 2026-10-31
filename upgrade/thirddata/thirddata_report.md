# Paper 8 — Third-dataset replication report (UCI Bank Marketing, probe 9)

Builder report for the Paper 8 upgrade pass ("A Type-Aware Bottom-Up
Discovery Pipeline for Cross-Domain Enterprise KPIs with
Spurious-Correlation Filtering", target: *Knowledge and Information
Systems*). This is the third real-dataset run of the pipeline, after
Olist (probe 5: 107 pairs → 13 promoted, 0 controls) and IBM Telco
(probe 7: 231 pairs → 88 promoted, 0 controls). The manuscript .tex
was not edited; nothing was submitted.

Every number below is taken from `probe9_bank_run.json` (the run's
own output) or `probe9_full_pair_table.csv` (the full 253-row pair
table), both in this folder.

## Order of operations and timestamps (UTC; EDT = UTC−4)

| Step | Timestamp |
|---|---|
| Pre-registration `bank_preregistration.md` written (file mtime) | 2026-10-27 19:04:54 UTC (15:04 EDT) |
| Official UCI archive download started | 2026-10-27 19:05:15 UTC |
| `probe9_bank_run.py` run — started | 2026-10-27 19:06:37.741536 UTC |
| Run finished (single run, exit 0) | 2026-10-27 19:06:47.517979 UTC |

The pre-registration was written **before** any Bank Marketing data
file was downloaded or inspected, from the public UCI documentation
page only (UCI Machine Learning Repository, "Bank Marketing",
dataset 222, fetched as page text the same day). It is an in-record
plan-before-run with no external registry — the same honesty
standard, stated in the pre-registration itself, as the Olist and
Telco runs.

## Dataset provenance and verification

- Source page: `https://archive.ics.uci.edu/dataset/222/bank+marketing`
- Archive downloaded (official UCI only, no mirror):
  `https://archive.ics.uci.edu/static/public/222/bank+marketing.zip`
  — SHA-256 `e0bf5f5de5b846e2f18e9d90606637267d46dfa260e0f17bb12e605db5efbeb4`
- The outer archive contains `bank-additional.zip`
  — SHA-256 `a607b5edab6c6c75ce09c39142a77702c38123bd5aa7ae89a63503bbe17d65cd`
- File used: `bank-additional/bank-additional-full.csv`
  — SHA-256 `74adfc578bf77a7ff4bb1ba4a9f8709d9e3c6907342959c2c8416847e0afb4d8`
- Verification against the documentation, before mining:
  semicolon-delimited, **41,188 rows × 21 columns**, column names
  and order exactly as pre-registered (`age` … `nr.employed`, `y`),
  zero missing cells. Matches the documented bank-additional-full
  shape (20 input attributes + `y`). **Verification: PASS.**

## Deviations from the pre-registration

**None.** `deviations_from_preregistration` in the results JSON is
the empty list. The run implemented the pre-registered plan exactly:
all 41,188 rows, no sampling, no filtering, no recoding; `y` as an
ordinary column; exactly the two pre-named controls (`rand_num`,
`RandomState(20261027)`; `age_perm`, permutation of `age`,
`RandomState(1)`); facade_v2 statistics imported unchanged from
`hidden_files/paper08_work/probes/`; BH q = 0.05; corrected-strength
floor τ = 0.20; split-half gate P3 with seed 7. The nested-archive
layout (bank-additional.zip inside bank+marketing.zip) is the
official UCI distribution format, not a deviation.

## Headline counts (from `probe9_bank_run.json`)

- Analysis table: **41,188 rows**, 21 documented columns + 2
  controls = **23 columns**
- Pair universe: **253 pairs** — C(23,2), as pre-registered
- BH rejections at q = 0.05: **179**
- BH ∧ corrected-strength floor (≥ 0.20): **34**
- Split-half removals: **1** — `previous` × `cons.price.idx`
  (full-sample corrected strength 0.2031; failed the ≥ 0.20
  requirement in at least one half)
- **Promoted final: 33**
- **Control pairs promoted: 0**

## Negative-control outcome

43 pairs involve a control column. Four were BH-rejected, all at
negligible corrected strengths — `pdays` × `rand_num` (0.0124) and
three `rand_num` categorical pairs (≤ 0.00021) — and **none** passed
the strength floor; **0 control pairs were promoted**. The
pre-registered expectation (no control in any promoted pair) held;
the run's error control was not falsified. Across the three field
runs the control record is now 0 promoted controls on Olist, Telco,
and Bank Marketing.

## Classification of the 33 promoted pairs (pre-written rules)

| Class | Count |
|---|---|
| (0) control | 0 |
| (a) leakage / mechanical — `duration` | 0 |
| (b) campaign-history coupling | 3 |
| (c) shared-context (macro indicators) | 22 |
| (d) definitional / identity | 0 |
| (e) substantive | 8 |
| **Total** | **33** |

Notes on classes that fired zero times: rule (a) had no promoted
pair to classify — the strongest `duration` pair in the full table
is `duration` × `y` at corrected strength (Hays ω²) 0.1642,
BH-rejected but below the 0.20 floor, and every other `duration`
pair is ≤ 0.0717. The pre-registered leakage rule was checked
against the promoted set and correctly found nothing; on this
run the strength floor, not the rule, kept `duration` out. Rule
(d) fired zero times, as expected (no deterministic pair was known
a priori and none was promoted unclaimed by an earlier rule).

**Class (b) — campaign-history coupling (3):**
`pdays` × `previous` (Pearson r, 0.5875), `pdays` × `poutcome`
(ω², 0.9067), `previous` × `poutcome` (ω², 0.7888) — the
deterministic "was previously contacted" structure the
pre-registration named (`previous` = 0 ⇒ `poutcome` =
"nonexistent", `pdays` = 999).

**Class (c) — shared-context macro (22):** 7 of the 10 pairs among
the five macro indicators themselves (`emp.var.rate` ×
`euribor3m`, r = 0.9722; `euribor3m` × `nr.employed`, r = 0.9452;
`emp.var.rate` × `nr.employed`, r = 0.9070; `emp.var.rate` ×
`cons.price.idx`, r = 0.7753; `cons.price.idx` × `euribor3m`,
r = 0.6882; `cons.price.idx` × `nr.employed`, r = 0.5220;
`cons.conf.idx` × `euribor3m`, r = 0.2777 — the three macro–macro
pairs not promoted all involve `cons.conf.idx` against another
macro indicator), plus 15 pairs of a macro indicator with
`pdays`, `previous`, `month`, `contact`, or `poutcome` (e.g.
`previous` × `nr.employed`, r = 0.5013; `nr.employed` × `month`,
ω² = 0.4202). Full per-pair values are in the JSON promoted list.

**Class (e) — substantive (8):**
`age` × `job` (ω², 0.2546), `job` × `education` (V*, 0.3595),
`housing` × `loan` (V*, 0.7079), `contact` × `month` (V*, 0.6091),
`contact` × `poutcome` (V*, 0.2424), `month` × `poutcome`
(V*, 0.2424), `month` × `y` (V*, 0.2740), `poutcome` × `y`
(V*, 0.3204). All promoted pairs have n = 41,188 and p ≈ 0
(p = 0.0 at the facade's numerical precision).


## Honest assessment — what this third run adds

This run adds a third, genuinely different domain — Portuguese
bank telemarketing campaign records, with a mix of client
sociodemographics, campaign-history fields, a known post-call
leakage variable, and dated national macro indicators — and the
pipeline's error control held for a third time (0 of 43
control-involving pairs promoted; the pre-registered expectation
met on all three datasets). Its most useful result is arguably a
cautionary one for the paper's own framing: of the 33 promoted
pairs, only 8 are substantive under rules written before the run,
while 22 are shared-context macro co-movement (the five national
indicators move together over calendar time and drag
period-correlated columns with them) and 3 are the deterministic
campaign-history coupling — so on this dataset a raw promoted
count would overstate discovered KPI structure roughly four-fold,
and the name-based classification layer, not the promotion gates
alone, is what makes the output interpretable. Two limits should
be stated plainly: the pre-registered `duration` leakage rule
never fired because no `duration` pair cleared the strength floor
(`duration` × `y` reached ω² = 0.1642), so this run demonstrates
the rule was in place and checked, not that it caught a promoted
leakage pair; and, as with the earlier runs, the pre-registration
is in-record only, with no external registry timestamp.

## Files in this folder

- `bank_preregistration.md` — the frozen plan-before-run
- `probe9_bank_run.py` — the run script (imports `facade_v2` from
  `hidden_files/paper08_work/probes/`)
- `probe9_bank_run.json` — full results (counts, promoted list with
  corrected strengths and classes, control outcomes, hashes)
- `probe9_full_pair_table.csv` — all 253 pairs with method,
  strengths, p-values, BH and floor flags
- `thirddata_report.md` — this report
- `bank+marketing.zip`, `bank-additional.zip`,
  `bank-additional/bank-additional-full.csv` — the downloaded
  official UCI archive and the extracted data file (hashes above)
- `probe9_run_stdout.log`, `probe9_run_stderr.log` — run logs
  (stderr empty; exit 0)
