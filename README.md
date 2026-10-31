# Bottom-Up KPI Discovery — replication package

Replication package for the paper:

> **A Type-Aware Bottom-Up Discovery Pipeline for Cross-Domain Enterprise KPIs
> with Spurious-Correlation Filtering**
> Ranadhir Ghosh, FIS Global, Jacksonville, Florida, USA
> (ranadhir.ghosh@fisglobal.com; ORCID 0009-0001-5439-4258)
> Submitted to *Knowledge and Information Systems* (Springer), 27 October 2026
> (submission ID 57e2ecf6-4cae-4a23-a146-4b35176d83d4).

The manuscript's Data and Code Availability statement promises that the
pipeline implementation, the probe scripts, the pre-registrations and the
field studies' full pair tables "are to be deposited in a public repository
with a persistent identifier before submission". This repository is that
deposit. It contains the pipeline implementation, every probe script behind
Sections 4–9 of the paper (including the certificate calibration and stress
test, the stability-selection comparison, and the third field run), the
pre-registrations, the full pair tables, and the deposited result files —
plus a validator that recomputes the headline numbers from the deposited
pair tables and run records and checks them against the printed values.

**Headline results** (all recomputed by `validate.py` from this deposit).
In pre-registered runs on three public warehouses — Olist e-commerce
(107 pairs), IBM Telco (231 pairs) and UCI Bank Marketing (253 pairs) —
the pipeline promotes 13, 88 and 33 pairs respectively and promotes no
negative-control pair on any of the three. A calibration audit shows the
strength-only rule of the original design averages 97.9 false discoveries
per null world at per-pair n = 60, while the conjoint rule (Benjamini–
Hochberg at q = 0.05, corrected-strength floor 0.20, split-half
reproducibility) averages 0.0. The family error certificate
(Theorem 3.6), evaluated on the three field families, is at most
4.2 × 10⁻⁶ (largest deposited value 4.14 × 10⁻⁶, Bank Marketing).

## Layout

| path | contents |
|---|---|
| `pipeline/facade_v2.py` | the pipeline implementation: the type-aware `correlate()` facade (Pearson r; ANOVA η²/ω²; Cramér's V with Bergsma's V* correction; exact p-values; BH-FDR; persistence test). numpy + pandas only |
| `probes/` | the original probe suite (calibration, FDR, persistence, split-half, text branch, the Olist and Telco field runs, the PhiK cross-check) |
| `probes/results/` | the deposited result JSONs and the Olist + Telco **full pair tables** |
| `upgrade/theorem/` | probe 8: the Theorem 3.6 family-error certificate — calibration (50,000 null pairs/cell), worlds, validation, report |
| `upgrade/cluster/` | probe 10: the cluster-level certificate extension, its pre-registration and report |
| `upgrade/fieldcert/` | probe 11: the certificate evaluated on the three field families (`field_certificates.json`) |
| `upgrade/thirddata/` | probe 9: the third field run (Bank Marketing) — script, run record, **full pair table**, report |
| `upgrade/baseline/` | the complementary-pairs stability-selection (CPSS) comparison: engine, protocol, synthetic + field runs, field pair tables |
| `preregistrations/` | the three field-run pre-registrations (Olist, Telco, Bank) + the cluster-extension pre-registration |
| `supplementary/bukpi_seed_results.md` | the seed-facade sanity run the probe suite was built on |
| `MANIFEST.md` | the deposit inventory: promised class → deposited files, exclusions, and the edits made to deposit copies |
| `validate.py` | the deposit validator (see below) |

## Validating the deposit

```sh
python validate.py
```

Standard library only. Exits 0 and prints `OVERALL PASS` only if all
25 checks pass: presence/parse of every artifact class (Stage A);
recomputation of the headline numbers from the deposited pair tables
and run JSONs — family sizes, BH-and-floor counts recomputed from the
tables' own columns, promoted sets, the 97.9 calibration figure, the
certificate bounds (Stage B); and near-miss negatives — the split-half
condition must be load-bearing (Telco 89 → 88, Bank 34 → 33) and the
planted random control column must never pass (Stage C).

## Regenerating the field runs

The datasets are public and are **not** redistributed here:

- **Olist** — Brazilian E-Commerce Public Dataset (Olist, 2018):
  the five CSVs `olist_orders_dataset.csv`, `olist_order_items_dataset.csv`,
  `olist_order_payments_dataset.csv`, `olist_order_reviews_dataset.csv`,
  `olist_customers_dataset.csv`.
- **Telco** — IBM Telco customer churn sample (`telco_churn.csv`, 7,043 rows).
- **Bank Marketing** — UCI dataset 222, official archive
  https://archive.ics.uci.edu/static/public/222/bank+marketing.zip
  (SHA-256 `e0bf5f5d…beb4`; inner `bank-additional.zip` SHA-256
  `a607b5ed…65cd`; `bank-additional-full.csv` SHA-256 `74adfc57…4d8`,
  all pinned in `upgrade/thirddata/probe9_bank_run.json`).

Put the Olist + Telco CSVs in one directory and point the field probes
at it via the `BUKPI_FIELD_DATA` environment variable (default
`data/field`); place `bank-additional-full.csv` beside probe 9. Then,
from the repository root:

```sh
python probes/probe5_olist_run.py     # -> probes/results/probe5_olist_run.json + full pair table
python probes/probe7_telco_run.py     # -> probes/results/probe7_telco_run.json + full pair table
python upgrade/thirddata/probe9_bank_run.py
python upgrade/theorem/probe8_certificate.py
python upgrade/cluster/probe10_cluster_certificate.py
python upgrade/fieldcert/probe11_field_certificates.py
```

Requirements: Python 3.12 with numpy and pandas (the calibration and
synthetic probes need nothing else). `probes/probe6_phik_crosscheck.py`
additionally needs the `phik` package — it is a cross-check of the
facade against an independent implementation, not a headline path.

## One disclosed defect in the preserved record

`probes/probe5b_promoted_signed.py` is deposited exactly as preserved,
and as preserved it is **truncated mid-statement** (an unclosed list at
its line 67) and does not compile. It was an auxiliary signed-strength
check on the Olist run, not a headline path; its deposited output
(`probes/results/probe5b_promoted_signed.json`, 13 entries — the signed
strengths of the 13 promoted pairs) is intact and is what the manuscript
record rests on. The file was not silently repaired: inventing the
missing lines would fabricate an artifact. Every other deposited script
compiles, and every headline number is recomputed by `validate.py`
from the pair tables and run records without probe 5b.

## License

Code and documentation: MIT (see `LICENSE`). The deposited pair tables
and result files derive from the public datasets cited above, under
their respective terms.
