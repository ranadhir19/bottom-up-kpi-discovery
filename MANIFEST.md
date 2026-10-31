# Deposit manifest — Paper 8 replication package

Deposit made 2026-10-31 to cure the availability-statement defect recorded
in the paper's Scoring System v2 Editor Screen backfill (S4 = 0: the
promised deposit had never been made). Source of the deposit: the paper's
preserved supplementary record in the author's programme workspace
(`paper08_work/` and `paper08_upgrade/`, preserved 2026-10-27).

## Promised vs deposited

The manuscript's Data and Code Availability statement promises four
artifact classes. All four were located in the preserved record and are
deposited:

| Promised class | Found in the preserved record | Deposited as |
|---|---|---|
| Pipeline implementation | `paper08_work/probes/facade_v2.py` — the production type-aware `correlate()` facade | `pipeline/facade_v2.py` |
| Probe scripts (calibration–field sections) | probes 1–7 + 1b/2b/3b/5b/6 in `paper08_work/probes/`, with their result JSONs | `probes/`, `probes/results/` |
| — incl. certificate calibration and stress test | probe 8 (theorem), probe 10 (cluster extension), probe 11 (field certificates) in `paper08_upgrade/` | `upgrade/theorem/`, `upgrade/cluster/`, `upgrade/fieldcert/` |
| — incl. stability-selection comparison | CPSS engine + protocol + synthetic/field runs in `paper08_upgrade/baseline/` | `upgrade/baseline/` |
| — incl. the third field run | probe 9 (Bank Marketing) in `paper08_upgrade/thirddata/` | `upgrade/thirddata/` |
| The three pre-registrations | `olist_preregistration.md`, `telco_preregistration.md` (work dir), `bank_preregistration.md` (thirddata) | `preregistrations/` (+ the cluster extension's pre-registration as a fourth, labelled as such) |
| Field studies' full pair tables | `probe5_full_pair_table.csv` (Olist, 109 rows), `probe7_full_pair_table.csv` (Telco, 231), `probe9_full_pair_table.csv` (Bank, 253) | `probes/results/`, `upgrade/thirddata/` (CPSS field pair tables in `upgrade/baseline/`) |

## Known gap (disclosed, not repaired)

`probes/probe5b_promoted_signed.py` is preserved **truncated mid-statement**
(unclosed list, its line 67) and does not compile — in the archival
original, not only in this copy. Its output
(`probes/results/probe5b_promoted_signed.json`, 13 entries) is intact.
The missing lines were not reconstructed: that would fabricate an
artifact. Probe 5b is an auxiliary check, not a headline path; all
headline numbers validate from the pair tables and run records.

## Deposit edits (the only differences from the archival originals)

Nine scripts carried absolute paths into the author's internal workspace.
In the deposit copies — and only in the deposit copies — those constants
were made configurable, each marked inline with a `DEPOSIT COPY` comment:

- `probes/probe4_text.py`, `probe5_olist_run.py`, `probe5b_promoted_signed.py`,
  `probe6_phik_crosscheck.py`, `probe7_telco_run.py`,
  `upgrade/baseline/run_field.py`: the field-data directory now comes from
  the `BUKPI_FIELD_DATA` environment variable (default `data/field`).
- `upgrade/thirddata/probe9_bank_run.py`: the probes path and its working
  directory are now resolved relative to the script's own location.
- `upgrade/baseline/run_synthetic.py`, `upgrade/baseline/validate_engine.py`:
  the probes base path is now resolved relative to the script's location.

No computational code, parameter, seed or result was altered. The
archival originals in the programme workspace are untouched.

## Exclusions (deliberate)

- **Raw datasets** (Olist/Telco CSVs; the UCI Bank archives): public and
  cited in the manuscript; not redistributed. Sources and the Bank
  archive's pinned SHA-256 hashes are in the README.
- `paper08_work/phik_venv/` (a 404 MB working virtualenv) and all
  `__pycache__` directories: environment, not artifacts.
- `paper08_work/results/bukpi/charts.png`: a binary chart from the seed
  run; its results table is deposited as
  `supplementary/bukpi_seed_results.md`.
- Internal build/upgrade reports and hypothesis-test reports: programme
  bookkeeping, not manuscript artifacts. The per-cycle technical reports
  that document the deposited probes (`theorem_report.md`,
  `cluster_report.md`, `fieldcert_report.md`, `thirddata_report.md`,
  `baseline_report.md`, `cpss_protocol.md`) ARE deposited with their probes.

## Verification at deposit time

- `python3 validate.py` → **OVERALL PASS** (25 checks; headline numbers
  recomputed from the deposited pair tables and run JSONs).
- All deposited scripts compile under Python 3.12 except the disclosed
  probe 5b truncation.
- A pre-push content scan of the exact deposited tree was run and
  cleared before publication (the only candidate hits, judged in
  context, are the author's own manuscript byline above); the scan
  record is held with the paper's submission record.
