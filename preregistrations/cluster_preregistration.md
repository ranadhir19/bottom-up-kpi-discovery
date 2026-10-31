# Paper 8 upgrade cycle 2 — cluster certificate: PRE-REGISTRATION

**Written before any cycle-2 run** (calibration, validation worlds, sensitivity).
Timestamp of writing: 2026-10-27 (America/New_York), prior to executing
`probe10_cluster_certificate.py`. All seeds, counts, designs and decision
rules below are fixed by this document. Any deviation discovered during the
run will be reported as a deviation, not silently absorbed.

## Question

Cycle 1 (probe 8) proved the family error certificate under i.i.d. rows and
showed it FAILING under clustered rows: 50 all-null worlds, n = 200,
20 clusters × 10 rows, ICC 0.5 — conjoint mean false 14.78 and P3-only mean
15.44 under a random row split, conjoint 8.94 and P3-only 20.72 under a
cluster-level split, against an i.i.d. certificate C^U = 0.849. Under the
cluster-level split the factorisation was restored (dependence factor 1.01)
but the i.i.d.-calibrated gammas were wrong under clustering.

Cycle 2 asks: does a **cluster-level certificate** — same theorem shape, with
per-half exceedance probabilities γ^cl calibrated **under the clustered
sampling at the half's cluster count**, and a split that assigns whole
clusters to halves — (a) factorise, and (b) dominate the observed false
counts on fresh clustered worlds?

## Designs (fixed)

**Clustered world generator — exactly the cycle-1 design**
(`probe8_certificate.clustered_columns`, reused by import, unmodified):
40 columns (24 numeric, 8 categorical-3, 8 categorical-2; 780 pairs:
276 num–num, 28 cat3×3, 28 cat2×2, 64 cat3×2, 192 num–cat g=3,
192 num–cat g=2); columns mutually independent; within a column, rows share
a cluster effect. Numeric: X = a_c + e, a ~ N(0, v), e ~ N(0, 1),
ICC = v/(v+1). Categorical: per-cluster probability vector
p_c ~ Dirichlet(3 × balanced), rows drawn from their cluster's p_c.

- **Primary design P:** 20 clusters × 10 rows (n = 200), v = 1 (ICC 0.5),
  Dirichlet total concentration 3 (as cycle 1). Cluster-split halves:
  10 clusters / 100 rows each.
- **Sensitivity design S:** 40 clusters × 5 rows (n = 200), same v = 1
  (ICC 0.5) and same Dirichlet concentration 3. Cluster-split halves:
  20 clusters / 100 rows each. (Same row count per half as P; only the
  cluster count / cluster size changes, isolating the K_h dependence the
  theorem predicts.) The alternative sensitivity (ICC 0.25 at 20 × 10) is
  NOT run in this cycle; design S is the pre-registered sensitivity.

**Promotion rule:** the manuscript's Definition 1 conjoint rule, unchanged:
(P1) Benjamini–Hochberg at q = 0.05 on full-sample p-values ∧ (P2)
full-sample corrected strength ≥ τ = 0.20 ∧ (P3) corrected strength
≥ τ = 0.20 on each half. Corrected strengths exactly as cycle 1 /
`facade_v2.correlate` (num–num |r|; cat–cat Bergsma V*; num–cat
max(ω², 0)). P3-only counts (pairs passing the split-half gate alone) are
recorded alongside, as in cycle 1, because E[#P3 passes] = Σ_j γ^cl_j1
γ^cl_j2 is the near-tight check of the factorisation.

## Calibration (fixed)

For each design (P, S) and each of the 6 branch cells, γ^cl is estimated
from **N = 50,000 independent null pairs per cell**. One calibration draw =
one half-sized clustered sample (P: 10 clusters × 10 rows; S: 20 clusters
× 5 rows) containing two mutually independent columns of the branch type,
each column with its own independent cluster effects drawn under the
design above; γ̂^cl = fraction of draws with corrected strength ≥ 0.20.

Justification of N = 50,000 (matches cycle 1): worst-case (γ = 0.5)
Monte Carlo SE is √(0.25/50000) ≈ 0.0022; at the expected primary num–num
γ^cl ≈ 0.26 the SE is ≈ 0.0020, so the simultaneous Clopper–Pearson upper
bound (below) exceeds γ̂ by ≲ 0.007 and the induced certificate
Σ n_cell (γ̂^U)² inflates the point prediction Σ n_cell γ̂² by only a few
percent — small against the world-mean SEs (≈ 0.5–0.7) that dominate the
validation comparison. For the smallest cells (num–cat, γ^cl ≈ 0.003),
50,000 draws give ≈ 150 expected exceedances, SE ≈ 0.00024 (≈ 8% relative).

Bounds: one-sided Clopper–Pearson. Per design, the certificate uses
**simultaneous** bounds at family level δ = 0.05 Bonferroni over that
design's K = 6 cells (per-cell level 1 − 0.05/6). Designs P and S are
separate certificates, each at simultaneous confidence 0.95; no joint
claim across the two designs is made. Point prediction C^cl,pt =
Σ_cell n_cell (γ̂^cl_cell)² and certificate C^cl,U = Σ_cell n_cell
(γ̂^cl,U_cell)², n_cell from the 780-pair family above.

## Validation worlds (fixed; fresh seeds, disjoint from cycle 1 and calibration)

- **Primary:** W_P = **60 worlds**, design P. Data seeds 930000 + w,
  cluster-split seeds 940000 + w, row-split seeds 950000 + w, w = 0…59.
  Cycle 1's cluster seeds were 860000/870000/880000 + w; no seed is reused.
  Each world is evaluated under (i) the **cluster-level split** (whole
  clusters to halves; the certified deployment) and (ii) descriptively,
  the random row split (NOT covered by the theorem; reported as a
  contrast diagnostic only — numerical non-violation there would not be
  claimed as certification).
- **Sensitivity:** W_S = **50 worlds**, design S. Data seeds 960000 + w,
  cluster-split seeds 970000 + w, w = 0…49. Cluster-level split only.

**Pre-registered checks (primary, cluster split):**
1. Factorisation: dependence factor = observed P3-only mean ÷ independence
   prediction from the worlds' own measured per-half marginals. Cycle 1
   found 1.01 under this split; the check is reported as measured, with
   ≈ 1.0 expected.
2. Point agreement: observed P3-only mean vs C^cl,pt, as a difference in
   combined standard errors (world-mean SE and delta-method calibration
   SE of C^cl,pt combined in quadrature). No pass/fail threshold is set;
   the z-value is reported.
3. Dominance: observed P3-only mean ≤ C^cl,U, and observed conjoint mean
   false count ≤ C^cl,U. (Pre-acknowledged: C^cl is expected ≫ 1 at this
   design, so the FWER/FDR certificate min(1, C^cl,U) = 1 is vacuous; the
   informative claim is E[V]. Observed FWER is still reported against it.)
4. Contrast: the cycle-1 i.i.d. certificate at n = 200 (C^U = 0.849, from
   `probe8_calibration.json`, simultaneous bounds) is recomputed from that
   file and reported on the same fresh worlds; the i.i.d. point prediction
   Σ γ̂² is also reported. Expected: both observed means exceed it by
   roughly an order of magnitude (the cycle-1 failure, reproduced on
   fresh seeds).
5. Within-world marginal check: per-cell per-half exceedance rates
   measured inside the validation worlds vs the calibration γ̂^cl
   (independent samples of the same quantity; agreement within sampling
   error expected).

Sensitivity checks: items 1–3 and 5 repeated for design S against the
design-S calibration/certificate. A failure in the sensitivity run is
reported as a failure; no re-tuning of designs, seeds, or counts after
seeing results.

## Statistic validation (fixed, run first)

Before calibration/worlds: (a) the cycle-1 scalar reference statistics
(imported from `probe8_certificate.py`) are checked against
`facade_v2.correlate` on **30 clustered null pairs** (10 num–num,
10 cat–cat across 3×3/2×2/3×2, 10 num–cat across g = 3/2; half-sized
primary-design samples; seed 91001), reporting max |Δstrength| and
max |Δp| (cycle 1 achieved 2.78×10⁻¹⁷ / 4.89×10⁻¹⁴ on i.i.d. pairs;
tolerances 10⁻⁸ / 10⁻⁶); (b) probe 10's batched calibration kernels are
checked against those scalar references on 2,000 clustered draws per
branch (seed 91002), reporting max |Δstrength| and exact agreement of
exceedance counts. The world engine is cycle 1's `eval_pairs`/`run_world`,
imported unmodified (already validated in cycle 1 at 2.22×10⁻¹⁶).

## Environment

Archived venv `hidden_files/paper08_work/phik_venv/bin/python`
(Python 3.12.3, numpy 2.5.4, scipy 1.18.1, pandas 3.1.0), as cycle 1.
`facade_v2` is imported with `np.math` restored to the stdlib `math`
module (numpy 2.5 removed the alias); `facade_v2.py` itself unmodified.

## Outputs

`probe10_calibration_cluster.json` (designs P+S calibration),
`probe10_worlds_cluster.json` (primary worlds, both splits),
`probe10_worlds_sensitivity.json` (design-S worlds),
`probe10_validation.json` (statistic validation), run log, and
`cluster_report.md`. Nothing in the manuscript is edited; nothing is
submitted anywhere.
