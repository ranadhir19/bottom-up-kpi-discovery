# CPSS head-to-head baseline — report

Paper 8 (KAIS submission): complementary-pairs stability selection
(CPSS; Shah & Samworth 2013, after Meinshausen & Bühlmann 2010) run
head-to-head against the manuscript's promotion rule
(P1∧P2∧P3 = BH q = 0.05 ∧ corrected-strength floor τ = 0.20 ∧
split-half), under the protocol frozen in `cpss_protocol.md` before
any run. Every number below traces to a JSON/CSV in this folder.

## Artefacts

| File | Contents |
|---|---|
| `cpss_protocol.md` | Pre-registered design (written first) |
| `cpss_engine.py` | Optimised numpy branch-statistics engine |
| `validate_engine.py`, `engine_validation.json` | Engine validation gate |
| `run_synthetic.py`, `synthetic_results.json`, `synthetic_run.log` | Synthetic worlds |
| `run_field.py`, `field_results.json`, `field_run.log` | Olist + Telco |
| `field_olist_pair_table.csv`, `field_telco_pair_table.csv` | Per-pair full-sample strength, p, BH, stability frequency, per-rule membership |

Implementation note: all runs use the system `python3` (numpy 1.26.4,
scipy 1.11.4, pandas 2.1.4). The `phik_venv` interpreter (numpy
2.5.4) cannot run `facade_v2.py` at all — numpy 2.x removed the
`np.math` attribute the facade uses — and the facade had to run
unmodified, since the paper rule is defined by it.

## 1. Engine validation (gate, protocol §7)

36 random pairs (12 per branch), each compared engine-vs-facade on
the full sample **and** on a random half-sample: 72 evaluations.

| Branch | Max abs difference in corrected strength |
|---|---|
| numeric–numeric | 0.0 |
| categorical–categorical | 1.11e-16 |
| numeric–categorical | 1.94e-16 |
| **Overall** | **1.94e-16** |

Gate (< 1e-6) passed by ten orders of magnitude
(`engine_validation.json`; runtime 0.87 s in-process).
Determinism spot-check: re-running synthetic planted world 0 in
isolation reproduced its per-rule counts and stability frequencies
exactly (max frequency difference 0.0).

## 2. Synthetic comparison (`synthetic_results.json`)

Probe-2 designs verbatim, n = 500, 10 null + 10 planted worlds,
CPSS B = 50 complementary pairs per world (protocol-declared
reduction from the field B = 100). One correction to the record:
probe 2's docstring says the planted world has "46 columns / 1029
other pairs", but its code as written produces **52 columns / 1326
pairs** (40 base + 12 plant columns); the verbatim code was run and
the actual counts are reported here. Null worlds: 40 columns /
780 pairs, as documented.

### Null worlds (all 780 pairs null in every world)

| Rule | Mean discoveries | Mean false discoveries | FWER |
|---|---|---|---|
| Paper conjoint | 0.0 | 0.0 | 0.00 |
| CPSS π = 0.90 | 0.0 | 0.0 | 0.00 |
| CPSS π = 0.80 | 0.0 | 0.0 | 0.00 |
| Hybrid (π = 0.90 ∧ BH) | 0.0 | 0.0 | 0.00 |

Separation margin: the largest stability frequency attained by any
null pair in any null world was **0.38** — less than half the
π = 0.80 threshold. Across planted worlds, 0 of 13,200 non-plant
pair-worlds reached frequency 0.80 (max 0.38).

### Planted worlds (6 plants per world)

| Rule | Mean discoveries | Mean false disc. | Power | Mean realised FDP |
|---|---|---|---|---|
| Paper conjoint | 4.4 | 0.0 | **0.733** | 0.000 |
| CPSS π = 0.90 | 3.8 | 0.0 | 0.633 | 0.000 |
| CPSS π = 0.80 | 4.3 | 0.0 | 0.717 | 0.000 |
| Hybrid (π = 0.90 ∧ BH) | 3.8 | 0.0 | 0.633 | 0.000 |

Per-plant mean stability frequencies (10 worlds): ρ = .50 → 1.000;
ρ = .35 → 1.000; V = .45 → 1.000; η² = .25 → 0.822 (≥ 0.90 in only
3/10 worlds); V = .25 → 0.764 (≥ 0.90 in 5/10); η² = .10 → 0.000
(below the τ floor by construction — no rule can promote it, as the
manuscript itself states). CPSS's power loss is concentrated exactly
in the two moderate plants whose half-sample (n = 250) corrected
strength fluctuates around τ. The hybrid is identical to CPSS
π = 0.90 in every world: BH never removes a pair that stability
selected, i.e. on these designs the stability frequency is the
binding constraint, not the p-value screen.

Runtime: 946.9 s total for 20 worlds (mean 35.8 s per null world,
58.9 s per planted world), dominated by the 100 half-sample
evaluations of up to 1326 pairs per world.

## 3. Field comparison (`field_results.json`)

Tables rebuilt exactly as probes 5/7; B = 100 complementary pairs.
Sanity anchor: the recomputed paper rule reproduces the stored
probe results **exactly** — 13/13 Olist pairs, 88/88 Telco pairs
(`paper_set_matches_probe_reference: true` for both datasets).

### Olist (25,000 orders; family = 107 pairs)

| Rule | Promoted | Overlap with paper (13) | Control pairs promoted |
|---|---|---|---|
| Paper conjoint | 13 | 13 | 0 of 27 |
| CPSS π = 0.90 | 12 | 12 | 0 of 27 |
| CPSS π = 0.80 | 12 | 12 | 0 of 27 |
| Hybrid | 12 | 12 | 0 of 27 |

No rule promotes any pair the paper's rule does not. All CPSS
variants miss exactly one paper pair: **freight_total —
installments_max** (full-sample corrected strength 0.2037,
stability frequency 0.680). Control pairs: 27 in the family; the
maximum stability frequency among them is **0.000** — no control
pair reached τ in any of the 200 halves.

### Telco (7,043 rows; family = 231 pairs)

| Rule | Promoted | Overlap with paper (88) | Control pairs promoted |
|---|---|---|---|
| Paper conjoint | 88 | 88 | 0 of 41 |
| CPSS π = 0.90 | 86 | 86 | 0 of 41 |
| CPSS π = 0.80 | 89 | 88 | 0 of 41 |
| Hybrid | 86 | 86 | 0 of 41 |

CPSS π = 0.90 (and the hybrid) miss two paper pairs, both
near-floor: **InternetService — Contract** (strength 0.2063,
frequency 0.815) and **SeniorCitizen — OnlineSecurity** (strength
0.2104, frequency 0.870). CPSS π = 0.80 recovers both and adds one
pair the paper's rule does not promote: **SeniorCitizen —
Dependents** (full-sample strength 0.2109, p = 2.8e-70, stability
frequency 0.885). This is precisely the pair probe 7's split-half
gate removed (`promoted_but_splithalf_failed` in
`probe7_telco_run.json`): it passes BH ∧ floor, fails the single
seed-7 split, yet clears τ in 88.5% of 200 CPSS halves. Control
pairs: 41 in the family; maximum control stability frequency
**0.000**; no rule promotes a control pair.

Field runtimes: Olist total 44.7 s, of which the CPSS loop is
38.0 s; Telco total 28.8 s, CPSS loop 17.9 s.

## 4. Honest assessment — where CPSS wins, where it loses

**Where CPSS beats the paper's rule.**
(1) *It replaces a one-draw lottery with a measured quantity.* The
paper's P3 is a single fixed split (seed 7): for near-floor pairs
its verdict is an artefact of that partition. The Telco
SeniorCitizen — Dependents pair is the concrete case — rejected by
the paper's rule on one split while clearing the floor in 88.5% of
200 halves. CPSS exposes exactly this gradation (frequencies 0.68 /
0.815 / 0.87 / 0.885 for the four contested field pairs) instead of
a binary coin flip. (2) *Error behaviour is at least as clean on
every design run*: zero discoveries in 10/10 null worlds, zero false
discoveries in planted worlds at both thresholds, zero control
promotions on both field datasets, with a wide empirical margin
(maximum null frequency 0.38 against thresholds of 0.80/0.90).
(3) *The hybrid result is informative*: BH adds nothing on top of
CPSS π = 0.90 anywhere, so a stability-frequency rule can stand
alone without the p-value machinery and its dependence caveats.

**Where CPSS loses.**
(1) *Power at the primary threshold.* Synthetic power 0.633 vs the
paper's 0.733 at π = 0.90 (0.717 at the π = 0.80 sensitivity, still
below the paper). The mechanism is structural: CPSS judges a pair
by its strength at n/2, where sampling noise is √2 larger, so
genuinely moderate associations (η² = .25, V = .25 plants) sit
below π = 0.90 in most worlds, while the paper's rule judges
strength at full n and asks the halves only for a single
reproducibility check. (2) *Recall at the floor boundary on real
data*: CPSS π = 0.90 drops 1 of 13 Olist and 2 of 88 Telco paper
pairs, all with full-sample strengths 0.204–0.210. One can defend
this ("barely above the floor at full n" ≈ "not stable"), but it is
a real recall cost, and the π = 0.80 sensitivity buys most of it
back only by also promoting a pair (SeniorCitizen — Dependents)
whose status is genuinely ambiguous. (3) *Compute.* CPSS evaluates
the full family 200 times: 38.0 s vs 6.7 s for everything else on
Olist, 17.9 s vs 10.9 s on Telco, and 946.9 s for the 20 synthetic
worlds — roughly two orders of magnitude more pair evaluations than
the conjoint rule, which needs one full pass plus two half-sample
evaluations per BH ∧ floor survivor (13 and 89 survivors here). At
warehouse scale that gap, not accuracy, is the practical argument
for the paper's cheaper gate. (4) *No transferable guarantee*: the
Shah–Samworth error bound does not apply to this adaptation
(protocol §1), so CPSS's cleanliness here is measured, not proved —
the same epistemic status the manuscript claims for its own rule.

**Bottom line for a referee response**: the published-procedure
head-to-head does not overturn the paper's rule. CPSS at its
primary threshold is more conservative (equal error control, lower
power, ~200× the resampling cost); at π = 0.80 it nearly matches
the paper's power and disagrees with it on exactly the near-floor
pairs where any threshold rule is arbitrary — and in the one field
case where they disagree at π = 0.80 (SeniorCitizen — Dependents),
the stability frequency arguably indicts the paper's single-split
gate rather than vindicating CPSS. The defensible claim is
complementarity: CPSS frequencies are a strictly richer diagnostic
to *report alongside* promoted pairs, at a compute price the
conjoint rule avoids.
