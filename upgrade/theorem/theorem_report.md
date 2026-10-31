# Paper 8 upgrade — error-control theorem and simulation validation

**Paper:** "A Type-Aware Bottom-Up Discovery Pipeline for Cross-Domain Enterprise KPIs with Spurious-Correlation Filtering" → *Knowledge and Information Systems*.
**Role of this document:** builder output for the theorem track of the upgrade (the item intended to lift the contribution past the R2 cap "assembly-plus-audit, no new method/theorem"). It contains (1) a family error-control theorem for the conjoint promotion rule with a complete proof, in LaTeX-ready form for the manuscript builder; (2) a Monte Carlo calibration of the per-half exceedance probabilities the theorem needs; (3) world-level validation of the resulting certificate, planted-world behaviour, and a cluster stress test. No manuscript file was edited, nothing was submitted, authorship was not touched. All work is in this folder (`hidden_files/paper08_upgrade/theorem/`).

**Verdict on the tasked statement, up front:** the target statement is **correct as stated** — no fatal flaw was found. Four refinements are needed for it to be rigorous, and they are built into the theorem below: (i) the two half-tail probabilities γ_{j1}, γ_{j2} must be allowed to differ (unequal half sizes or shapes) — the product form already permits this; (ii) the FWER/FDR conclusions sharpen to min(1, C); (iii) the factorisation needs the half statistics to be functions of their own half's rows only — no full-sample-fitted component (this excludes the pipeline's text pairs as currently implemented; see Remark R3); (iv) the calibration corollary transfers to field data only under an explicit cell-dominance assumption (C1) — balanced-marginal calibration cells do not automatically dominate unbalanced field marginals. The simulation validation is exact on point (iv) because the simulated worlds use the calibrated marginals; the field transfer remains an assumption, stated as such.

---

## 1. The theorem (LaTeX-ready)

The block below is written to drop into the manuscript (it uses the manuscript's notation: Definition 1 promotion rule with P1 Benjamini–Hochberg at q = 0.05, P2 corrected-strength floor τ = 0.20, P3 split-half reproducibility; S_{jh} is pair j's corrected strength computed on half h). It generalises the manuscript's Proposition 1 from a single pair to the whole family.

```latex
\begin{assumption}[Sampling and split]
\label{ass:theorem}
Fix a pair family $\mathcal{F}$ of size $M$, a floor $\tau > 0$, and
the promotion rule of Definition~\ref{def:promotion}. For pair
$j \in \mathcal{F}$ let $n_j$ be the number of analysed rows and
$S_{j1}, S_{j2}$ the corrected strengths computed on the two halves
of a disjoint partition of those rows, of fixed sizes $n_{j1}$ and
$n_{j2}$, $n_{j1} + n_{j2} = n_j$.
\begin{enumerate}
\item[(A1)] The analysed rows are independent and identically
distributed draws from a row distribution; the null $H_{0j}$ states
that pair $j$'s two columns are independent under that distribution.
\item[(A2)] The partition into halves is fixed in advance, or drawn
independently of the data.
\item[(A3)] Each $S_{jh}$ is a function of half $h$'s rows only: no
quantity fitted on the full sample or on the other half (no shared
preprocessing, clustering, imputation, or normalisation) enters it.
\end{enumerate}
Write $\gamma_{jh} = \Pr_{H_{0j}}[S_{jh} \ge \tau]$, pair $j$'s own
per-half exceedance probability at half size $n_{jh}$.
\end{assumption}

\begin{theorem}[Family error control from the split-half gate alone]
\label{thm:family}
Under Assumption~\ref{ass:theorem}, for every null pair $j$,
\[
\Pr[\text{$j$ is promoted}] \;\le\; \gamma_{j1}\,\gamma_{j2}.
\]
Consequently, with $\mathcal{H}_0 \subseteq \mathcal{F}$ the set of
true null pairs, $V$ the number of promoted null pairs, $R$ the
total number promoted, and
$C = \sum_{j \in \mathcal{H}_0} \gamma_{j1}\gamma_{j2}$,
\[
\mathbb{E}[V] \;\le\; C, \qquad
\mathrm{FWER} = \Pr[V \ge 1] \;\le\; \min(1, C), \qquad
\mathrm{FDR} = \mathbb{E}\!\left[\tfrac{V}{\max(R,1)}\right]
\;\le\; \min(1, C).
\]
No assumption about dependence across pairs is made, and no
property of\/ {\rm (P1)} or {\rm (P2)} is used: the bounds hold for
any promotion rule that has {\rm (P3)} as a necessary condition.
\end{theorem}

\begin{proof}
Fix a null pair $j$. Promotion requires (P3), so
$\{\text{$j$ promoted}\} \subseteq \{S_{j1} \ge \tau,\,
S_{j2} \ge \tau\}$ and
$\Pr[\text{$j$ promoted}] \le \Pr[S_{j1} \ge \tau,\,
S_{j2} \ge \tau]$.
By (A1) the rows are independent; by (A2) the two halves are
disjoint sets of rows whose composition does not depend on the
data; by (A3) $S_{j1}$ and $S_{j2}$ are functions of disjoint
collections of independent rows, hence independent random
variables. The joint event therefore factorises:
$\Pr[S_{j1} \ge \tau,\, S_{j2} \ge \tau]
= \Pr[S_{j1} \ge \tau]\,\Pr[S_{j2} \ge \tau]
= \gamma_{j1}\gamma_{j2}$.
Linearity of expectation over $j \in \mathcal{H}_0$ — which needs
no independence across pairs — gives
$\mathbb{E}[V] = \sum_{j \in \mathcal{H}_0} \Pr[\text{$j$ promoted}]
\le C$.
For the family-wise rate, the union bound (equivalently Markov's
inequality applied to $V$) gives
$\Pr[V \ge 1] = \Pr[\bigcup_{j \in \mathcal{H}_0}
\{\text{$j$ promoted}\}] \le \sum_{j \in \mathcal{H}_0}
\Pr[\text{$j$ promoted}] \le C$, and a probability is at most $1$.
For the false discovery rate, the false discovery proportion
satisfies $V/\max(R,1) \le \mathbf{1}\{V > 0\}$ pointwise: if
$V = 0$ both sides are $0$, and if $V \ge 1$ then $R \ge V \ge 1$
so the proportion is at most $1$. Taking expectations,
$\mathrm{FDR} \le \Pr[V \ge 1] = \mathrm{FWER} \le \min(1, C)$.
\end{proof}

\begin{corollary}[Calibrated certificate]
\label{cor:calibration}
Partition the (branch, half-size) combinations that occur in
$\mathcal{F}$ into $K$ cells $c$, each with a true per-half
exceedance $\gamma_c$ under a specified calibration null
distribution. Draw $N_c$ independent null pairs from cell $c$'s
calibration distribution, independently of the field sample and of
the other cells, and let $\hat\gamma^{U}_{c}$ be the one-sided
Clopper--Pearson upper confidence bound for $\gamma_c$ at level
$1 - \delta/K$. Assume in addition:
\begin{enumerate}
\item[(C1)] (cell dominance) every field null pair $j$ and half $h$
satisfies $\gamma_{jh} \le \gamma_{c(j,h)}$, where $c(j,h)$ is the
pair's cell;
\item[(C2)] the calibration draws are independent of the field
sample.
\end{enumerate}
Then with probability at least $1 - \delta$ over the calibration
draws, the certificate
$C^{U} = \sum_{j \in \mathcal{F}} \hat\gamma^{U}_{c(j,1)}
\hat\gamma^{U}_{c(j,2)}$ (the sum may be taken over all pairs,
since $\mathcal{H}_0$ is unknown) satisfies: conditional on the
calibration, $\mathbb{E}[V] \le C^{U}$,
$\mathrm{FWER} \le \min(1, C^{U})$ and
$\mathrm{FDR} \le \min(1, C^{U})$ for the field sample.
\end{corollary}

\begin{proof}
Each Clopper--Pearson bound covers its cell parameter with
probability at least $1 - \delta/K$; by the union bound over the
$K$ cells, the event
$E = \{\gamma_c \le \hat\gamma^{U}_{c} \text{ for all } c\}$ has
probability at least $1 - \delta$ over the calibration draws.
Condition on a calibration outcome in $E$. The values
$\hat\gamma^{U}_{c}$ are then fixed numbers, the field sample is
independent of them by (C2), and by (C1) and $E$,
$\gamma_{j1}\gamma_{j2} \le
\hat\gamma^{U}_{c(j,1)}\hat\gamma^{U}_{c(j,2)}$ for every null
pair. Theorem~\ref{thm:family} applied to the field sample gives
the three bounds with $C \le C^{U}$.
\end{proof}
```

**Remarks (for the manuscript's discussion, and load-bearing for §4 below).**

- **R1 (the bound is a P3-only bound, and is conservative exactly where P1∧P2 work).** The proof never uses BH or the full-sample floor. The certificate therefore upper-bounds the *split-half gate's own* false-pass count, and the conjoint rule can only do better. Section 3 quantifies the slack: at n = 60 the certificate is 28.84 while the conjoint rule's measured mean false count is 0.00 — the gap is the work of P1∧P2, which this theorem deliberately does not credit.
- **R2 (relation to the manuscript's Proposition 1).** With equal half sizes and a common marginal bound γ_τ ≥ max(γ_{j1}, γ_{j2}), the per-pair statement reduces to the manuscript's γ_τ². The theorem is the family-level generalisation with pair-specific γs.
- **R3 (where the assumptions can fail).** (a) *Clustered rows* violate (A1): §4 shows the certificate failing by a factor of ≈17 in clustered worlds, decomposed into marginal inflation and cross-half dependence. (b) *Data-dependent splits* violate (A2). (c) *Full-sample-fitted statistics* violate (A3): the pipeline's **text pairs** are outside the theorem's scope as implemented, because the text reduction (clustering) is fitted on the full commented subset before the split, so the two halves' label statistics share the fitted reduction; the theorem covers text pairs only if the reduction is fitted on separate data or refitted within each half. (d) *Missingness:* the theorem is about the analysed rows; if pairwise complete-case selection depends on the pair's own values in a way that distorts the null, (A1) can fail for the analysed subsample. (e) The statement is superpopulation: probabilities are over i.i.d. row draws, not over re-splits of one fixed table.
- **R4 (vacuity at small n).** When C > 1 the FWER/FDR parts are vacuous and only E[V] ≤ C carries information. At n = 60, C ≈ 28.8 for a 780-pair family: the split-half gate alone is *not* an error-control device at that size; the conjoint rule's small-n protection comes from P1∧P2, as the manuscript's own component tests already showed.

---

## 2. Simulation setup and statistic validation

All simulations are in `probe8_certificate.py`, run with the archived venv interpreter (`hidden_files/paper08_work/phik_venv/bin/python`; Python 3.12.3, numpy 2.5.4, pandas 3.1.0, scipy 1.18.1). Corrected strengths match `probes/facade_v2.py` exactly: numeric–numeric |r|; categorical–categorical Bergsma V\*; numeric–categorical max(ω², 0). p-values: exact t pivot, χ², and F tails via scipy.

**Validation of the fast statistics (file `probe8_validation.json`).** Before any run, the reference implementations were checked against `facade_v2.correlate` on 30 random pairs (10 numeric–numeric, 10 categorical–categorical across 3×3 / 2×2 / 3×2, 10 numeric–categorical across g = 3 / 2; seed 81001): maximum absolute difference **2.78×10⁻¹⁷** in strength and **4.89×10⁻¹⁴** in p-value. The batched world engine was then checked against the reference implementations on all 780 pairs of a random world (seed 81002): maximum differences 2.22×10⁻¹⁶ (strength) and 7.77×10⁻¹⁶ (p). (Import note: the archived venv's numpy 2.5 removed the `np.math` alias that `facade_v2` uses; the harness restores it as the stdlib `math` module — which is what the alias was — before importing. `facade_v2.py` itself is unmodified.)

**Designs.** Null and planted worlds follow `probe2_fdr.py`: base worlds have 40 columns (24 numeric Gaussian, 8 categorical with 3 balanced levels, 8 with 2; 780 pairs: 276 numeric–numeric, 28 cat-3×3, 28 cat-2×2, 64 cat-3×2, 192 numeric–cat g = 3, 192 numeric–cat g = 2). The conjoint rule is applied exactly as in Definition 1 (BH q = 0.05 on full-sample p-values ∧ full-sample corrected strength ≥ 0.20 ∧ split-half ≥ 0.20 on each half), with **one fixed row permutation per world** shared by all pairs (split seeds recorded in the JSONs), matching the field runs' fixed-seed split. P3-only counts (pairs passing the split-half gate alone) are recorded alongside, because they are the quantity the theorem's proof actually bounds before P1∧P2 are applied.

---

## 3. Calibration: per-half exceedance γ (file `probe8_calibration.json`)

50,000 independent null pairs per cell (the task floor was 20,000), seed per cell recorded in the JSON. Cells are (branch × total n); totals are even, so both halves have size n/2 and one estimate serves both halves. γ̂ is the empirical exceedance rate of the corrected strength at τ = 0.20 on a half-size sample; γ̂_U^CP95 is the one-sided Clopper–Pearson 95% upper bound; γ̂_U^sim is the simultaneous bound used for certificates (per-cell level 1 − 0.05/30, Bonferroni over all 30 cells, simultaneous confidence 0.95).

| Half m (total n) | Branch | k / 50,000 | γ̂ | γ̂_U^CP95 | γ̂_U^sim |
|---|---|---|---|---|---|
| 30 (60) | num–num | 14,434 | 0.28868 | 0.29203 | 0.29466 |
| 30 (60) | cat 3×3 | 9,106 | 0.18212 | 0.18498 | 0.18724 |
| 30 (60) | cat 2×2 | 7,242 | 0.14484 | 0.14745 | 0.14951 |
| 30 (60) | cat 3×2 | 10,866 | 0.21732 | 0.22038 | 0.22278 |
| 30 (60) | num–cat g=3 | 894 | 0.01788 | 0.01889 | 0.01969 |
| 30 (60) | num–cat g=2 | 348 | 0.00696 | 0.00760 | 0.00812 |
| 50 (100) | num–num | 8,207 | 0.16414 | 0.16689 | 0.16905 |
| 50 (100) | cat 3×3 | 4,794 | 0.09588 | 0.09807 | 0.09980 |
| 50 (100) | cat 2×2 | 4,300 | 0.08600 | 0.08809 | 0.08974 |
| 50 (100) | cat 3×2 | 6,863 | 0.13726 | 0.13982 | 0.14183 |
| 50 (100) | num–cat g=3 | 88 | 0.00176 | 0.00210 | 0.00238 |
| 50 (100) | num–cat g=2 | 29 | 0.00058 | 0.00079 | 0.00097 |
| 100 (200) | num–num | 2,263 | 0.04526 | 0.04682 | 0.04805 |
| 100 (200) | cat 3×3 | 900 | 0.01800 | 0.01901 | 0.01981 |
| 100 (200) | cat 2×2 | 1,212 | 0.02424 | 0.02540 | 0.02633 |
| 100 (200) | cat 3×2 | 2,514 | 0.05028 | 0.05192 | 0.05321 |
| 100 (200) | num–cat g=3 | 0 | 0 | 5.99×10⁻⁵ | 1.28×10⁻⁴ |
| 100 (200) | num–cat g=2 | 0 | 0 | 5.99×10⁻⁵ | 1.28×10⁻⁴ |
| 250 (500) | num–num | 90 | 0.00180 | 0.00214 | 0.00243 |
| 250 (500) | cat 3×3 | 5 | 0.00010 | 0.00021 | 0.00031 |
| 250 (500) | cat 2×2 | 44 | 0.00088 | 0.00113 | 0.00134 |
| 250 (500) | cat 3×2 | 125 | 0.00250 | 0.00290 | 0.00323 |
| 250 (500) | num–cat g=3 | 0 | 0 | 5.99×10⁻⁵ | 1.28×10⁻⁴ |
| 250 (500) | num–cat g=2 | 0 | 0 | 5.99×10⁻⁵ | 1.28×10⁻⁴ |
| 500 (1000) | all six branches | 0 | 0 | 5.99×10⁻⁵ | 1.28×10⁻⁴ |

Notes. (1) The numeric–categorical γ collapses fastest: at half size 100 the ω² null exceedance at τ = 0.20 was 0 in 50,000 pairs for both g (the CP bound, 1.28×10⁻⁴, is what enters certificates — conservative by construction where γ̂ = 0). (2) At half size 500 every cell had k = 0, so all certificates at total n = 1000 would rest on the CP bound 1.28×10⁻⁴ per half; no n = 1000 worlds were run (the tasked world sizes are 60/100/500), so this row is reported as calibration only. (3) The point-estimate certificates Σγ̂² are also computed in §4 for comparison: the theorem's expectation is Σγ² with the *true* γ, so the point version is the sharper prediction and the simultaneous version is the certificate proper.

---

## 4. World validation (files `probe8_worlds_null.json`, `probe8_worlds_planted.json`)

### 4.1 All-null worlds — the certificate holds for the conjoint rule at every n

50 worlds per n (data seeds 820000 + 100n + w, split seeds 830000 + 100n + w, w = 0…49). "Certificate" C^U = Σ over the 780 pairs of (γ̂_U^sim,cell)² from §3; C^pt is the same sum with point estimates γ̂. In all-null worlds FDR = FWER.

| n | Mean false V̄ (SE) | FWER (= FDR) | P3-only mean (SE) | C^pt = Σγ̂² | Certificate C^U | V̄ ≤ C^U? | P3 mean vs C^U |
|---|---|---|---|---|---|---|---|
| 60 | 0.00 (0.00) | 0.00 | 26.98 (0.61) | 27.61 | 28.84 | yes | 26.98 ≤ 28.84 ✓ (ratio 1.07) |
| 100 | 0.02 (0.02) | 0.02 | 10.16 (0.44) | 9.11 | 9.68 | yes (484× slack) | 10.16 vs 9.68 — nominally above by 0.48 = 1.1 SE |
| 500 | 0.00 (0.00) | 0.00 | 0.00 (0.00) | 0.00132 | 0.00236 | yes | 0.00 ≤ 0.00236 ✓ |

Reading, plainly:

- **The conjoint rule's false count respects the certificate at all three sizes**, with large slack: zero false promotions in all 50 worlds at n = 60 and n = 500, and exactly one false promotion in one world at n = 100 (mean 0.02 against a certificate of 9.68). This is the expected picture — the certificate ignores P1∧P2 (Remark R1), and at these sizes BH∧floor do nearly all of the conjoint rule's work.
- **The P3-only diagnostic is the near-tight check of the factorisation itself**, since E[#P3 passes] = Σ_j γ_{j1}γ_{j2} exactly (with true γ). At n = 60 the observed 26.98 sits between the point prediction 27.61 (within 1.1 SE) and the simultaneous certificate 28.84: the product bound is essentially exact there, i.e. the two halves really are independent and the calibration really does measure γ. At n = 500 no pair passed P3 in any world, consistent with a per-world expectation of ≈0.001–0.002.
- **The one blemish, reported rather than smoothed over:** at n = 100 the observed P3-only mean (10.16, SE 0.44) is nominally *above* the simultaneous certificate (9.68) by 0.48, i.e. 1.1 standard errors of the world mean, and above the point prediction (9.11) by 2.4 SE. The theorem bounds an expectation; the observed mean is a Monte Carlo estimate of it (50 worlds, with column-sharing overdispersion already inside the SE), and the calibration side also carries Monte Carlo error, so a ≈1 SE excursion above the simultaneous bound is within the noise the design tolerates — but on a strict reading, in this run the simultaneous certificate did not dominate the observed P3 mean at n = 100. The conjoint rule itself was nowhere near the bound (0.02). A tester re-running with different seeds should expect the P3 diagnostic to land within roughly ±1.5 SE of C^pt.

### 4.2 Planted worlds — FDP 0, power 0.75, certificate over null pairs respected

20 worlds at n = 500 of the probe-2 planted design as written in `probe2_fdr.py` (see the side finding in §6: the code builds 52 columns / 1,326 pairs / 1,320 null pairs, not the docstring's 46 / 1,035). Conjoint rule throughout.

- **Mean FDP = 0.000** — no false promotion in any of the 20 worlds (mean false count 0.00), against a certificate over the 1,320 null pairs of **C^U = 0.00380**.
- **Mean power = 0.750** (90 of 120 plant-instances). Per plant: ρ = 0.50 → 20/20; ρ = 0.35 → 20/20; V = 0.45 → 20/20; η² = 0.25 → 16/20; V = 0.25 → 14/20; η² = 0.10 → 0/20 (below the floor by construction, as the manuscript states — the rule working as defined, not lost power in the reportable region). For reference, probe 2's P1∧P2-only rule (R4, no split-half) reported power 0.808 on its own worlds/seeds; the split-half gate's power cost at n = 500 is therefore ≈0.06, concentrated in the two borderline plants — consistent with probe 3b's per-pair split-half power measurements.

---

## 5. Cluster stress test (file `probe8_cluster.json`) — the certificate FAILS under clustered rows, in two separable ways

Design: 50 all-null worlds at n = 200, rows in 20 clusters of 10. Columns are mutually independent (each column's cluster effects are drawn independently), but within a column, clustered rows are dependent: numeric columns X = a_cluster + e with a ~ N(0, 1), e ~ N(0, 1) (intra-cluster correlation 0.5); categorical columns draw a per-cluster probability vector from Dirichlet(3 × balanced) and then rows from it. The same 50 worlds are evaluated under two splits: a **random row split** (clusters straddle the halves — the deployment default the manuscript's fixed permutation implements) and a **cluster-level split** (whole clusters assigned to halves). The i.i.d. certificate at n = 200 from §3 is C^U = 0.849.

| Split | Mean false V̄ (SE) | FWER | P3-only mean (SE) | Certificate | Violated? |
|---|---|---|---|---|---|
| Random row split | 14.78 (0.62) | 1.00 | 15.44 (0.59) | 0.849 | **Yes — 17.4×** |
| Cluster-level split | 8.94 (0.45) | 1.00 | 20.72 (0.64) | 0.849 | **Yes — 10.5× (conjoint), 24.4× (P3-only)** |

Decomposition (per-half exceedance rates measured inside the clustered worlds, vs the i.i.d. γ̂ at half size 100; "independence prediction" = Σ_cell n_cell γ̂²_{c,cell} using the *clustered* marginals, i.e. what E[#P3] would be if the halves were independent given the clustered marginals):

| Cell | i.i.d. γ̂ | Clustered γ̂, random split | Clustered γ̂, cluster split |
|---|---|---|---|
| num–num | 0.0453 | 0.1689 (3.7×) | 0.2620 (5.8×) |
| cat 3×3 | 0.0180 | 0.0496 (2.8×) | 0.1039 (5.8×) |
| cat 2×2 | 0.0242 | 0.0525 (2.2×) | 0.0829 (3.4×) |
| cat 3×2 | 0.0503 | 0.0919 (1.8×) | 0.1342 (2.7×) |
| num–cat g=3 | 0 (bound 1.28×10⁻⁴) | 0.0006 | 0.0030 |
| num–cat g=2 | 0 (bound 1.28×10⁻⁴) | 0.0001 | 0.0010 |

- *Random split:* independence prediction from the clustered marginals is 8.56, observed P3 mean 15.44 → **dependence factor 1.80**. Both failure channels are open: the marginals are inflated (clustering shrinks the effective sample), *and* the halves are positively dependent in the tail (both halves contain the same clusters, so a chance alignment of cluster means inflates both half statistics together).
- *Cluster-level split:* independence prediction 20.59, observed 20.72 → **dependence factor 1.01**. The factorisation is restored, exactly as the theorem's proof requires (disjoint cluster sets ⇒ independent halves). **But the i.i.d.-calibrated bound is not restored**: the marginals are *more* inflated than under the random split (each half now holds only 10 clusters, so cluster means dominate each half statistic), and the certificate built from i.i.d. γ still fails by an order of magnitude.

Plain conclusion: **with clustered rows, do not use the i.i.d. certificate — under either split.** A cluster-level split fixes the independence the proof needs, but the γ's themselves must then be calibrated under the clustered design (e.g. by cluster-respecting simulation or a cluster bootstrap); with i.i.d. γ the certificate is violated in every one of the 50 worlds (FWER 1.00 under both splits). This does not contradict the theorem — assumption (A1) fails in these worlds by construction — but it is the realistic deployment case for warehouse rows clustered by customer or region, which is why the manuscript's existing i.i.d.-rows caveat next to Proposition 1 should be strengthened to name both channels (marginal inflation and, under row-level splits, cross-half dependence) if the theorem is adopted.

---

## 6. Side finding (not part of the theorem; flagged for the manuscript builder)

`probe2_fdr.py`'s planted design as coded builds **52 columns / 1,326 pairs** (40 base + 12 plant columns), while its own docstring and the manuscript's Table (promotion-rules table) describe the planted worlds as **46 columns / 1,035 pairs**. The 6 plants and all per-rule comparisons are unaffected in kind — §4.2 follows the code — but the printed pair count in the manuscript does not match the shipped probe. This probe did not edit the manuscript; the discrepancy is recorded here for the builder/tester to disposition.

---

## 7. What this does and does not establish

**It establishes** a genuine, if elementary, theorem where the paper previously had a single-pair proposition plus measurements: family-wide control of E[V], FWER, and FDR for *any* rule containing the split-half gate, with no cross-pair dependence assumption and no reliance on BH's dependence conditions — the conjoint rule's error control is thereby decoupled from the BH-under-dependence caveat that the manuscript's own null worlds exhibit (BH alone measured FDR = FWER = 0.15 > q there). The calibration corollary turns that theorem into a computable certificate, and the simulations show the certificate's mechanism is real (P3-only counts track Σγ̂² closely at n = 60 and n = 500) and that the conjoint rule sits far inside it at every tested size, with planted-world FDP 0 and power 0.75 at n = 500. **It does not establish** a tight bound — at small n the certificate is vacuous for FWER/FDR and the conjoint rule's actual protection comes from P1∧P2, which the theorem deliberately does not analyse; it does not cover the text pairs as implemented (full-sample-fitted reduction, Remark R3); it does not transfer to field data with marginal shapes unlike the calibrated cells (assumption C1) or with clustered rows (§5 — where it fails badly and measurably); and it does not rescue the certificate at n = 100 from a ≈1 SE Monte Carlo excursion on the P3-only diagnostic. The honest one-line summary for the paper: *the split-half gate alone certifies family error control at warehouse scale (n ≳ 500 per pair), the conjoint rule is far more conservative than that certificate everywhere tested, and the certificate must not be quoted for clustered rows or text pairs.*

## 8. Reproducibility

- **Script:** `probe8_certificate.py` (this folder). One full run produced all JSONs; sections B and D were then re-run once with identical seeds to add per-world/per-cell diagnostics — headline numbers reproduced exactly.
- **Outputs:** `probe8_validation.json` (facade + engine agreement, 30 + 780 pairs), `probe8_calibration.json` (30 cells × 50,000 pairs), `probe8_worlds_null.json` (50 worlds × 3 sizes), `probe8_worlds_planted.json` (20 worlds), `probe8_cluster.json` (50 worlds × 2 splits), `probe8_run.log` (full-run log).
- **Seeds:** validation 81001/81002; calibration 810000 + 7919·(cell index) per cell (recorded per cell in the JSON); null worlds data 820000 + 100n + w, splits 830000 + 100n + w; planted data 840000 + w, splits 850000 + w; cluster data 860000 + w, row splits 870000 + w, cluster splits 880000 + w.
- **Runtimes** (full run, 190.3 s total): validation 7.5 s; calibration 98.0 s; null worlds 48.6 s; planted worlds 13.8 s; cluster stress 22.4 s. Diagnostics re-run of B + D: 62.4 s.
- **Replicate counts:** calibration 50,000 null pairs per cell (30 cells); worlds 50 (null, per n), 20 (planted), 50 (cluster, both splits on the same worlds).
