# Paper 8 upgrade cycle 2 — the cluster certificate

**Paper:** "A Type-Aware Bottom-Up Discovery Pipeline for Cross-Domain Enterprise KPIs with Spurious-Correlation Filtering" → *Knowledge and Information Systems*.
**Role of this document:** builder output for cycle 2. The cycle-1 family error certificate (Theorem `thm:family`, calibrated in §`sec:certificate` of the manuscript) was scored 7/10 by a fresh Hypothesis Tester against an 8/10 bar, with the named gap: the certificate is never evaluated on a field dataset, and it **fails under clustered rows** — the realistic deployment case — which the manuscript reports as a pure scope failure (cluster-level split: factorisation restored, dependence factor 1.01, but conjoint mean false count 8.94 against an i.i.d. certificate of 0.849, because i.i.d.-calibrated γ's are wrong under clustering). This cycle converts that failure into a proved, validated guarantee: a cluster-level extension of the theorem, a calibration of γ^cl under the clustered design, and a pre-registered validation on fresh clustered worlds. **No manuscript file was edited, authorship was not touched, nothing was submitted anywhere.** All work is in this folder (`hidden_files/paper08_upgrade/cluster/`); the design, seeds, counts and decision rules were fixed in `cluster_preregistration.md` **before** any run.

**Verdict, up front: the cluster certificate works, on fresh seeds, exactly where the i.i.d. certificate failed.** At the cycle-1 clustered design (20 clusters × 10 rows, ICC 0.5), with γ^cl calibrated under clustering at 10 clusters per half (50,000 null pairs per branch cell):

- **P3-only mean false count 21.23 (SE 0.81)** over 60 fresh worlds, against a point prediction Σ(γ̂^cl)² = **20.90** (agreement: +0.39 combined SE) and a simultaneous cluster certificate **C^cl,U = 21.70** — the bound **dominates**, by contrast with the i.i.d. certificate (0.849), which the same worlds violate 25.0× on the P3-only count. The dependence factor is **1.02**, confirming the factorisation under the cluster-level split.
- **Full conjoint rule: mean false count 9.37 (SE 0.51)**, FWER 1.00 — dominated by C^cl,U = 21.70 on the E[V] scale (the FWER/FDR certificate min(1, C^cl,U) = 1 is vacuous at this design; see §4 and the limits in §6), and 11.0× above the i.i.d. certificate on the same worlds. The narrative the paper needs is now measured on fresh seeds: *wrong calibration fails by an order of magnitude; right calibration certifies.*
- **Sensitivity (40 clusters × 5, same ICC):** γ^cl falls as the theorem predicts (num–num 0.264 → 0.157), the certificate tracks (C^cl,U = 7.93), and on 50 fresh worlds the P3-only mean is 7.04 (SE 0.41) and the conjoint mean 2.00 (SE 0.25), both dominated. No failure was found and none was tuned away; the sensitivity dependence factor (0.94) is reported as measured.

---

## 1. The cluster theorem (LaTeX-ready)

Written in the manuscript's notation (Definition 1 promotion rule with P1 Benjamini–Hochberg at q = 0.05, P2 corrected-strength floor τ = 0.20, P3 split-half reproducibility; S_{jh} is pair j's corrected strength on half h). Suggested labels: `ass:cluster`, `thm:cluster`, `cor:cluster-calibration`. The statement **handles variable cluster sizes** — sizes are part of the cluster draw and the halves are counted in clusters — so no fixed-size scoping is needed in the theorem itself (the calibration-matching caveat for variable sizes is real and is stated in the corollary discussion and §6).

```latex
\begin{assumption}[Clustered sampling and cluster-level split]
\label{ass:cluster}
Fix a pair family $\mathcal{F}$ of size $M$, a floor $\tau > 0$, and
the promotion rule of Definition~\ref{def:promotion}. The analysed
sample for pair $j \in \mathcal{F}$ consists of $K_j$ clusters
$Z_1, \dots, Z_{K_j}$, where a cluster $Z_c$ is the full record of
one sampling unit: a (possibly random) number $m_c \ge 1$ of rows,
each row carrying the pair's two column values.
\begin{enumerate}
\item[(B1)] The clusters are independent and identically
distributed draws from a cluster distribution. No assumption is
made about dependence \emph{within} a cluster: a cluster's rows
may be arbitrarily dependent, and the cluster size may be random.
The null $H^{\mathrm{cl}}_{0j}$ states that pair $j$'s two columns
are independent \emph{as clustered processes}: the cluster draw
governing the first column is independent of the cluster draw
governing the second, each column retaining its own marginals and
its own within-cluster dependence.
\item[(B2)] The split assigns \emph{whole clusters} to halves: the
two halves are disjoint sets of clusters, of fixed counts
$K_{j1}$ and $K_{j2}$ with $K_{j1} + K_{j2} = K_j$, and the
assignment is fixed in advance or drawn independently of the data.
No cluster straddles the halves.
\item[(B3)] Each $S_{jh}$ is a function of half $h$'s clusters
only: no quantity fitted on the full sample or on the other half
(no shared preprocessing, clustering, imputation, or
normalisation) enters it.
\end{enumerate}
Write $\gamma^{\mathrm{cl}}_{jh} =
\Pr_{H^{\mathrm{cl}}_{0j}}[S_{jh} \ge \tau]$, pair $j$'s per-half
exceedance probability \emph{under the clustered sampling} at
$K_{jh}$ clusters. Note that $\gamma^{\mathrm{cl}}_{jh}$ is a
functional of the whole cluster law at that cluster count --- the
cluster-size distribution, the column marginals, and the
within-cluster dependence --- not of the row count alone.
\end{assumption}

\begin{theorem}[Family error control under clustered rows]
\label{thm:cluster}
Under Assumption~\ref{ass:cluster}, for every null pair $j$,
\[
\Pr[\text{$j$ is promoted}] \;\le\;
\gamma^{\mathrm{cl}}_{j1}\,\gamma^{\mathrm{cl}}_{j2}.
\]
Consequently, with $\mathcal{H}_0 \subseteq \mathcal{F}$ the set of
true null pairs, $V$ the number of promoted null pairs, $R$ the
total number promoted, and
$C^{\mathrm{cl}} = \sum_{j \in \mathcal{H}_0}
\gamma^{\mathrm{cl}}_{j1}\gamma^{\mathrm{cl}}_{j2}$,
\[
\mathbb{E}[V] \;\le\; C^{\mathrm{cl}}, \qquad
\mathrm{FWER} = \Pr[V \ge 1] \;\le\; \min(1, C^{\mathrm{cl}}),
\qquad
\mathrm{FDR} = \mathbb{E}\!\left[\tfrac{V}{\max(R,1)}\right]
\;\le\; \min(1, C^{\mathrm{cl}}).
\]
No assumption about dependence across pairs is made, and no
property of\/ {\rm (P1)} or {\rm (P2)} is used: the bounds hold for
any promotion rule that has {\rm (P3)} as a necessary condition.
\end{theorem}

\begin{proof}
Fix a null pair $j$. Promotion requires (P3), so
$\{\text{$j$ promoted}\} \subseteq \{S_{j1} \ge \tau,\,
S_{j2} \ge \tau\}$. Condition on the split assignment $A$ (if the
assignment is fixed, this conditioning is trivial). By (B2), $A$
is independent of the data and selects disjoint sets of clusters
of fixed counts $K_{j1}$ and $K_{j2}$. By (B1) the clusters are
i.i.d., and a subcollection selected independently of the cluster
values is itself a collection of i.i.d.\ clusters; the two halves'
subcollections are disjoint, hence independent. Therefore, given
$A$, each $S_{jh}$ --- a function of half $h$'s clusters only, by
(B3) --- has the distribution of the corrected strength computed
on $K_{jh}$ fresh clusters, so
$\Pr[S_{jh} \ge \tau \mid A] = \gamma^{\mathrm{cl}}_{jh}$, and
$S_{j1}$ and $S_{j2}$ are conditionally independent given $A$.
The joint event therefore factorises conditionally, and taking
expectations over $A$,
$\Pr[S_{j1} \ge \tau,\, S_{j2} \ge \tau]
= \gamma^{\mathrm{cl}}_{j1}\gamma^{\mathrm{cl}}_{j2}$.
The family conclusions are now identical to the proof of
Theorem~\ref{thm:family}: linearity of expectation over
$j \in \mathcal{H}_0$ (no independence across pairs) gives
$\mathbb{E}[V] \le C^{\mathrm{cl}}$; the union bound gives
$\mathrm{FWER} \le C^{\mathrm{cl}}$, and a probability is at most
$1$; and $V/\max(R,1) \le \mathbf{1}\{V > 0\}$ pointwise gives
$\mathrm{FDR} \le \mathrm{FWER} \le \min(1, C^{\mathrm{cl}})$.
\end{proof}

\begin{corollary}[Calibrated cluster certificate]
\label{cor:cluster-calibration}
Partition the (branch, clusters-per-half, cluster-design)
combinations that occur in $\mathcal{F}$ into $K$ cluster cells
$c$, where a cell records: the branch (statistic type and level
counts); the number of clusters per half $K_h$; the cluster-size
distribution; the column marginal family; and the within-cluster
dependence law (for example an intra-cluster correlation for
numeric columns, a concentration parameter for categorical
cluster probabilities). Cell $c$ has a true per-half exceedance
$\gamma^{\mathrm{cl}}_c$ under a specified calibration clustered
null distribution. Draw $N_c$ independent null pairs from cell
$c$'s calibration distribution --- each draw a fresh set of
$K_h$ i.i.d.\ clusters carrying two columns that are independent
as clustered processes --- independently of the field sample and
of the other cells, and let $\hat\gamma^{\mathrm{cl},U}_c$ be the
one-sided Clopper--Pearson upper confidence bound for
$\gamma^{\mathrm{cl}}_c$ at level $1 - \delta/K$. Assume in
addition:
\begin{enumerate}
\item[(C1$^{\mathrm{cl}}$)] (cluster-cell dominance) every field
null pair $j$ and half $h$ satisfies
$\gamma^{\mathrm{cl}}_{jh} \le \gamma^{\mathrm{cl}}_{c(j,h)}$,
where $c(j,h)$ is the pair's cell;
\item[(C2$^{\mathrm{cl}}$)] the calibration draws are independent
of the field sample.
\end{enumerate}
Then with probability at least $1 - \delta$ over the calibration
draws, the certificate
$C^{\mathrm{cl},U} = \sum_{j \in \mathcal{F}}
\hat\gamma^{\mathrm{cl},U}_{c(j,1)}
\hat\gamma^{\mathrm{cl},U}_{c(j,2)}$ (the sum may be taken over
all pairs, since $\mathcal{H}_0$ is unknown) satisfies: conditional
on the calibration, $\mathbb{E}[V] \le C^{\mathrm{cl},U}$,
$\mathrm{FWER} \le \min(1, C^{\mathrm{cl},U})$ and
$\mathrm{FDR} \le \min(1, C^{\mathrm{cl},U})$ for the field sample.
\end{corollary}

\begin{proof}
Each Clopper--Pearson bound covers its cell parameter with
probability at least $1 - \delta/K$; by the union bound over the
$K$ cells, the event $E = \{\gamma^{\mathrm{cl}}_c \le
\hat\gamma^{\mathrm{cl},U}_c \text{ for all } c\}$ has probability
at least $1 - \delta$ over the calibration draws. Condition on a
calibration outcome in $E$. The values
$\hat\gamma^{\mathrm{cl},U}_c$ are then fixed numbers, the field
sample is independent of them by (C2$^{\mathrm{cl}}$), and by
(C1$^{\mathrm{cl}}$) and $E$,
$\gamma^{\mathrm{cl}}_{j1}\gamma^{\mathrm{cl}}_{j2} \le
\hat\gamma^{\mathrm{cl},U}_{c(j,1)}
\hat\gamma^{\mathrm{cl},U}_{c(j,2)}$ for every null pair.
Theorem~\ref{thm:cluster} applied to the field sample gives the
three bounds with $C^{\mathrm{cl}} \le C^{\mathrm{cl},U}$.
\end{proof}
```

**What (C1^cl) matching means, exactly** (prose for the manuscript, following the corollary). γ^cl is not determined by branch and row count: at the *same* 100 rows per half, 10 clusters of 10 at ICC 0.5 and 20 clusters of 5 at the same ICC give different γ^cl (0.264 vs 0.157 for numeric–numeric in §5 below). Cluster-cell dominance is therefore a design-matching requirement, not a sample-size requirement. A calibration cell can dominate a field pair only if it matches the field design on branch, clusters per half K_h, and cluster-size distribution, matches its column marginals, **and** its within-cluster dependence induces a per-half exceedance at least as large as the field's. In the parametric families used here (additive Gaussian cluster effects; Dirichlet cluster probabilities) the operational route is to calibrate at a dependence level at or above the field's — ICC at least the field ICC, Dirichlet concentration at most the field's — with sizes and marginals matched; monotonicity of γ^cl in the dependence parameter is a property of the family in use and should be verified numerically for that family, not assumed as a general theorem. Where the field design is only partially known, dominance can be claimed only by calibrating at a worst-case envelope over a stated range of designs, and the certificate inherits that envelope's conservatism. Where no cell dominates, the corollary does not apply and no certificate should be quoted.

**Remarks.** (R1) With all cluster sizes equal to 1, Assumption `ass:cluster` reduces to Assumption `ass:theorem` and Theorem `thm:cluster` to Theorem `thm:family`: the cluster theorem strictly generalises the i.i.d. one. Cluster size never enters the proof — the halves are counted in clusters, which is why variable sizes are covered. (R2) A **random row split is not covered**: halves then share clusters, (B2) fails, and §4 measures the resulting tail dependence factor at 1.90. (R3) The text pairs remain outside scope as implemented, for the (B3) reason of cycle 1 (reduction fitted on the full sample before the split). (R4) When C^cl > 1 — the strongly clustered, few-cluster regime — the FWER/FDR parts are vacuous and only E[V] ≤ C^cl carries information; FWER-level certification requires C^cl,U < 1, i.e. materially more clusters per half (§5 shows the direction: doubling clusters per half at fixed row count cut C^cl,U from 21.70 to 7.93; it did not reach below 1 at this ICC and family size).


---

## 2. Simulation setup and statistic validation

All simulations are in `probe10_cluster_certificate.py`, run with the archived venv interpreter (`hidden_files/paper08_work/phik_venv/bin/python`; Python 3.12.3, numpy 2.5.4, pandas 3.1.0, scipy 1.18.1) — the cycle-1 environment. Worlds are generated by cycle 1's `clustered_columns`, **imported unmodified**, and evaluated by cycle 1's `eval_pairs`/`run_world`, also imported unmodified, so the world design and the statistic definitions are exactly cycle 1's: 40 columns (24 numeric, 8 categorical-3, 8 categorical-2; 780 pairs: 276 num–num, 28 cat 3×3, 28 cat 2×2, 64 cat 3×2, 192 num–cat g = 3, 192 num–cat g = 2); columns mutually independent; numeric X = a_c + e with a ~ N(0, 1), e ~ N(0, 1) (ICC 0.5); categorical per-cluster probability vectors ~ Dirichlet(3 × balanced), rows drawn from their cluster's vector. The conjoint rule is Definition 1 exactly (BH q = 0.05 ∧ full-sample strength ≥ 0.20 ∧ split-half ≥ 0.20 on each half); P3-only counts are recorded alongside, as in cycle 1.

**Validation of the statistics (file `probe10_validation.json`).** (a) Cycle 1's scalar reference statistics were checked against `facade_v2.correlate` on **30 clustered null pairs** (10 num–num, 10 cat–cat across 3×3 / 2×2 / 3×2, 10 num–cat across g = 3 / 2; half-sized primary-design samples of 10 clusters × 10 rows; seed 91001): maximum absolute difference **1.04×10⁻¹⁶** in strength and **1.91×10⁻¹⁴** in p-value — the same order as cycle 1's i.i.d. validation (2.78×10⁻¹⁷ / 4.89×10⁻¹⁴). (b) Probe 10's batched calibration kernels were checked against those scalar references on **2,000 clustered draws per branch** (seeds 91002 + 131·branch index): maximum absolute strength difference **6.08×10⁻¹⁶**, and exceedance counts agreed **exactly in all six branches**. (Import note as cycle 1: the venv's numpy 2.5 removed the `np.math` alias `facade_v2` uses; the harness restores it as the stdlib `math` module before importing; `facade_v2.py` itself is unmodified.)

**Pre-registration.** `cluster_preregistration.md` was written before the run and fixes: the two designs (primary P: 20 × 10; sensitivity S: 40 × 5, chosen over the ICC-0.25 alternative because it isolates the clusters-per-half dependence at identical ICC, Dirichlet concentration, and rows per half); 50,000 calibration pairs per cell with its justification; 60 primary worlds and 50 sensitivity worlds, all with fresh seeds disjoint from cycle 1 (860000–880000 range) and from calibration (910000–970000 range, allocated per section); τ = 0.20, δ = 0.05, Bonferroni over each design's 6 cells; and the checks to report. There were **no deviations**.

## 3. Calibration of γ^cl (file `probe10_calibration_cluster.json`)

One calibration draw = one half-sized clustered sample (primary: 10 clusters × 10 rows; sensitivity: 20 clusters × 5 rows) with two mutually independent columns of the branch type, each with its own independent cluster effects. γ̂^cl is the exceedance rate of the corrected strength at τ = 0.20; γ̂^cl,U_sim is the simultaneous Clopper–Pearson upper bound (per-cell level 1 − 0.05/6; each design is a separate 0.95-simultaneous certificate). N = 50,000 per cell, as cycle 1: worst-case Monte Carlo SE ≈ 0.0022, and the induced delta-method SE on the point certificate is 0.29 (primary) / 0.14 (sensitivity) — small against the world-mean SEs that dominate §4's comparisons.

| Design (half) | Branch | k / 50,000 | γ̂^cl | CP95 upper | γ̂^cl,U_sim | i.i.d. γ̂ (m = 100, cycle 1) |
|---|---|---|---|---|---|---|
| P: 10 cl × 10 | num–num | 13,178 | 0.26356 | 0.26682 | 0.26831 | 0.04526 |
| P: 10 cl × 10 | cat 3×3 | 4,906 | 0.09812 | 0.10034 | 0.10135 | 0.01800 |
| P: 10 cl × 10 | cat 2×2 | 3,575 | 0.07150 | 0.07342 | 0.07430 | 0.02424 |
| P: 10 cl × 10 | cat 3×2 | 7,158 | 0.14316 | 0.14576 | 0.14695 | 0.05028 |
| P: 10 cl × 10 | num–cat g=3 | 127 | 0.00254 | 0.00294 | 0.00313 | 0 (bound 1.28×10⁻⁴) |
| P: 10 cl × 10 | num–cat g=2 | 35 | 0.00070 | 0.00093 | 0.00104 | 0 (bound 1.28×10⁻⁴) |
| S: 20 cl × 5 | num–num | 7,872 | 0.15744 | 0.16014 | 0.16138 | — |
| S: 20 cl × 5 | cat 3×3 | 2,358 | 0.04716 | 0.04875 | 0.04948 | — |
| S: 20 cl × 5 | cat 2×2 | 2,360 | 0.04720 | 0.04879 | 0.04952 | — |
| S: 20 cl × 5 | cat 3×2 | 4,709 | 0.09418 | 0.09636 | 0.09735 | — |
| S: 20 cl × 5 | num–cat g=3 | 16 | 0.00032 | 0.00049 | 0.00057 | — |
| S: 20 cl × 5 | num–cat g=2 | 6 | 0.00012 | 0.00024 | 0.00030 | — |

Certificates (780-pair family): **primary C^cl,pt = Σ n_cell (γ̂^cl)² = 20.90** (delta SE 0.29), **C^cl,U = Σ n_cell (γ̂^cl,U_sim)² = 21.70**; **sensitivity C^cl,pt = 7.53** (delta SE 0.14), **C^cl,U = 7.93**. For contrast, the cycle-1 i.i.d. certificate at total n = 200, recomputed from `probe8_calibration.json` with the same pair counts, is C^pt_iid = 0.753 and **C^U_iid = 0.849** (reproduces cycle 1's printed 0.849 exactly). The clustered γ's are 3–6× their i.i.d. counterparts in the categorical/numeric cells at the same row count — that inflation, invisible to any row-count-only calibration, is the entire mechanism of the cycle-1 failure.

## 4. Validation on fresh clustered worlds (file `probe10_worlds_cluster.json`)

60 fresh all-null worlds, design P (data seeds 930000 + w, cluster-split seeds 940000 + w, row-split seeds 950000 + w, w = 0…59). In all-null worlds FDR = FWER.

| Quantity (cluster-level split) | Observed (SE) | Point prediction Σ(γ̂^cl)² | Cluster certificate C^cl,U | i.i.d. certificate C^U (n = 200) |
|---|---|---|---|---|
| P3-only mean false count | **21.23 (0.81)** | 20.90 | **21.70 — dominates ✓** | 0.849 — violated 25.0× |
| Conjoint mean false count V̄ | **9.37 (0.51)** | — (≤ P3 expectation) | **21.70 — dominates ✓** | 0.849 — violated 11.0× |
| FWER (= FDR, all-null) | 1.00 (0.00) | — | min(1, C^cl,U) = 1 — holds, vacuously | min(1, 0.849) = 0.849 — violated |

Reading, plainly:

- **The factorisation holds and the bound dominates — both halves of the cycle-1 failure are closed.** The P3-only mean (21.23) sits +0.39 combined standard errors from the point prediction (20.90; world SE 0.81 and calibration delta SE 0.29 combined in quadrature) and 0.57 world-SEs *below* the simultaneous certificate (21.70). The dependence factor — observed P3 mean ÷ independence prediction from the worlds' own measured marginals (20.75) — is **1.02**, replicating cycle 1's cluster-split factor of 1.01 on fresh seeds. The conjoint rule's mean (9.37) is dominated with the expected slack: as in cycle 1, the certificate credits only P3, and BH ∧ floor remove roughly 56% of the gate's false passes here.
- **The i.i.d. certificate fails on the same fresh worlds, as pre-registered:** conjoint 9.37 vs 0.849 (11.0×; cycle 1 measured 10.5× on its seeds) and P3-only 21.23 vs 0.849 (25.0×; cycle 1: 24.4×). Same worlds, same rule, same split — the only difference is which sampling model the γ's were calibrated under. Wrong calibration fails; right calibration certifies.
- **FWER honesty:** every one of the 60 worlds had at least one conjoint false promotion (FWER 1.00), because at 10 clusters per half and ICC 0.5 the gate's per-pair exceedances are large (γ^cl up to 0.26) and the family is wide. The cluster certificate's FWER/FDR statement, min(1, 21.70) = 1, is therefore **vacuous at this design** — it holds trivially and certifies nothing about FWER. What is certified, non-vacuously, is E[V] ≤ 21.70 against an observed conjoint mean of 9.37. A deployment needing FWER < 1 needs C^cl,U < 1, which needs more clusters per half (or a narrower family), not a better proof. The manuscript should carry the cluster result as an E[V] certificate at this design, exactly as cycle 1's small-n i.i.d. certificates were E[V]-only in substance.
- **Within-world marginal check** (pre-registered check 5): per-cell per-half exceedance rates measured inside the validation worlds vs the independent calibration γ̂^cl — num–num 0.2619 vs 0.2636; cat 3×3 0.0923 vs 0.0981; cat 2×2 0.0711 vs 0.0715; cat 3×2 0.1497 vs 0.1432; num–cat g=3 0.00330 vs 0.00254; num–cat g=2 0.00117 vs 0.00070. All agree within ≈2 naive standard errors (naive SEs understate the truth here because pairs sharing a column are dependent, so the world-side rates are overdispersed); the two num–cat cells sit at the ≈2 SE end, with absolute differences ≤ 0.0008 whose combined contribution to the certificate sum is ≈ 0.001 against a certificate of 21.70 — immaterial to every comparison above, and reported rather than smoothed.
- **Random row split (descriptive only, NOT certified):** on the same 60 worlds, conjoint mean 15.30 (SE 0.60), P3-only mean 16.22 (SE 0.57), and the dependence factor against the row-split marginals' independence prediction (8.56) is **1.90** — the cross-half tail dependence cycle 1 measured at 1.80, reproduced. Both row-split means happen to fall below C^cl,U = 21.70 numerically; that is slack in a bound computed for a different split, **not** certification — γ^cl is defined for cluster-split halves, the factorisation the theorem needs measurably fails (factor 1.90), and against the certificate that *would* be quoted for row splits from marginals alone the P3 count fails by 1.90×. The manuscript's deployment instruction stands: cluster-level split or no certificate.

## 5. Sensitivity: 40 clusters × 5 rows (files `probe10_calibration_cluster.json`, `probe10_worlds_sensitivity.json`)

Design S holds total rows (200), rows per half (100), ICC (0.5), Dirichlet concentration (3), and family fixed, and changes only the cluster geometry: 40 clusters × 5 rows, halves of 20 clusters. If γ^cl were a row-count quantity, nothing would change; the theorem says it is a cluster-count quantity.

- **γ^cl falls in every cell** (table in §3): num–num 0.26356 → 0.15744; cat 3×2 0.14316 → 0.09418; cat 3×3 0.09812 → 0.04716; cat 2×2 0.07150 → 0.04720; num–cat g=3 0.00254 → 0.00032; g=2 0.00070 → 0.00012. Certificate: C^cl,pt = 7.53, **C^cl,U = 7.93** (vs 21.70 at 20 × 10).
- **50 fresh worlds** (data seeds 960000 + w, cluster-split seeds 970000 + w): P3-only mean **7.04 (SE 0.41)** — −1.15 combined SEs from the point prediction 7.53 and 2.20 world-SEs below the certificate 7.93: **dominates ✓**. Conjoint mean **2.00 (SE 0.25)**, dominated ✓. FWER **0.78 (SE 0.06)** — the certificate's FWER bound is still formally min(1, 7.93) = 1 (vacuous), but the conjoint rule's own FWER has dropped below 1 because P1∧P2's slack grows as the gate tightens; reported as an observation about the conjoint rule, not as a certified FWER claim. World-marginal rates again match calibration (num–num 0.1565 vs 0.1574).
- **Dependence factor 0.94** (independence prediction from the worlds' own marginals: 7.53, essentially identical to the calibration point prediction 7.53; observed 7.04). A factor below 1 is sampling noise in the ratio (the observed mean is 1.2 world-SEs below its own marginal prediction), not evidence of negative dependence; it is reported as measured, per the pre-registration's no-tuning rule.
- The i.i.d. certificate (0.849) also fails on these worlds (P3-only 8.3×, conjoint 2.4×): fewer, larger clusters are not the only geometry where row-count calibration misleads.

**Conclusion of the sensitivity run:** the certificate *tracks* — recalibrating γ^cl at the new cluster count moves the certificate from 21.70 to 7.93 and the observed false counts move with it, staying inside the bound in both geometries. The certificate is a function of the cluster design, and must be recomputed when the design changes; there is no geometry-free cluster certificate.


## 6. Honest limits — what the cluster certificate does NOT cover

The following paragraph is written for direct adaptation into the manuscript (e.g. as the closing paragraph of the certificate section or in the limitations), and is also the builder's own statement of scope:

> The cluster certificate covers exactly the design it is calibrated for, and its limits are the assumption's limits. **First, the cluster structure must be known in advance:** every analysed row must carry a cluster identifier, the split must be made on those identifiers, and the clusters must be independent draws — if rows are grouped wrongly, if the true clustering is latent, or if clusters are themselves dependent (spatially or through a network), (B1)/(B2) fail and no certificate is claimed. If the two columns of a pair are clustered by *different* units (say by customer for one and by region for the other), the assumption as stated does not apply; a split on the common refinement of the two clusterings would be needed, and we do not treat that case. **Second, calibration matching is an assumption, (C1^cl), not a fact the field data verify:** the validation above is exact because the simulated worlds use the calibrated cluster law; a field deployment whose intra-cluster dependence is stronger, whose cluster sizes are more variable, or whose marginals are heavier-tailed than the calibration cells can have field γ^cl above the calibrated bounds, and the certificate can then fail in the field exactly as the i.i.d. certificate failed here. The theorem does handle variable cluster sizes (the proof counts clusters, not rows), but a fixed-size calibration cell does not thereby dominate a variable-size field design — the calibration must reproduce the field's size distribution, or dominate it by an argued envelope. **Third, the certificate remains a (P3)-only bound:** it is conservative for the conjoint rule by the full slack of (P1)∧(P2) (here, an observed conjoint mean of 9.37 against a bound of 21.70), and at strongly clustered, few-cluster designs it certifies an expected false count only — the FWER/FDR bounds are min(1, C^cl,U) and are vacuous whenever C^cl,U ≥ 1, as in both geometries validated here. **Fourth**, the statement is superpopulation (over fresh cluster draws), not a statement about re-splits of one fixed table; the pipeline's text pairs remain outside it as implemented, by (B3); and pairwise complete-case selection that depends on a pair's own values can break the i.i.d.-clusters assumption for the analysed subsample, as in the i.i.d. case. Finally, γ^cl must be **recalibrated whenever the cluster design changes** — cluster count per half, size distribution, or dependence level: §5 shows the certificate moving from 21.70 to 7.93 under a geometry change at fixed row count, so a certificate carried across designs is not a certificate.

## 7. What this does and does not establish

**It establishes** the missing half of cycle 1's cluster finding. Cycle 1 showed the i.i.d. certificate failing under clustered rows and identified the two channels (marginal inflation; cross-half dependence under row splits) with the cluster-level split repairing the second channel only. This cycle proves that the first channel is also repairable, by calibration rather than by assumption: with γ defined and calibrated under the clustered sampling at the half's cluster count, the *same* product bound is a theorem (Thm `thm:cluster`, which strictly generalises the i.i.d. theorem and needs no fixed cluster sizes), the Clopper–Pearson machinery carries over unchanged (Cor `cor:cluster-calibration`), and the resulting certificate dominates on fresh worlds in two cluster geometries, with the P3-only diagnostic tracking the point prediction within ≈1 combined SE in both. The cycle-1 "pure scope failure" paragraph in the manuscript can now be replaced by: a theorem, a calibration table, a validation table, and a deployment rule (cluster-level split + matched calibration, or no certificate). **It does not establish** field transfer: no field dataset with known cluster ground truth was evaluated (the Hypothesis Tester's other named gap — a field evaluation — is only partly answered: the theorem now *covers* clustered field data under (C1^cl), but dominance for any specific field table remains an assumption to be argued per deployment, as §6 states). It does not give FWER-level control at realistic clustered geometries (C^cl,U ≫ 1 in both validated designs; the honest claim is E[V]). It does not cover row-level splits, latent or cross-classified clustering, dependent clusters, or the text pairs. And the validated dependence levels are the simulated ones (Gaussian additive effects at ICC 0.5; Dirichlet(3) categorical probabilities): heavier-tailed or more strongly clustered field designs need their own calibration cells before any number is quoted.

## 8. Reproducibility

- **Script:** `probe10_cluster_certificate.py` (this folder). One run produced all outputs; no re-runs, no seed changes, no post-hoc design changes.
- **Pre-registration:** `cluster_preregistration.md`, written before the run (designs, seeds, counts, checks).
- **Outputs:** `probe10_validation.json` (facade + kernel agreement), `probe10_calibration_cluster.json` (2 designs × 6 cells × 50,000 pairs), `probe10_worlds_cluster.json` (60 primary worlds, both splits, per-world counts), `probe10_worlds_sensitivity.json` (50 design-S worlds), `probe10_run.log` (full-run log).
- **Seeds:** validation 91001 (facade pairs), 91002 + 131·(branch index) (kernel checks); calibration 910000 + 7919·(branch index) (design P) and 920000 + 7919·(branch index) (design S), branch order num–num, cat 3×3, cat 2×2, cat 3×2, num–cat g=3, num–cat g=2, recorded per cell in the JSON; primary worlds data 930000 + w, cluster splits 940000 + w, row splits 950000 + w (w = 0…59); sensitivity worlds data 960000 + w, cluster splits 970000 + w (w = 0…49). Cycle-1 seeds (calibration 810000-range; cluster worlds 860000/870000/880000) were not reused.
- **Runtimes** (single run, 16.8 s total): validation 1.3 s; calibration (both designs) 6.8 s; primary worlds 6.1 s; sensitivity worlds 2.4 s. (Much faster than cycle 1's 190 s full run because the calibration kernels are batched and only the two needed designs are calibrated; the world engine is cycle 1's, unchanged.)
- **Replicate counts:** calibration 50,000 null pairs per cell (12 cells); worlds 60 (primary, each under both splits) and 50 (sensitivity, cluster split).
- **Environment:** `hidden_files/paper08_work/phik_venv/bin/python` — Python 3.12.3, numpy 2.5.4, scipy 1.18.1, pandas 3.1.0 (identical to cycle 1).
