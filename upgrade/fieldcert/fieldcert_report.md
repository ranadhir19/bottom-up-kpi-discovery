# Paper 8, Upgrade Cycle 2 — the family error certificate, evaluated on the field families

**Task.** The fresh Hypothesis Tester scored the upgraded paper 7/10 (bar 8/10) with one named residual: the family error certificate (Theorem `thm:family`; Corollary `cor:calibration`) is *"validated only in synthetic worlds and never actually evaluated on a field dataset."* This report evaluates the certificate on the paper's three field families — Olist (probe 5), Telco (probe 7), Bank (probe 9) — exactly as the Corollary defines it, and adds the small-*n* contrast that keeps the field evaluation honest.

**Deliverables in this folder:** `probe11_field_certificates.py` (the computation; archived venv interpreter, as for probe 8), `field_certificates.json` (every count, cell value and result below), this report. No manuscript file was edited, nothing was submitted, authorship was not touched.

**Headline results.** With $f = 1.2793\times10^{-4}$ the $m=500$ simultaneous Clopper–Pearson bound (identical in all six calibrated branches, §2), every field pair contributes $f^2 = 1.6366\times10^{-8}$, so $C^U = M f^2$:

| Family | $M$ | $C^U$ | Implied bounds $\mathbb{E}[V]$, FWER, FDR |
|---|---|---|---|
| Olist | 107 | $1.75\times10^{-6}$ | each $\le 1.75\times10^{-6}$ |
| Olist, (A3)-covered pairs only | 105 | $1.72\times10^{-6}$ | each $\le 1.72\times10^{-6}$ |
| Telco | 231 | $3.78\times10^{-6}$ | each $\le 3.78\times10^{-6}$ |
| Bank | 253 | $4.14\times10^{-6}$ | each $\le 4.14\times10^{-6}$ |

Observed on all three: **0 known false promotions** — the only known nulls on field data are the negative-control pairs (27 Olist, 41 Telco, 43 Bank), and none was promoted (promoted sets: 13, 88, 33). The true $\mathcal{H}_0$ is unknown on field data, so this is a **bound evaluation, not a validation against ground truth** — stated in the LaTeX text itself, §5.

**The informative contrast.** The same Olist family composition at small per-pair $n$ gives $C^U = 5.19$ at $n=60$ and $1.73$ at $n=100$ (FWER/FDR bounds vacuous at both), $4.09\times10^{-4}$ at $n=500$, and $1.75\times10^{-6}$ at $n=1000$ — the last *calibrated*, where the field value of the same size is *extrapolated* (§2). The certificate therefore has no bite at the small per-pair sizes where bottom-up discovery is usually attempted — there (P1)∧(P2) must do the work, exactly as Sections `sec:fdr` and `sec:certificate` found in synthetic worlds — and is effectively zero at warehouse scale, where the split-half gate alone already certifies the family. That, not the smallness of the field numbers by itself, is what the field evaluation establishes.

---

## 1. Method

**Definition applied.** Corollary `cor:calibration`: $C^U = \sum_{j \in \mathcal{F}} \hat\gamma^{U}_{c(j,1)}\,\hat\gamma^{U}_{c(j,2)}$, summed over all pairs ($\mathcal{H}_0$ unknown), with $\hat\gamma^{U}$ the simultaneous one-sided Clopper–Pearson bounds of probe 8 (50,000 null pairs per cell; Bonferroni over the 30 cells; $\delta = 0.05$; archived in `hidden_files/paper08_upgrade/theorem/probe8_calibration.json`, printed as Table `tab:gamma`).

**Per-pair inputs, from the artifacts only.** Branch and analysed $n_j$ are read per pair from the archived full pair tables (`probe5_full_pair_table.csv`, family = the 105 column pairs + the 2 TF–IDF text pairs, selected exactly as probe 5 selects them; `probe7_full_pair_table.csv`; `probe9_full_pair_table.csv`). Branch is the table's `method` (`pearson_r` → numeric–numeric; `eta_squared` → numeric–categorical; `cramers_v` → categorical–categorical). Half sizes follow the probes' split convention: $n_{j1} = \lfloor n_j/2 \rfloor$, $n_{j2} = n_j - n_{j1}$. Categorical cardinalities (for the fine cell distributions and the dimension-matched sensitivity) were recomputed from the raw field data following the probes' own table constructions (probe 5's Olist rebuild replicated line-for-line for construction only; Telco and Bank read directly).

**Checks run by the script (all asserted, all passed).**
1. *Half-size check (tasked):* **no field pair has a half size ≤ 500.** Smallest half anywhere: 3,516 (Telco pairs involving `TotalCharges`, $n_j = 7{,}032$). So no pair needed its own small-$m$ cell; every pair uses the §2 extrapolation.
2. *Monotonicity check:* $\hat\gamma^{U}$ is non-increasing in $m$ over $m = 30, 50, 100, 250, 500$ in **all six** calibrated branches (verified cell by cell from the JSON; the numeric–categorical branches are flat at the floor from $m=100$ on).
3. *Floor check:* at $m = 500$ all six branches recorded 0 exceedances in 50,000 pairs, so all six $m=500$ bounds equal the single CP floor value (check 2 of the script).
4. Family sizes reproduce the manuscripts': $M = 107, 231, 253$; branch counts reproduce the combinatorial counts from the probes' column lists (e.g. Bank: $\binom{12}{2}=66$ numeric–numeric, $12\times11=132$ numeric–categorical, $\binom{11}{2}=55$ categorical–categorical).

## 2. The extrapolation, stated plainly

Every field half size (3,516–20,594) exceeds the largest calibrated half size ($m = 500$). The per-half bound used for all field pairs is the $m=500$ simultaneous bound,

$$f \;=\; 1 - (0.05/30)^{1/50000} \;\approx\; 1.2793\times10^{-4},$$

which is the Clopper–Pearson upper bound at **0 exceedances in 50,000 pairs** — a confidence bound, not an estimate — and is the same number in every branch because every branch's $m=500$ cell had count 0.

*Justification.* At fixed $\tau = 0.20$, the exceedance probability of these statistics decreases with sample size, and the calibrated bounds do so in every one of the six branches across the full calibrated range $m = 30 \dots 500$ (check 2 above). Taking the $m=500$ bound as the per-half bound for larger $m$ is therefore conservative *if* that decay continues.

*What is **not** claimed.* (i) Monotonicity beyond $m = 500$ is an **empirical regularity plus the CP-floor argument, not a proved fact** — no theorem in the paper bounds $\gamma$ at uncalibrated sizes. (ii) Assumption **(C1) cell dominance remains an assumption for field marginals**: the calibration cells are balanced-marginal worlds with at most $3\times3$ tables and $g \le 3$ groups, while the field families contain tables up to $27 \times 7$ (Olist) and $12 \times 10$ (Bank) and numeric–categorical pairs with $g$ up to 27 (§3) — shapes for which no calibrated cell exists, whose true $\gamma$ at field half sizes could in principle exceed $f$. (iii) (A1) independent rows is likewise assumed, not verified, on field data; the cluster stress test of Section `sec:certificate` shows what its failure costs. The LaTeX text in §5 carries all three caveats in its body.

## 3. Arithmetic path (every number traceable to the artifacts)

Per-pair contribution at field scale: $f^2 = (1.2793040931\times10^{-4})^2 = 1.6366190\times10^{-8}$. Hence $C^U = M \times 1.6366190\times10^{-8}$ — the branch mix is recorded in full because the Corollary is defined per pair, and it drops out of the field total only because all $m=500$ cells coincide at the floor.

**Olist** ($M = 107$; pair table: 105 pairs at $n_j = 25{,}000$, halves 12,500; 2 text pairs at $n_j = 9{,}378$, halves 4,689).
- Branches: numeric–numeric 55 $=\binom{11}{2}$; numeric–categorical 44 $= 11 \times 4$ (by $g$: $g{=}2$: 11, $g{=}4$: 11, $g{=}7$: 11, $g{=}27$: 11); categorical–categorical 8 (column pairs $2{\times}4, 2{\times}7, 2{\times}27, 4{\times}7, 4{\times}27, 7{\times}27$, one each; text pairs $3{\times}3$ and $3{\times}4$, one each).
- $C^U = 107 \times 1.6366190\times10^{-8} = 1.7512\times10^{-6}$. Excluding the 2 text pairs — outside the Theorem's scope as implemented, since the reduction is fitted on the full commented subset before the split, violating (A3) — $C^U = 105 \times f^2 = 1.7184\times10^{-6}$ for the covered family. (Including the text pairs in the sum, as the Corollary's all-pairs form does, keeps the number a valid bound for the covered pairs; it does not certify the text pairs.)
- Observed: 13 promoted; 27 negative-control pairs, 0 promoted.

**Telco** ($M = 231$; 210 pairs at $n_j = 7{,}043$, halves 3,521/3,522; 21 pairs involving `TotalCharges` at $n_j = 7{,}032$, halves 3,516).
- Branches: numeric–numeric 10 $=\binom{5}{2}$; numeric–categorical 85 $= 5 \times 17$ (by $g$: $g{=}2$: 35, $g{=}3$: 45, $g{=}4$: 5); categorical–categorical 136 $=\binom{17}{2}$ ($2{\times}2$: 21, $2{\times}3$: 63, $2{\times}4$: 7, $3{\times}3$: 36, $3{\times}4$: 9).
- $C^U = 231 \times 1.6366190\times10^{-8} = 3.7806\times10^{-6}$.
- Observed: 88 promoted; 41 negative-control pairs, 0 promoted.

**Bank** ($M = 253$; all pairs at $n_j = 41{,}188$, halves 20,594).
- Branches: numeric–numeric 66 $=\binom{12}{2}$; numeric–categorical 132 $= 12 \times 11$ (by $g$: $g{=}2$: 24, $g{=}3$: 48, $g{=}4$: 12, $g{=}5$: 12, $g{=}8$: 12, $g{=}10$: 12, $g{=}12$: 12); categorical–categorical 55 $=\binom{11}{2}$ (dimensions from $2{\times}2$ to $12{\times}10$; full distribution in the JSON).
- $C^U = 253 \times 1.6366190\times10^{-8} = 4.1406\times10^{-6}$.
- Observed: 33 promoted; 43 negative-control pairs, 0 promoted.

**Small-$n$ contrast (Olist composition: 55 / 44 / 8; per pair, total $n$, half $m = n/2$; calibrated cells, no extrapolation).** Headline uses the *branch envelope* — within each pipeline branch the largest simultaneous bound among its calibrated cells at $m$ (numeric–categorical: $\max(g{=}2, g{=}3)$; categorical–categorical: $\max(3{\times}3, 2{\times}2, 3{\times}2)$), the conservative choice inside the calibrated set:

| total $n$ | $m$ | $\hat\gamma^U$ NN | $\hat\gamma^U$ NC (env.) | $\hat\gamma^U$ CC (env.) | $C^U$ (envelope) | $C^U$ (dimension-matched) |
|---|---|---|---|---|---|---|
| 60 | 30 | 0.29466 | 0.019689 | 0.22278 | **5.1896** $= 55(0.29466)^2 + 44(0.019689)^2 + 8(0.22278)^2$ | 5.1715 |
| 100 | 50 | 0.16905 | 0.0023840 | 0.14183 | **1.7330** $= 55(0.16905)^2 + 44(0.0023840)^2 + 8(0.14183)^2$ | 1.7228 |
| 500 | 250 | 0.0024302 | 0.00012793 | 0.0032287 | **$4.0894\times10^{-4}$** | $3.9862\times10^{-4}$ |
| 1000 | 500 | 0.00012793 | 0.00012793 | 0.00012793 | **$1.7512\times10^{-6}$** $= 107 f^2$ | $1.7512\times10^{-6}$ |

The dimension-matched sensitivity (numeric–categorical: $g{=}2$ cell for $g=2$, $g{=}3$ cell for $g \ge 3$; categorical–categorical: exact cell where the dimensions are among $(2,2), (2,3), (3,3)$, envelope otherwise) differs from the envelope by under 3% at every size, so the contrast does not depend on the assignment rule. At $n = 60$ and $n = 100$ the implied FWER/FDR bounds are $\min(1, C^U) = 1$: **vacuous** — only $\mathbb{E}[V] \le C^U$ carries information there, and it allows several false promotions in expectation from the gate alone.

## 4. What this establishes, and what it does not

**Establishes.** (i) The Corollary's certificate is computable on the paper's actual field families from archived artifacts, and at field scale it is effectively zero ($10^{-6}$ order): conditional on the Corollary's assumptions, the split-half gate alone bounds the expected number of false promotions in each field family below $5\times10^{-6}$, with FWER and FDR bounded by the same numbers. (ii) The observed field outcomes (no control promoted on any dataset) are consistent with the bounds, at the only points where field ground truth exists. (iii) The contrast ladder locates the certificate's bite: it is vacuous for FWER/FDR at total $n \le 100$ for a family of this size and mix, first informative at $n = 500$, and overwhelming at warehouse scale — quantifying, on the paper's own family composition, the division of labour the paper argues for: at small per-pair $n$ the conjoint rule's protection comes from (P1)∧(P2); at warehouse scale (P3) alone already certifies the family.

**Does not establish.** (i) Any validation against field ground truth: $\mathcal{H}_0$ is unknown on field data, most promoted pairs are genuinely associated (that is the paper's own classification finding), and consistency at the control pairs is the only empirical contact. (ii) Coverage of the text pairs as implemented ((A3) violated; the 105-pair Olist subtotal is the covered statement). (iii) Anything under clustered or otherwise non-i.i.d. rows ((A1)); the stress test's failure stands as the boundary. (iv) A proved bound at field half sizes: the field numbers rest on the §2 extrapolation and on (C1) for field marginal shapes — both flagged, neither discharged. A revision that wants the field certificate proved rather than extrapolated must calibrate $\gamma$ at half sizes in the thousands and under field-like (unbalanced, wide-table) marginals; that is a larger simulation programme than this evaluation, and the LaTeX text says so by presenting the field numbers under their assumptions rather than as a validation.

## 5. LaTeX-ready text (for the manuscript builder; not inserted)

Placement suggestion: a new subsection at the end of Section `sec:field` (after `sec:cpss`), so the field runs it evaluates are already defined; it refers back to Section `sec:certificate` and Table `tab:gamma`. Label: `sec:fieldcert`; table label: `tab:fieldcert`. All caveats are in the body text, as tasked — there are no footnotes in this block.

```latex
\subsection{The certificate on the field families}
\label{sec:fieldcert}

The certificate of Corollary~\ref{cor:calibration} has so far
been exercised only in the synthetic worlds of
Section~\ref{sec:certificate}. It can also be computed, pair by
pair, on the three field families of Section~\ref{sec:field},
and doing so is instructive in both directions. Every field
pair's analysed size is far above the calibration range: the
smallest half size in any field family is $3{,}516$ (Telco
pairs involving total charges, $n_j = 7{,}032$), against a
largest calibrated half size of $m = 500$, so no field pair
uses an interpolated cell and the computation turns on one
explicit extrapolation. In every calibrated branch the
simultaneous bound $\hat\gamma^{U}$ of Table~\ref{tab:gamma}
is non-increasing in the half size across $m = 30$ to $500$,
and at $m = 500$ every branch recorded $0$ exceedances in
$50{,}000$ null pairs, so all six bounds coincide at the
Clopper--Pearson floor $f = 1.28 \times 10^{-4}$. We take $f$
as the per-half bound at all larger half sizes. Two limits of
that step belong here rather than in a footnote. The decay of
$\gamma$ beyond $m = 500$ is an empirical regularity of the
calibrated cells together with the floor argument, not a
proved fact; and cell dominance (C1) remains an assumption for
field marginals, whose tables run to $27 \times 7$ and whose
numeric--categorical pairs have up to $27$ groups, shapes no
calibrated cell contains. The numbers below are therefore the
certificate \emph{under its stated assumptions, with the
$m = 500$ bound extrapolated}, and the field runs' rows are
observational records, not verified independent draws, so
(A1) is assumed there exactly as (C1) is.

With that scope, the computation is short: each pair
contributes $f^2 = 1.64 \times 10^{-8}$, so
$C^{U} = M f^2$ for a family of size $M$, whatever its branch
mix (Table~\ref{tab:fieldcert}). The certificates are
$1.75 \times 10^{-6}$ on Olist ($M = 107$; $55$
numeric--numeric, $44$ numeric--categorical and $8$
categorical--categorical pairs), $3.78 \times 10^{-6}$ on
Telco ($M = 231$) and $4.14 \times 10^{-6}$ on Bank
($M = 253$): conditional on the calibration, the expected
number of false promotions attributable to the gate in each
field family is below $5 \times 10^{-6}$, and the
family-wise and false discovery rates are bounded by the same
numbers. The Olist figure includes its two text pairs in the
sum, as the Corollary's all-pairs form permits; since the
text reduction as implemented is fitted before the split, the
Theorem does not cover those pairs, and the covered subtotal
over the remaining $105$ pairs is $1.72 \times 10^{-6}$. The
observed outcomes are consistent with the bounds at the only
points where field ground truth exists: the negative-control
pairs ($27$ on Olist, $41$ on Telco, $43$ on Bank) are known
nulls, and none was promoted on any dataset. Beyond the
controls the true null set is unknown on field data, so this
is an evaluation of the bound, not a validation against
ground truth, and we report it as such.

The contrast half of Table~\ref{tab:fieldcert} keeps the
field numbers honest. Recomputing the certificate for the
Olist family composition at small total sizes from the
calibrated cells (branch-envelope cells; a dimension-matched
assignment changes no value by more than $3\%$) gives
$C^{U} = 5.19$ at $n = 60$ and $1.73$ at $n = 100$, where the
family-wise and false discovery bounds are vacuous and only
$\mathbb{E}[V] \le C^{U}$ carries information, $4.09 \times
10^{-4}$ at $n = 500$, and $1.75 \times 10^{-6}$ at
$n = 1000$ --- the same value as the field certificate, but
there calibrated rather than extrapolated. The certificate's
bite is therefore at small $n$ in reverse: at the per-pair
sizes where bottom-up discovery is usually attempted, the
gate alone certifies almost nothing and the conjoint rule's
protection is the work of (P1) and (P2), as
Sections~\ref{sec:fdr} and~\ref{sec:certificate} found in
synthetic worlds; at warehouse scale the gate alone already
certifies the family, and the conjoint rule's clean field
record is overdetermined. That quantification, on the
paper's own family composition, is what the field evaluation
adds: it shows where each component of the rule is the one
doing the work.

\begin{table}[t]
\centering
\caption{The calibrated certificate $C^{U}$ of
Corollary~\ref{cor:calibration} on the field families
(upper panel) and, for contrast, on the Olist family
composition at small per-pair sizes from the calibrated
cells of Table~\ref{tab:gamma} (lower panel). The bounds are
$\mathbb{E}[V] \le C^{U}$ and
$\mathrm{FWER}, \mathrm{FDR} \le \min(1, C^{U})$. Field
half sizes all exceed the largest calibrated half size
($500$); the field rows use the $m = 500$ simultaneous bound
$f = 1.28 \times 10^{-4}$ per half, extrapolated as described
in the text, under which each pair contributes
$f^2 = 1.64 \times 10^{-8}$. At $n = 60$ and $n = 100$ the
family-wise and false discovery bounds are vacuous
($\min(1, C^{U}) = 1$).}
\label{tab:fieldcert}
\begin{tabular}{lrrr}
\toprule
Family & $M$ & Half sizes & $C^{U}$ \\
\midrule
\multicolumn{4}{l}{\emph{Field families (per-pair $n$ as run)}} \\
Olist ($n_j = 25{,}000$; text $9{,}378$) & $107$ & $4{,}689$--$12{,}500$ & $1.75 \times 10^{-6}$ \\
Olist, covered pairs only (no text) & $105$ & $4{,}689$--$12{,}500$ & $1.72 \times 10^{-6}$ \\
Telco ($n_j = 7{,}043$ / $7{,}032$) & $231$ & $3{,}516$--$3{,}522$ & $3.78 \times 10^{-6}$ \\
Bank ($n_j = 41{,}188$) & $253$ & $20{,}594$ & $4.14 \times 10^{-6}$ \\
\midrule
\multicolumn{4}{l}{\emph{Olist composition at small $n$ (calibrated cells)}} \\
Total $n = 60$ & $107$ & $30$ & $5.19$ \\
Total $n = 100$ & $107$ & $50$ & $1.73$ \\
Total $n = 500$ & $107$ & $250$ & $4.09 \times 10^{-4}$ \\
Total $n = 1000$ & $107$ & $500$ & $1.75 \times 10^{-6}$ \\
\bottomrule
\end{tabular}
\end{table}
```

## 6. Reproducibility

- **Script:** `probe11_field_certificates.py` (this folder). One run produced `field_certificates.json`; all numbers in this report are read from that JSON.
- **Interpreter:** `hidden_files/paper08_work/phik_venv/bin/python` (Python 3.12, pandas 3.1.0), the probe-8 environment. Runtime seconds.
- **Inputs:** probe 8 calibration JSON; probe 5/7/9 full pair tables and run JSONs; raw field CSVs for cardinalities only (Olist construction replicated from `probe5_olist_run.py`; no statistic was re-estimated and no field run was repeated).
- **Cross-checks inside the script (asserted):** 30 calibration cells; monotonicity of $\hat\gamma^{U}$ in $m$ in all six branches; all $m=500$ cells at 0 exceedances with a single shared bound; family sizes 107/231/253; every field half $> 500$; branch counts equal to the combinatorial counts of the probes' column lists.
- **Files:** `probe11_field_certificates.py`, `field_certificates.json`, `fieldcert_report.md` (this file), all in `hidden_files/paper08_upgrade/fieldcert/`.
