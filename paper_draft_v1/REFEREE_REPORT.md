# Referee report on "Convergence of the Glowinski decomposition–coordination method for the g-Laplacian: assumptions on the solution, the discrete problem, the Lavrentiev phenomenon, and numerical experiments" (new/main.tex, 26 pp.)

Scope of the check:
- I read every proof of main.tex line by line.
- I compared every attribution with materials/REFERENCES_A/B (and their status flags), BRIEF.md, REVIEW.md, NOTES.md and TODO_VERIFY.md.
- I recomputed the numerics with Python from materials/convergence_order.csv and new/results/*.
- I re-ran the Lyapunov script for t²/2+t⁴/4 and reproduced the reported numbers exactly.

Line numbers refer to main.tex.

---

## (a) Verdict

**Major revision (mathematically sound; scope, wording and attribution must be fixed).**

The mathematical core is correct. I found no error in any of the following proofs:
- Proposition 2.3 (saddle points);
- Lemmas 2.4–2.7;
- the Lyapunov identity, its signs, and the cluster-point and Opial steps of Theorem A;
- Theorem B(a), and Theorem B(b) including its strengthenings (L^s, s<2, and the strong multiplier convergence);
- the Hahn–Banach multiplier lemma and Prop. 2.15(i) (for p<2, (S) ⇔ minimum attained on H¹₀ ⇔ u_p ∈ H¹₀);
- Prop. 2.13 (Glowinski's (iii) and (iv) fail for p<2);
- Theorem C;
- Lemma 5.1.

No integral-functional theorem of Rockafellar is needed: the pointwise ∂F lemma is proved from scratch, correctly.

All tables agree with the raw data to the printed digits.

The must-fix problems are of four kinds:
1. Scope errors, where theorems are applied outside their stated hypotheses: the nonhomogeneous manufactured test, and point evaluation of a.e.-defined data in Theorem C.
2. A few false side statements: the t^p log^a superlinearity condition, the "isotropic ⇒ uniformly elliptic" sentence, and the description of (S) for double phase.
3. Attributions stronger than the sources: Table 1 "yes" for Cianchi–Maz'ya 2011 in d = 2, Breit–Diening–Schwarzacher for discontinuous coefficients, and the "0.8" order read off Del Pezzo–Martínez.
4. Internal inconsistencies in the numerics bookkeeping: the count of non-converged runs, and unmarked fits that use non-converged runs.

In addition, the paper never proves, and should say explicitly that it does not prove, that the P1 minimisers u_h converge to u, even for the radial test problems. Lemma 5.1 gives only the lower bound J(v_h) ≥ J(u).

On novelty, a referee at a good numerical-analysis journal will remark on three points. Theorem A is a careful but standard Fejér/Opial argument; the authors themselves cite Moradifam–Nachman and Bartels–Milicevic. Theorem C is standard finite-dimensional ADMM. The genuinely new pieces are Theorem B(b), Prop. 2.15 and the careful (S)/Lavrentiev discussion. The paper should say this plainly in the introduction. As written, "Contributions (1) and (4)" read as more novel than they are.

---

## (b) ERRORS (must fix)

**E1. The manufactured test is outside the continuous theory as stated, and (A2) is violated by the shifted functional (l. 580, l. 725).**
- §6.4 asserts "(S) holds by Lemma 3.2 and Theorem B(a) applies" for u = e^{x₁+x₂} with nonhomogeneous Dirichlet data.
- Sections 2–3 are formulated only for V = H¹₀, and Lemma 3.2 assumes u ∈ H¹₀.
- Remark 4.1 handles nonhomogeneous data by the shift F̃_h(q) = F_h(q+∇ũ_h). But F̃_h(0) = F_h(∇ũ_h) ≠ 0, so (A2) ("F ≥ 0 = F(0)") fails as stated.
- (A2) is used for properness and coercivity in Lemma 2.6(ii) and for the admissible initialisation 0 ∈ ∂F(0). Nothing else needs F(0) = 0.

Fix:
1. Weaken (A2).
2. Add a continuous analogue of Remark 4.1.
3. Note that the code's initialisation ν₀ = η₁ = 0 corresponds, in the shifted variables, to q₀ = −∇ũ_h, η₁ = 0 ∈ ∂F̃_h(q₀) = ∂F_h(0). This is admissible, and by Remark 2.9 any initialisation is fine anyway.

```latex
\item\label{A2} $F\colon H\to[0,+\infty]$ is convex, lower semicontinuous and proper ($\dom F\neq\emptyset$);
```
In Lemma 2.6(ii), replace "proper ($\Psi_n(0)<\infty$)" with "proper (since $\dom F\ne\emptyset$)". In §2.3, replace the sentence on the admissible initialisation with:
"If moreover $F(0)=0$, then $0\in\partial F(0)$ and $\nu_0=\eta_1=0$ is admissible; by Remark~\ref{rem:init} the condition $\eta_1\in\partial F(\nu_0)$ is not needed anyway."

Add after Remark 4.1:
```latex
\begin{remark}[Nonhomogeneous data at the continuous level]
If $u=u_D$ on $\partial\Omega$ with $u_D\in H^1(\Omega)$ and $F(\nabla u_D)<\infty$, write $u=u_D+v$, $v\in H^1_0$, and
$\tilde F(q)=F(q+\nabla u_D)$, $\tilde\ell(v)=(f,v)$. Then $(V,H,B,\tilde F,\tilde\ell)$ satisfies \ref{A1}--\ref{A3},
and Proposition~\ref{prop:saddle} and Theorems~\ref{thm:A}--\ref{thm:B} apply in the variables $(v,q,\eta)$; Lemma~\ref{lem:explicit}
holds with ``$u\in H^1_0$'' replaced by ``$u-u_D\in H^1_0$''.
\end{remark}
```
(Theorem B's proofs use only ∂F̃(q) = ∂F(q+∇u_D) = {T(·, q+∇u_D)}, so they go through unchanged.)

**E2. Theorem C evaluates a.e.-defined data at points (Def. 2.1, l. 543–553).**
- A Young integrand is defined only for a.e. x. Then G(x_T,·) and f(x_T) are meaningless: the barycentres form a null set, and f ∈ L² has no point values.

Fix: add to §4.1, before F_h:
```latex
For the discretisation we assume that $G(x,\cdot)$ is a convex Young function for \emph{every} $x\in\overline\Omega$
(as in all our examples) and that $f$ has well-defined values at the barycentres (e.g.\ $f\in C(\overline\Omega)$);
otherwise $F_h$ and $\ell_h$ are to be defined with any representative, or $\ell_h$ by exact integration.
```

**E3. False superlinearity condition for t^p log^a(1+t) (l. 565).**
- "when it is a convex Young function with p+a>1" is wrong. For example p = 1/2, a = 1 gives t^{1/2} log(1+t) = o(t) at ∞, which is not superlinear.
- Superlinearity at ∞ needs p > 1, or p = 1 and a > 0.

Replace with:
```latex
to $t^p\log^a(1+t)$ whenever it is a convex Young function with $p>1$, or $p=1$ and $a>0$,
```

**E4. The ranges stated for the Lieberman condition are sampled ranges, not the true ones (l. 478).**
- I recomputed with mpmath on [10⁻⁸, 10⁸]. For t² log(1+t), tg'/g → 2 as t → 0 and → 1⁺ as t → ∞, so the range is (1,2], not [1.05,2]: 1.05 is just the value at t = 10⁸.
- For t³/log(1+t), tg'/g → 1 as t → 0 and → 2⁻ as t → ∞, so the range is [1,2), not [1,1.95].
- The condition still holds, with δ = 1, g₀ = 2, so the conclusion is unchanged.

Replace with:
```latex
for $t^2\log(1+t)$ ($tg'/g\in(1,2]$) and $t^3/\log(1+t)$ ($tg'/g\in[1,2)$), as one checks from the asymptotics at $0$ and $\infty$ and the monotonicity of $tg'/g$ (verified numerically on $[10^{-8},10^8]$).
```
Also note that h = g/t is nondecreasing for both (tg'/g ≥ 1). I checked this numerically; it is consistent with the "B(a)" entries in Table 3.

**E5. False sentence on the isotropic case (l. 508).**
- "These are not of the isotropic form G(|ξ|), which is uniformly elliptic in the sense of (3.1)…" is wrong. Isotropic G need not satisfy (3.1): e^t−1−t and e^t−1 are isotropic, and the paper itself says they violate (3.1) (l. 499).

Replace with:
```latex
These integrands are anisotropic; they do not concern isotropic energies $G(\abs\xi)$ satisfying \eqref{eq:lieberman}, for which
Fact~\ref{fact:lieberman} gives bounded gradients. (For isotropic $G$ violating \eqref{eq:lieberman}, such as $e^t-1-t$, we know of no counterexample and no positive result; see the paragraph \emph{Energies not covered}.)
```

**E6. Incomplete description of (S) for double phase (l. 497).**
- "it requires u ∈ H¹₀ and, if q>2, a|∇u|^{q−1} ∈ L²" omits the p-phase. The flux is |∇u|^{p−2}∇u + a|∇u|^{q−2}∇u, and a may vanish.

Replace with:
```latex
it requires $u\in H^1_0$ and $\abs{\nabla u}^{p-1}+a\abs{\nabla u}^{q-1}\in L^2(\Omega)$ up to the boundary
(automatic for the first term if $p\le2$).
```

**E7. Table 1, row 2: "yes" is not supported in the paper's own dimension (l. 520, l. 483).**
- For Cianchi–Maz'ya 2011 the paper has only the abstract. The range n ≥ 3 in REFERENCES_A is MEM.
- The JEMS 2014 Thm 3.5 that the row also cites is restricted to n ≥ 3 and 2 ≤ p < n.
- The paper works in d = 2 (§4, §6), so "yes" is not supported for the case that matters.

Change the entry to:
```latex
Autonomous $G$, $\Omega$ convex, $f\in L^{n,1}$ & yes, conditionally (range of $n$ unverified; $n\ge3$ for \cite{CianchiMazya2014JEMS}) & \cite{CianchiMazya2011} (abstract only), \cite[Thm.~3.5]{CianchiMazya2014JEMS} ($p$-type, $2\le p<n$, $n\ge3$).\\
```
and add "(possibly only for $n\ge3$)" to Fact 3.4(a).

**E8. The abstract and Contribution (2) state Theorem B(a) without its hypotheses (l. 57, l. 84).**
- "strong convergence in H¹₀ when g(x,t)/t is nondecreasing" omits three hypotheses: G(x,·) ∈ C¹ with g(x,0) = 0; the uniform positivity h(x,t) ≥ h₀(t) > 0; and (S).
- Positivity is not implied by monotonicity. For G(t) = (t−1)₊², h ≡ 0 on (0,1]; F is then not strictly convex, the minimiser need not be unique, and the theorem's uniqueness conclusion fails.

Replace in the abstract with:
```latex
(B) strong convergence in $H^1_0$ when $G(x,\cdot)\in C^1$, $g(x,0)=0$ and $g(x,t)/t$ is nondecreasing and bounded below by a positive function of $t$ uniformly in $x$,
```
Make the same change in Contribution (2).

**E9. Attributions in §6.9 that the sources do not support (l. 798).**
- (i) "The reduced orders for the discontinuous coefficient are consistent with the regularity-limited theory (e.g. [BreitDieningSchwarzacher2015])". BDS 2015 treats Hölder exponents p ∈ C^{0,α} and quasi-norm errors, under F(·,∇v) ∈ W^{1,2}. It says nothing about jump coefficients, L² errors of u, or nodal maxima. Delete the sentence, or replace it with "The reduced orders for the discontinuous coefficient are not covered by any estimate we know."
- (ii) "for discontinuous Galerkin discretisations, L² errors from which an order of about 0.8 can be read off [DelPezzoMartinez2014DCM]". REFERENCES_B D1 says explicitly that this is the librarian's estimate ("my estimate … not stated by the authors"). It is derived from an error reduction of about 4.1 between m = 100 and m = 2916, under an unverified meaning of m. Either recompute it from their table, stating the assumption on m (if m is the number of elements, the order is log 4.1 / log √29.16 ≈ 0.84), or delete it. Suggested replacement: "report $L^2$ errors on a sequence of meshes without fitting an order \cite{DelPezzoMartinez2014DCM}".

**E10. The bookkeeping of non-converged runs is inconsistent (l. 682, l. 660, l. 801, Table 7).**
- Limitation (iv) says "Ten runs did not reach the tolerance". Table 7 contains two more non-converged runs: manufactured r = 1 (5000*) and dead core r = 1 (6000*). The correct count is twelve.
- Table 3 reports the expdead fit (1.97 / 1.01 / 1.11) without marking that two of its three meshes did not converge. The pto1 fits do the same: all eight runs are non-converged.

Fix:
1. In (iv): "Twelve runs did not reach the tolerance: the two finest meshes for $e^t-1$ in Table~\ref{tab:iterations}, all eight runs for $p\in\{1.1,1.2\}$, and the runs with $r=1$ in Table~\ref{tab:r}; …".
2. Mark the expdead row of Table 3 with $^\dagger$ and add to the caption: "$^\dagger$: two of the three meshes used in the fit did not reach the tolerance (Table~\ref{tab:iterations})."

**E11. Notation clash in a displayed formula (l. 618) and in Fig. 1.**
- "The exact radial solutions $u(r)=\int_r^1 g^{-1}(cs/2)\,ds$" uses r, which everywhere else is the penalty parameter. The Fig. 1 legend "a = 1_{r>1/2}" does the same.

Write instead:
```latex
$u(x)=\int_{\abs x}^1 g^{-1}(cs/2)\,ds$
```
and use "$\abs x>1/2$" in the figure legend (scripts/make_figures.py).

---

## (c) GAPS and overstatements (should fix)

**G1. Item 2 of the brief: what the radial remark and Lemma 5.1 actually prove.**

*What is correct.*
- Remark 5.2 is correct for the stated data: unit disk, d = 2, f = c > 0. I checked each step:
  - η = −cx/2 gives −div η = c, using d = 2; the remark is tied to the disk, so this is fine.
  - Where τ > 0, η = g₀(|x|,τ)·(−x/|x|) = −cx/2 is the gradient of Φ at ∇u = −τ x/|x|.
  - In the dead core of e^t−1, |η| = c|x|/2 ≤ 1 = g(0), so η ∈ ∂Φ(0) (Lemma 2.5(i)).
  - The two-phase interface is null.
- So (S) holds by Lemma 3.2, including for e^t−1−t, e^t−1 and the indicator coefficient.
- Lemma 5.1 is correct. Integrability of Φ(x,∇u) is automatic: Φ(x,∇u) ≤ η·∇u ∈ L^∞, by the subgradient inequality at 0. Say this in the proof, since the lemma does not assume sup_x G(x,M) < ∞.
- What Lemma 5.1 proves: the W^{1,1}-minimum is attained by a Lipschitz function, so it equals the infimum over Lipschitz functions (no gap between W^{1,1} and W^{1,∞}), and J(v_h) ≥ J(u) for every conforming v_h.

*Gap (a): the smooth-function version of "no gap".*
- §5 defines m_𝓗 as the infimum over "smooth (or Lipschitz)" functions. For the Lipschitz class the statement is immediate. For C^∞ (or C_c^∞ perturbations) one also needs u_k smooth with J(u_k) → J(u).
- This holds here: dilate to compact support in B₁, then mollify, so that ‖∇u_k‖_∞ ≤ C‖∇u‖_∞ and ∇u_k → ∇u a.e. Dominated convergence then applies with majorant sup_x G(x,C‖∇u‖_∞) < ∞, which is true for all examples, including the indicator two-phase.
- Add one sentence to that effect, or define the gap with the Lipschitz class only.

*Gap (b): convergence of u_h → u is not proved.*
- Nothing in the paper proves that the P1 minimisers u_h converge to u. Lemma 5.1 gives only the lower bound, and Remark 4.2 explicitly disclaims the approximation property.
- Yet §5 ("free of this obstruction"), Contribution (5) and the "theory" column of Table 3 invite the reader to think it has been shown.
- Either state explicitly in §5 and in "What is not proved": "We do not prove that $u_h\to u$ as $h\to0$; Lemma~\ref{lem:nogap} gives only $J(v_h)\ge J(u)$ for conforming $v_h$."
- Or add the following short result. It is true and easy for the autonomous cases, because there J_h = J exactly on V_h extended by 0: G is sampled exactly, and the barycentre rule is exact for constant f on P1.

```latex
\begin{proposition}\label{prop:uh-conv}
Let $G=G(t)$ be autonomous, $\Omega=B_1\subset\R^2$, $f=c$, and let $u,\eta$ be as in Remark~\ref{rem:radial-S}, with $\tau$ continuous on $[0,1]$.
Let $(\mathcal T_h)$ be shape-regular triangulations of inscribed polygons $\Omega_h$ with boundary nodes on $\partial B_1$, and let $u_h$
minimise $J_h$ (extended by $0$ to $B_1$). Then $J(u_h)\to J(u)$ and
$\int_{B_1}\big[\Phi(\nabla u_h)-\Phi(\nabla u)-\eta\cdot\nabla(u_h-u)\big]\dx\to0$; consequently $\nabla u_h\to\nabla u$ in measure, and in $L^2$ if $g(t)/t$ is nondecreasing and positive.
\end{proposition}
\begin{proof}
$J_h=J$ on $V_h$. By Lemma~\ref{lem:nogap}, $J(u_h)\ge J(u)$; by minimality $J(u_h)\le J(I_hu)$, and $J(I_hu)\to J(u)$ by dominated convergence, since $\norm{\nabla I_hu}_\infty\le C\norm{\nabla u}_\infty$ (shape regularity), $\nabla I_hu\to\nabla u$ a.e.\ ($u\in C^1(\overline{B_1}\setminus\{0\})$) and $\abs{B_1\setminus\Omega_h}\to0$. The integrand of the Bregman remainder is nonnegative, and its integral equals $J(u_h)-J(u)$ (Lemma~\ref{lem:nogap}). The last assertions follow from strict convexity, respectively from $\Phi(b)-\Phi(a)-T(a)\cdot(b-a)=\int_0^1\psi(a+s(b-a),a)\,\frac{ds}s$ and Lemma~\ref{lem:psi}(i), as in the proof of Theorem~\ref{thm:B}(a).
\end{proof}
```
- For the indicator two-phase problem one needs in addition the consistency of barycentre sampling on the O(h)-measure set of triangles cut by |x| = 1/2, which requires a uniform L^∞ bound on ∇u_h. That bound is not available, so say that this case remains unproved.

**G2. Why p ≥ 2 was "required": for p > 2, Glowinski's (iv) is not obviously satisfied either (l. 79, §2.6).**
- For p > 2, F = +∞ on L² \ L^p, so F cannot be Gâteaux differentiable on H.
- Also, F − (c/2)‖·‖² is not convex near 0, so no obvious splitting F = F₀ + F₁ with F₀ differentiable and uniformly convex on H exists.
- Hence the reading of (iv) taken from GPT16 does not literally hold for p > 2 either. Either the book's hypothesis is weaker (e.g. differentiability on dom F), or the application in [ABR, Thm 4.3] needs comment.
- The paper says "for G = t^p/p this holds only for p ≥ 2". This is literally true of uniform convexity and superlinearity, but it suggests that the theorem applies for p ≥ 2.
- Add: "For $p>2$, $F$ is not finite on $L^2$, so (iv) must be understood with $F_0$ differentiable on $\dom F$ (or similar); this has to be checked in \cite{Glowinski1984ChVI}." Also add this to TODO_VERIFY item 1.

**G3. GPT16 "strongly" is inferred from context.** REFERENCES_B C4 records the conclusion "(u^n,p^n) → (u,Bu) in V×H (strongly, by context)". In §2.6 write "converge (strongly, as we read \cite{GlowinskiPanTai2016})", or flag it.

**G4. Boyd et al., App. A (MEM) is cited as fact (Rem. 2.10).** Either verify it, or write "and, as far as we recall, in the finite-dimensional analysis of \cite[App.~A]{Boyd2011ADMM}". It is listed in TODO_VERIFY, but the text is unhedged.

**G5. Remark 5.2, "Since τ is bounded on [0,1]".** This is not implied by the stated hypotheses when G₀ depends on s (τ(s) = g₀(s,·)⁻¹(cs/2) could blow up). Add the hypothesis "and $\inf_s g_0(s,t)\to\infty$ as $t\to\infty$ (true in all our examples, where $g_0(s,t)\ge g_{\min}(t)$ with $g_{\min}(t)=t^{p-1}$, $t$, or $e^t-1$)". Measurability of τ should also be noted (it is a monotone inverse).

**G6. Remark 2.11, "no argument based only on the general ADMM/Douglas–Rachford theory can give strong convergence".**
- This is defensible, because Bùi–Combettes' operators are subdifferentials and DR on (∂f, ∂g) is ADMM with B = I on the conjugate problem. But the paper does not explain that mapping.
- Either add that sentence, or soften to "the general ADMM/Douglas–Rachford theory gives only weak convergence".

**G7. The p → 1 runs (§6.3).**
- "The errors are nevertheless meaningful" is not supported. The paper's own r-study (§6.6) shows that a small primal residual does not control the iteration error.
- Add evidence (the change in the L² error between iterations 3000 and 6000, or a run to 2·10⁴ iterations on one mesh), or weaken to "the errors are reported with the caveat …".
- Delete "consistent with Theorem B(b), whose estimate degenerates as p → 1". Theorem B(b) contains no rate, so it cannot be "consistent" with an iteration count.

**G8. The stopping criterion.** The text identifies that ‖ν_n − ν_{n−1}‖ (the dual residual) is not monitored. Monitoring both residuals is the standard ADMM practice (Boyd et al. §3.3), and it would remove the r = 100 anomaly. Either recommend it or rerun with it.

**G9. Lyapunov check (§6.7).**
- The identity is purely algebraic once ν* = ∇u* and the multiplier equation hold. Its verification to 10⁻¹⁶ is a check of the implementation, not of the theorem. Say so.
- For t²/2+t⁴/4 and e^t−1−t, E_n/E₀ reaches 10⁻²⁴ after about 44 iterations (from lyapunov_curves.npz). The remaining ~950 of the 1000 checked iterations are pure round-off.
- Report the range of n over which E_n/E₀ > 10⁻²⁰, or restrict the table to it.
- "E_n never increases" should read "never increases by more than 3×10⁻²⁸ E₀", since Table 8 shows positive maxima.
- The "linear decay" is measured against a proxy taken from iterate L. Say that the rate near n ≈ L is artificial; the figure is cut at 10⁻²⁴, which mitigates this.

**G10. Observation (5) and Remark 4.2 on iteration counts.**
- "grows with h⁻¹" reads as proportional to h⁻¹. The data grow sublinearly: for p = 1.5, log₂ of the successive ratios is 0.94, 0.64, 0.41, 0.52.
- The counts also grow for t² log(1+t) (36 → 60) and for e^t−1.
- Write "grows as $h$ decreases (roughly like $h^{-1/2}$ on the finest meshes)" and add t² log(1+t) to the list.

**G11. Observation (2)** lists the reduced nodal-maximum orders for t⁴/4 (1.17) and t² log(1+t) (1.37), but omits e^t−1 (1.11), which is also highlighted in Fig. 1. Add it.

**G12. §6.9 heuristic.** "g(x,t)/t bounded above and below on the range of |∇u|, which is the case away from the origin": this is false in the whole dead core |x| < 1/4 of e^t−1, where g(t)/t → ∞ as t → 0, and for t^{1.5} near the origin. Add "(except in the dead core of $e^t-1$)".

**G13. "What is not proved" should list:**
- (iv) convergence u_h → u as h → 0, also for the test problems (see G1);
- (v) the d = 2 restriction of Prop. 2.15 (or p ≥ 2d/(d+2)), and of §4.

**G14. The Table 3 "theory" column** refers to the continuous iteration, while all runs are of the discrete iteration (to which Theorem C applies in every case). Add to the caption: "the discrete iteration is covered by Theorem~\ref{thm:C} in all cases".

**G15. Theorem B(a), proof.**
- The replacement of h₀ by ess inf_x h(x,t) is unnecessary and slightly delicate (an uncountable family of null sets).
- The proof only uses h(x,|e_n|/2) ≥ h(x,δ/2) ≥ h₀(δ/2) for a.e. x, at the single value t = δ/2, together with monotonicity in t. Delete the first sentence of the proof of (a), or rephrase it accordingly.

**G16. Fact 3.6 (De Filippis–Mingione Thm 1.3).** REFERENCES_A records the additional condition "q/p < p if n = 2". Include it, since d = 2 is the case of interest.

**G17. Bibliography versus the REFERENCES files.**
- No key is missing or unused.
- All entries match the REFERENCES BibTeX in authors, titles and journals.
- The following data are flagged MEM there but printed as fact:
  - GPT16 pp. 19–94;
  - Svaiter vol. and pages;
  - Bùi–Combettes pp. 1118–1120;
  - BOS2022 vol. 151, pp. 779–805;
  - BGS2022 pp. 209–240;
  - BDS2020 Art. 180;
  - BCM2018 Art. 62;
  - DFM2021 vol. 242;
  - Hästö–Ok vol. 24(4);
  - Cianchi–Maz'ya JEMS 16;
  - BCDM2024 Art. 35;
  - Gabay pp. 299–331;
  - Acerbi–Mingione 2002;
  - Zhikov 1986/1995.
- Antonini 2026 is the only entry with a DOI, and that DOI is unverified. Either give DOIs uniformly (they are in the REFERENCES files) or none.
- All of these are in TODO_VERIFY, but the authors must clear them before submission.

---

## (d) Minor points and typos

1. **Notation clashes.**
   - H is both L²(Ω;ℝ^d) and H(x,t) = t^p + a t^q (Fact 3.6); 𝓗 is Zhikov's class. Use ℋ_{dp} or A(x,t) for the double phase.
   - n is the dimension in §3 ("Throughout, n denotes the dimension"), but the iteration index everywhere else, and d is the dimension in §§1, 2, 5. Use d throughout.
   - h is both the mesh size and h(x,t) = g/t. The clash is acute in §6.4: "($h=1+(1+x_1)t^2\ge1$)" with h = 1/n. Rename the latter, e.g. k(x,t).
   - T is both a triangle and the map T(x,ξ) (both appear in §4).
   - F is both the energy and F(a) = (φ'(|a|)/|a|)^{1/2} a in §6.9; use V(a), as in Diening–Kreuzer.
   - p is used for elements of H_h in "(p,q)_h" (l. 543).
   - ψ is used in Lemma 2.7 and in Theorem B.
   - t is both the Young variable and the order in "error ≈ Ch^t" (l. 618).
2. l. 620: "on these meshes it is superclose to $u_h$". It is I_h u that is superclose to u_h. Also, the quantity is the exact L² norm of the P1 function I_h u − u_h (fem.norm_L2 uses the exact midpoint rule), not a "discrete L² norm". Rephrase: "The $L^2$ norm of $I_hu-u_h$ is not the error: $I_hu$ is superclose to $u_h$".
3. l. 620: "spurious orders around 2.5". The nodal order for t⁴/4 is 2.08; write "between 2.1 and 2.5".
4. Theorem A, last line: "$\bar u=u$ is the solution": u is undefined there. Write "$\bar u$ is the unique minimiser".
5. Prop. 2.15(ii) and the following paragraph: under (ii) the converse also holds (for p > 2, (S) ⇔ ∇u_p ∈ L^{2p−2}). Stating it makes "(S) is a higher integrability requirement" precise.
6. Lemma 2.14: closedness of R is not used. Harmless.
7. Table 7 caption: for manufactured r = 300 the error is 4.4135·10⁻⁵, so "4.413·10⁻⁵ in all converged runs" should read "4.413–4.414".
8. §6.2: the target h is used instead of the actual mesh size; this is stated in the limitations, but give the actual h_max on the finest mesh.
9. §6.1: the 9·10⁻⁹ accuracy check of the exact radial solutions is not backed by any file in results/. Add the script output.
10. l. 508: Marcellini 1989/1991 contain mainly positive results. Cite Marcellini 1991 §6 specifically for the counterexample, as REFERENCES_A B3 does.
11. hyperref "Token not allowed in a PDF string" warnings come from math in section titles. Use `\texorpdfstring`.
12. TODO placeholders (addresses, funding) remain.

---

## (e) Items checked and confirmed correct

**Proofs.**
1. Lemma 2.2 (lsc via subsequence, a.e. convergence and Fatou; weak lsc via Mazur).
2. Prop. 2.3 (both directions; the multiplier step uses (u*+sw, ν*), as REVIEW E1 required; independence of r).
3. Lemma 2.4 (pointwise ∂F), proved directly with rational z, the sets Ω^z_M and arbitrary A ⊂ Ω^z_M. This uses boundedness of Ω (so that z 1_A ∈ L²) and measurability of G(·,t). No Rockafellar citation is needed. The consequence (2.3), "∂F(q) = {T} if T ∈ H and ∅ otherwise", follows correctly.
4. Lemma 2.5 (i)–(iii).
5. Lemma 2.6 (Riesz; existence and uniqueness of the prox; its optimality condition proved by a directional argument; (2.6)).
6. Lemma 2.7 (i)–(iv), including the nonsmooth threshold and the case ŵ·z < 0 in (iv).
7. Corollary 2.8.
8. Theorem A, Step 1: the four-term expansion and both polarisation identities recomputed. Step 2: m_n ≥ 0 and the cross term ≥ 0 by monotonicity, using η_n ∈ ∂F(ν_{n−1}) for n ≥ 1. Step 3: the limit of (η_{n+1},ν_n), then the subgradient inequality plus weak lsc. Step 4: |x_n − z|² = 2E_n(z), and the Opial argument.
9. Remark 2.9.
10. Lemma 2.12(i): the identity was expanded and verified; the lower bound holds.
11. Lemma 2.12(ii): the eigenvalues of DT; integrability at the zero crossing; the constant 2^p.
12. Theorem B(a): the splitting estimate, ν* = ν̄, and the uniqueness argument via the subgradient equality.
13. Theorem B(b): Hölder with exponents 2/p and 2/(2−p); interpolation for L^s, s < 2; ‖η_{n+1}−η̄‖² ≤ 4^p|Ω|^{1−2(p−1)/p}‖e_n‖_p^{2(p−1)}, using 2(p−1) ≤ p. These strengthenings are correct.
14. Prop. 2.13: q_M = M 1_{A_M} e with |A_M| = M^{−2+ε} gives ‖q_M‖ = M^{ε/2} and F/‖q‖ = M^{p−2+ε/2}/p → 0 for ε < 2(2−p). With |A_M| = M^{−2} the midpoint gap (M^{p−2}/p)(1/2−2^{−p}) → 0. The splitting argument is valid because F finite ⇒ F₁ finite.
15. Lemma 2.14 (Hahn–Banach): the finite directional derivative, sublinearity, the extension dominated by Ψ'(q₀;·), and continuity from local boundedness.
16. Prop. 2.15(i): the multiplier via Lemma 2.14 with Ψ = F − (λ,·); density of H¹₀ in W^{1,p}₀; continuity of J on W^{1,p}₀ using W^{1,p}₀ ↪ L² (d = 2, or p ≥ 2d/(d+2)); identification with u_p by uniqueness. (ii) is also correct.
17. Lemmas 3.1 and 3.2.
18. The Cianchi–Maz'ya 2018 route to (S) (W^{1,2} ⊂ L^s in 2D, and b(t) ≥ a(1)t^{1+i_a}). It is correct modulo the flagged notion of solution.
19. Theorem C (i): the L¹–L² comparison with a_min; coercivity; strict convexity via strict monotonicity of γ. (ii): local boundedness of a finite convex function in finite dimension, and Lemma 2.14. (iii): weak = strong convergence in finite dimension.
20. Lemma 5.1.
21. Remark 5.2: the flux identification and d = 2 (modulo G1 and G5).

**Attributions consistent with REFERENCES A/B.**
- Lieberman (hedged), Antonini Cor. 1.4 (L^q, q > n, equals the REFERENCES' L^d, d > n).
- Cianchi–Maz'ya 2018 Thm 2.3; JEMS Thm 3.5.
- Fan (hedged).
- CM2015a Thm 1.2 (strict); BCM2018 Thms 1.1, 1.3, 1.4; CM2015b; Hästö–Ok.
- ELM and FMM (with the corrected meaning of "regularity").
- Zhikov 1995 Ex. 3.1.
- BDS Thms 30 and 32 (with the arXiv-numbering caveat); BGS Thm 2.3.
- Carstensen–Ortner; BOS Thm 1.
- Glowinski–Pan–Tai §3.2 hypotheses (i)–(iv) and the golden-ratio bound; ABR Thm 4.3.
- Del Pezzo–Martínez (H4) remark; Eckstein–Bertsekas Thm 8; Boyd §3.2.
- Bartels–Milicevic: D_j² = 2rE_n for τ = r.
- Bùi–Combettes counterexample; Moradifam–Nachman.
- Barrett–Liu, Liu–Barrett Thm 4.3, Diening–Kreuzer (1.4), Del Pezzo–Martínez 2015, Glowinski–Marrocco.

**Numerics recomputed from the raw data.** All of the following match the paper:
- Table 3 orders over the three finest meshes: p1.5 2.01/1.02/1.94; p2 2.00/1.02/2.02; p4 2.00/1.00/1.17; pq24 2.00/1.01/2.02; pq1.5_3 2.01/1.02/1.98; plog2_1 2.00/1.01/1.37; plog3_-1 2.00/1.02/2.02; exp 2.00/1.01/2.02; twophase 1.90/0.98/0.79; expdead 1.97/1.01/1.11.
- The pairwise nodal orders for t⁴/4 (0.91, 1.55, 1.81, 0.52).
- The finest-mesh errors and all iteration counts (Table 4).
- Table 2 (sc.txt): true and nodal orders 2.00/2.53, 2.00/2.08, 2.00/2.44, 2.00/2.46; true/nodal ratios 25.2, 6.9, 16.3, 18.1.
- Table 5 and §6.3: fits 1.98/1.03 (p = 1.2) and 1.84/1.02 (p = 1.1); u_max = 5.2·10⁻³ and 8.9·10⁻⁵ agree with the closed form ((p−1)/p)(c/2)^{1/(p−1)}.
- Table 6: pairwise orders 2.015, 2.004, 2.001 and 1.001, 1.000, 1.000.
- Table 7: both columns, against rman.txt and rdead.txt.
- Table 8 against lyapunov.json; the normalisation by E₀ was checked in the script. I re-ran the t²/2+t⁴/4 case and reproduced min cross = −4.79·10⁻²⁹ E₀ at n = 54, where E_n/E₀ = 1.9·10⁻²⁸, which confirms "at iterations where E_n itself is at that level".
- §6.8 against two_phase.txt.
- Fig. 1 highlights exactly the four reduced-order cases. Fig. 4 uses the best r (30), consistent with exp_dead_core.py.
- The dead-core radius is 2/c = 1/4.
- All non-converged runs are disclosed, apart from the count inconsistency in E10.
