# Referee report on `convergence_note.tex`

*Convergence of the decomposition–coordination method for the g-Laplacian at the continuous level (working note, draft)*

## Overall verdict

The core of the note is correct. That covers Proposition 1 (⇐), the well-posedness of the steps, (alg.u2), the Lyapunov identity (lyap) with both sign claims, Step 3 (maximal monotone closure), the Opial step, and Theorem B(a) and B(b). With minor corrections it is a sound, self-contained proof of:

> weak convergence of ADMM/ALG2 (ρ = r) under convexity plus a saddle point, and strong convergence under the uniform-monotonicity conditions (a) and (b).

The problems are of three kinds:

1. **One genuine (easily fixed) error in a proof:** the (⇒) direction of Proposition 1 uses a test pair that does not give the multiplier equation.
2. **Several overstated or incorrect statements in the discussion parts:**
   - Remark 1 has the finiteness of F backwards. It also misdescribes what (S) means for p < 2 versus p > 2.
   - Remark 2 claims that a discrete saddle point exists "for every convex Φ", which is false without coercivity.
   - Remark 4 asserts that a statement is false without giving a counterexample.
   - The paragraph after Assumption (S) misstates the data needed for Lipschitz regularity.
3. **Small imprecisions:**
   - ∂F(q) = {T(·,q)} should be an inclusion.
   - The condition on (ν₀, η₁) is unnecessary.
   - The pointwise prox formula fails in the nonsmooth case.
   - Some claims are redundant.

None of these affects the main theorems once corrected.

---

## 1. Confirmed steps (checked line by line)

**Setting.**
- F is convex, since t ↦ Φ(x,t) is convex and nondecreasing and |·| is convex.
- F is measurable, since Φ is Carathéodory.
- F is proper, since F(0) = 0.
- F is lsc on H. Take a subsequence realising the liminf, then a further subsequence converging a.e. Continuity of Φ(x,·) (given by Carathéodory, and in any case automatic for a convex nondecreasing finite function on [0,∞)) and Fatou then give the result.

**Proposition 1 (⇐).** Correct. Since ∇u* = ν*, the term (η*, ∇v − q − ∇u* + ν*) equals (η*, ∇v − q), and the computation closes.

**Proposition 1 (⇒), the other parts.**
- The left inequality forces ∇u* = ν*.
- The second part (v = u*, q arbitrary, then q ↦ ν* + s(q − ν*), convexity, s → 0) is correct.
- F(ν*) < ∞ is automatic: if it were infinite, the right inequality at (v, q) = (0, 0) would fail.
- "u* minimises (P)" follows from (S): F(∇v) − F(∇u*) ≥ (η*, ∇(v − u*)) = (f, v − u*).
- The independence of r is clear from (S).

**Well-posedness of the steps.**
- (alg.u) is the Euler–Lagrange equation of v ↦ 𝓛_r(v, ν_{n−1}, η_n). Lax–Milgram applies because the right-hand side (f,w) − (η_n − rν_{n−1}, ∇w) lies in V*.
- (alg.nu): the functional is proper, lsc and r-strongly convex, so it has a unique minimiser.
- The sum rule used for the optimality condition is legitimate, because the quadratic term is continuous on H.
- Hence η_{n+1} ∈ ∂F(ν_n) for n ≥ 1, and for n = 0 by the choice of η₁.

**(alg.u2).** η_n + r(∇u_n − ν_{n−1}) = η_{n+1} + r(ν_n − ν_{n−1}). Correct.

**(lyap), recomputed by hand.**
- Step 1 expansion: with w = u_n − u* we get ∇w = s_n + ν_n − ν*, and the four-term expansion is correct.
- First term: (η_{n+1} − η*, s_n) = (1/2r)(‖η_{n+1} − η*‖² − ‖η_n − η*‖²) + (r/2)‖s_n‖². Correct, using r s_n = η_{n+1} − η_n.
- Second term: r(ν_n − ν_{n−1}, ν_n − ν*) = (r/2)(‖ν_n − ν*‖² − ‖ν_{n−1} − ν*‖²) + (r/2)‖ν_n − ν_{n−1}‖². Correct.
- Collecting terms gives exactly (lyap), with E_n and E_{n−1} as defined; E₀ uses (ν₀, η₁).
- Note that the identity is purely algebraic. It uses only (alg.u2), ν* = ∇u* and (η*, ∇w) = (f, w). The subgradient property is used only for the signs.

**Signs.**
- m_n ≥ 0 by monotonicity of ∂F.
- r(s_n, ν_n − ν_{n−1}) = (η_{n+1} − η_n, ν_n − ν_{n−1}) ≥ 0 by monotonicity, since η_{n+1} ∈ ∂F(ν_n) and η_n ∈ ∂F(ν_{n−1}).
- The summability, boundedness and residual statements follow.

**Step 3.**
- ∇ is weakly continuous V → H, so ∇ū = ν̄.
- Passing to the limit in (alg.u2) gives (η̄, ∇w) = (f, w).
- The identity (η_{n+1}, ν_n) = (f, u_n) − (η_{n+1}, s_n) − r(ν_n − ν_{n−1}, ∇u_n) is correct, and its limit along the subsequence is (η̄, ν̄).
- The maximal-monotone closure argument is correct: Rockafellar's theorem applies because F is proper, convex and lsc on a Hilbert space.
- On ν_{n_k−1}: it is not needed, but it has the same weak limit because ν_n − ν_{n−1} → 0 strongly.

**Step 4 (Opial).**
- |x_n − z|² = 2E_n.
- The "two cluster points" argument is correct.
- Every weak cluster point of x_n is a saddle pair, after extracting a further subsequence for u_n.
- η_n ⇀ η̄ also holds, since η_{n+1} − η_n → 0.
- u_n ⇀ ū in V is immediate, because ∇ is an isometry onto a closed subspace.

**Theorem B(a).**
- The identity ψ = ½(h_a + h_b)|a−b|² + ½(h_a − h_b)(|a|² − |b|²) is correct; I expanded it by hand and checked it numerically.
- The second term is ≥ 0 because h is nondecreasing.
- ψ ≥ ½ h(max{|a|,|b|})|a−b|² ≥ ½ h(|a−b|/2)|a−b|² is correct.
- On {|e_n| > δ}, monotonicity of h(x,·) together with the hypothesis h(x, δ/2) ≥ h₀(δ/2) gives the bound. **There is no issue here:** the hypothesis "h(x,t) ≥ h₀(t) > 0 for all t > 0" is exactly what is used, at t = δ/2.
- The splitting estimate ∫|e_n|² ≤ δ²|Ω| + (2/h₀(δ/2)) m_n is correct, because ψ ≥ 0 everywhere.
- h₀ need not be measurable. One may take h₀(t) = ess inf_x h(x,t), which is automatically nondecreasing.
- All listed examples satisfy the hypotheses:
  - t^p/p with p ≥ 2;
  - t^p/p + a t^q/q with p, q ≥ 2 and a ≥ 0;
  - t² log(1+t), for which h = 2 log(1+t) + t/(1+t);
  - e^t − 1 − t, for which h = (e^t − 1)/t ≥ 1.

**Theorem B(b).**
- The inequality ψ ≥ (p−1)(|a|+|b|)^{p−2}|a−b|² for 1 < p ≤ 2 is correct. Proof: write T(a) − T(b) = ∫₀¹ DT(b + s(a−b))(a−b) ds. The smallest eigenvalue of DT(ξ) is (p−1)|ξ|^{p−2}, and |b + s(a−b)| ≤ |a| + |b| with p − 2 < 0. The singularity at 0 is integrable because p > 1. A random test (2·10⁵ pairs over 12 orders of magnitude) gives min ψ / bound ≥ 1.07.
- The Hölder step with exponents 2/p and 2/(2−p) is correct.
- The final rate is ‖e_n‖_p^p ≤ C (m_n/(p−1))^{p/2}. The factor (p−1)^{−p/2} goes into C.
- Strong W^{1,p}₀ convergence follows, since s_n → 0 in L² ⊂ L^p.

---

## 2. Errors and gaps, with corrections

### E1 (error in a proof). Proposition 1, (⇒): derivation of (η*, ∇w) = (f, w)

**Problem.** With (v, q) = (u* + sw, ν* + s∇w) one has ∇v − q = 0. Both the η*-term and the penalty then vanish, and the right inequality only gives F(ν* + s∇w) − s(f, w) ≥ F(ν*). That is the minimality of u*, not the multiplier equation (and the left side may be +∞).

**Correction.** Move v only and keep q = ν*:

```latex
The right inequality with $(v,q)=(u^*+sw,\nu^*)$, $s\in\R$, gives
$s\big[(\eta^*,\nabla w)-(f,w)\big]+\frac{r s^2}{2}\|\nabla w\|^2\ge0$ for all $s$,
hence $(\eta^*,\nabla w)=(f,w)$ (divide by $s>0$ and $s<0$ and let $s\to0$).
```

Also state that F(ν*) < ∞, from the right inequality at (0, 0).

### E2 (incorrect statements). Remark 1 ("What this says about p < 2") and the last sentence of the (S) paragraph

**(i) Finiteness is backwards.**
- For Φ = t^p/p with p ≤ 2, F is finite everywhere on H (|q|^p ≤ 1 + |q|²) and therefore continuous.
- For p > 2, F = +∞ on L² \ L^p.
- So "p ≥ 2 … guarantees a finite Lagrangian" is wrong: it is p ≤ 2 that gives a finite and continuous F on H.
- What p ≥ 2 does give is that J = +∞ off W^{1,p}₀ ⊂ H¹₀. Hence the W^{1,p}₀ minimiser exists, lies in V, and the H¹₀ problem is attained.

**(ii) What (S) means depends on p.**
- **p < 2:** if u ∈ H¹₀, then |η| = |∇u|^{p−1} ∈ L^{2/(p−1)} ⊂ L², because 2/(p−1) > 2 and Ω is bounded. F is continuous, so the sum rule ∂(F∘∇) = ∇*∂F∘∇ holds. Therefore, for p < 2, **(S) ⇔ the W^{1,p}₀ minimiser belongs to H¹₀**. The alternative "or |∇u|^{p−1} ∉ L²" in Remark 1 is vacuous.
- **p > 2:** (S) ⇔ |∇u|^{p−1} ∈ L², i.e. ∇u ∈ L^{2p−2}. This is a genuine higher-integrability requirement, not implied by u ∈ W^{1,p}₀ (which only gives η ∈ L^{p'}, with p' < 2). So (S) is **not automatic for p ≥ 2** either; it is equally an assumption in Glowinski's theorem.

**Suggested replacement** for the relevant sentences:

```latex
For $\Phi=t^p/p$, (S) means: the minimiser $u\in W^{1,p}_0$ satisfies $u\in H^1_0$ and
$|\nabla u|^{p-1}\in L^2$. If $p<2$ the second condition follows from the first
($|\nabla u|^{p-1}\in L^{2/(p-1)}\subset L^2$), and $F$ is finite and continuous on $H$;
if $p>2$ the first is automatic ($W^{1,p}_0\subset H^1_0$) and the second is a genuine
regularity requirement ($\nabla u\in L^{2p-2}$), since $F=+\infty$ on $L^2\setminus L^p$.
Thus neither range of $p$ is privileged: what matters is (S).
```

**(iii) Realism of (S).** Global regularity results of Cianchi–Maz'ya type state the second-order estimate a(|∇u|)∇u ∈ W^{1,2}(Ω) for f ∈ L², 1 < p < ∞ (and for more general Φ), on convex domains or domains with suitably regular boundary. Such a result would give η* ∈ W^{1,2} ⊂ L². For p < 2 it would also give ∇u ∈ L^s for every s < ∞ in 2D, hence u ∈ H¹₀. If that holds, **(S) is satisfied for all p > 1 on convex polygons with f ∈ L²**.

I believe this is in Cianchi & Maz'ya (ARMA 2018, "Second-order two-sided estimates in nonlinear elliptic problems"), but the exact hypotheses must be checked before citing. This would be the natural way to make Theorem A unconditional on convex domains.

**(iv) Open case in Remark 1.** When (S) fails, ADMM with ρ = r is Douglas–Rachford applied to the dual, and the fixed points correspond to saddle points. For averaged maps with no fixed point, ‖z_n‖ → ∞ (Pazy 1971). This suggests ‖η_n‖ → ∞ in that case, so a different argument would be needed. This should be stated as a remark with a reference to check, not as a claim.

### E3 (incorrect statement). Remark 2: "a saddle point always exists for every convex Φ"

**Problem.** Continuity of F_h gives the multiplier once a discrete minimiser exists. Existence of the minimiser requires coercivity, and convexity is not enough:
- For Φ(t) = t (convex, nondecreasing, Φ(0) = 0, admitted by the setting), J_h(v) = Σ|T||∇v_T| − (f,v)_h is unbounded below when f is large.
- For Φ ≡ 0 and f ≠ 0, J_h is unbounded below for every mesh.

**Correction.**

```latex
If $\Phi(x_T,t)/t\to\infty$ as $t\to\infty$ for every $T$ (true for $t^p/p$, $p>1$,
and for $e^t-1$, $e^t-1-t$), then $J_h$ is coercive on $V_h$ ($\nabla$ is injective),
a minimiser exists, and since $F_h$ is finite (hence continuous) the sum rule gives
$\eta_h\in\partial F_h(\nabla u_h)$ with $\nabla_h^{*}\eta_h=f_h$: (S) holds.
```

Replace "every convex Φ" with "every convex Φ with superlinear growth (more generally, whenever the discrete problem has a minimiser)". Also state explicitly that the discrete inner products (weighted by |T|, together with the discretised (f,v)) must be the ones used in the algorithm. The proof then applies verbatim.

### E4 (overstated). Remark 4: "the generic statement … is not correct for a general convex F"

**Problem.** No counterexample is given. The general theory only proves weak convergence, but that does not show the strong-convergence claim is false within the class "G linear, B injective with closed range, saddle point exists".
- The known counterexample to strong convergence (Güler 1991) concerns the proximal point algorithm.
- ADMM is a proximal point method for a maximal monotone operator that is not, in general, a subdifferential.

So the claim is plausible, but unproved here.

**Correction.**

```latex
Our argument yields only weak convergence for a general convex $F$; strong convergence
is obtained under a uniform-monotonicity property of $\partial F$ (Theorem~\ref{thm.B}).
We do not know whether the statement ``$B$ injective with closed range + saddle point
$\Rightarrow$ strong convergence'' holds for general convex $F$; the general ADMM/
Douglas--Rachford theory guarantees only weak convergence in infinite dimensions.
```

### G1. ∂F(q) = {T(·,q)} at the start of Section 4

This should read: ∂F(q) = {T(·,q)} if T(·,q) ∈ H, and ∂F(q) = ∅ otherwise. This is Rockafellar's characterisation of the subdifferential of integral functionals on L². It applies here because Φ(x,|·|) is a normal convex integrand, F is proper, and the conjugate integral is proper (Φ*(x,0) = 0). Cite it.

Alternatively, η_{n+1} = T(·,ν_n) follows directly from the pointwise prox. For η* ∈ ∂F(ν*), the pointwise identification is needed, so the citation is required.

### G2. The paragraph after Assumption (S)

- "g(x,·) locally bounded" should be **uniform in x**, e.g. ess sup_x g(x, M) < ∞ for every M (or g(·, ‖∇u‖_∞) ∈ L²).
- One must also say that η* ∈ ∂F(∇u) follows from the pointwise subgradient inequality, integrated. That needs the C¹ case or a measurable selection.
- Global Lipschitz regularity for the p-Laplacian on convex domains requires more than f ∈ L² in 2D (for example f ∈ L^q with q > 2, or the Lorentz space L^{2,1}). Fix the phrase "with smooth data", or better use the second-order regularity route in E2(iii), which needs only f ∈ L².

### G3. Initialisation

- Because F ≥ 0 = F(0), we have 0 ∈ ∂F(0) **always**. So ν₀ = η₁ = 0 is admissible for every Φ, including e^t − 1, with no differentiability assumption; delete the parenthetical condition.
- The condition η₁ ∈ ∂F(ν₀) is not needed for convergence anyway. Without it, (lyap) still holds and only the n = 1 cross term may be negative, which changes E₁ ≤ E₀ into E₁ ≤ E₀ + |r(s₁, ν₁ − ν₀)|. All conclusions are unaffected.

### G4. Pointwise formula for (alg.nu)

The formula "g(x,t) + rt = r|w|" assumes Φ(x,·) ∈ C¹ with g(x,0) = 0. In general:
- t = 0 if r|w| ≤ ∂⁺_tΦ(x,0);
- otherwise r|w| ∈ ∂_tΦ(x,t) + rt.

For e^t − 1 the threshold is ∂⁺_tΦ(x,0) = 1. Measurability of ν_n follows from the measurability of the prox of a Carathéodory integrand; alternatively, ν_n is the unique H-minimiser, so the pointwise minimiser coincides with it a.e.

### G5. Theorem B: redundant claims and possible strengthenings

- "η_{n+1} − η_n → 0 strongly" in (a) is already part of Theorem A. It may have been meant as η_n → η̄.
- **For (b), strong convergence of the multipliers is available.** For 1 < p ≤ 2 the map T is (p−1)-Hölder: |T(a) − T(b)| ≤ 2^{2−p}|a − b|^{p−1}. So

  ```latex
  \|\eta_{n+1}-\bar\eta\|^2\le C\int_\Omega|\nu_n-\bar\nu|^{2(p-1)}\,dx\to0,
  ```

  because 2(p−1) ≤ p and ν_n → ν̄ in L^p on a bounded Ω. Hence η_n → η̄ strongly in H.
- Also by interpolation (L^p convergence plus L² boundedness), ν_n → ν̄ and ∇u_n → ∇ū in L^q for every q < 2.
- In (a), η_n → η̄ strongly holds when the Nemytskii map T is continuous L² → L², e.g. g(x,t) ≤ C(1+t). It is not available in general for p > 2.
- On "unique if Φ strictly convex": in case (a) uniqueness of ν̄ (hence of ū) is automatic, since the proof gives ν_n → ν* for every saddle point. In addition, h > 0 and h nondecreasing imply that g is strictly increasing.
- Optional: the factor ½ in ψ ≥ ½h(|a−b|/2)|a−b|² can be dropped. In the numerical check the ratio ψ / (½h(|a−b|/2)|a−b|²) never fell below 2, with equality approached for a = −b. Not needed for the argument.

### G6. Minor points

- Step 3: Rellich is not needed, because (f, u_{n_k}) → (f, ū) already follows from weak convergence, f ∈ L² ⊂ V*. Where Rellich is used (L² convergence in Theorem A), it holds for H¹₀ on any bounded open set (extend by zero); the Lipschitz assumption is not needed for it.
- Identity in (a) at a = 0: set h(x,0) := h(x,0+). The formula holds for any value of h(x,0), but monotonicity of the second term requires h(x,0) ≤ h(x,|b|).
- In (P), the "i.e." PDE form is formal: it needs Φ ∈ C¹, and without (S) the minimum over V may not be attained (p < 2).
- Remark 3: the condition 0 < ρ < (1+√5)/2 · r is the Fortin–Glowinski / Gabay condition. Fine, but give a reference.
- Status paragraph and bibliography: the hypotheses and numbering of Glowinski's theorem remain unverified, as the note itself says. Keep that caveat.

---

## 3. Numerical sanity check (finite-dimensional ADMM)

Script: `/tmp/claude-0/-home-claude/13ce5d30-34f9-5095-a3f8-3470cf425272/scratchpad/admm_check.py`

**Setup.**
- V = R⁶, H = (R²)¹⁰, B a random 20×6 matrix (injective), r = 1.3.
- F(q) = Σ_i Φ(|q_i|) on 2D blocks, with Φ ∈ {t^1.5/1.5, t^1.2/1.2, t²/2, t³/3, e^t − 1 − t, e^t − 1}.
- ν₀ = η₁ = 0. The ν-step is an exact blockwise prox, solved by bisection, including the threshold at 0 for e^t − 1.
- The saddle point is taken as the limit after 4000 iterations, made exactly feasible (ν* = Bu*, Bᵀη* = f).

**Results over the first 200 iterations.**

| Φ | max relative residual of (lyap) | min m_n | min cross term |
|---|---|---|---|
| p = 1.5 | 1e-17 | ≥ 0 | ≥ 0 |
| p = 1.2 | 1e-17 | 4e-5 | ≥ 0 |
| p = 2 | 6e-17 | ≥ −4e-32 | ≥ −4e-32 |
| p = 3 | 2e-16 | ≥ 0 | ≥ −2e-31 |
| e^t − 1 − t | 2e-18 | ≥ −1e-32 | ≥ −4e-32 |
| e^t − 1 (mixed zero / nonzero blocks in ν*, ∂F multivalued) | 2e-16 | 2e-14 | 9e-17 |

- The negative values are at the 1e-32 roundoff level.
- The residual of (alg.u2) is ≤ 5e-15 and the final primal residual is ≤ 1e-9.
- So (lyap) holds to machine precision, and both sign claims hold, including in the multivalued case.

Script: `/tmp/claude-0/-home-claude/13ce5d30-34f9-5095-a3f8-3470cf425272/scratchpad/ineq.py`

This also checked, on random pairs:
- the identity for ψ in (a), with h = (e^t − 1)/t (error 6e-14);
- the bound ψ ≥ ½h(|a−b|/2)|a−b|² (smallest ratio 2.0);
- the p < 2 inequality (smallest ratio 1.07 at p = 1.9).

---

## 4. Summary

| Item | Status |
|---|---|
| Prop. 1 (⇐) | correct |
| Prop. 1 (⇒) | **error in the multiplier step (E1)**; easy fix |
| Well-posedness, η_{n+1} ∈ ∂F(ν_n), (alg.u2) | correct |
| (lyap), signs | correct (verified by hand and numerically); η₁ ∈ ∂F(ν₀) is not needed and is always satisfiable with 0 (G3) |
| Step 3 | correct (Rellich superfluous) |
| Opial | correct |
| Thm B(a) | correct; the h₀ uniformity is exactly what is used |
| Thm B(b) | correct; can add η_n → η̄ strongly and L^q convergence for q < 2 (G5) |
| ∂F = {T} | should be an inclusion, with a citation to Rockafellar (G1) |
| (S) paragraph | data/uniformity imprecise (G2); a better route via second-order regularity (E2 iii) |
| Remark 1 | **finiteness backwards; misstates (S) for p < 2 vs p > 2 (E2)** |
| Remark 2 | **needs coercivity (E3)** |
| Remark 4 | **unproved claim of falsity (E4)** |

**Verdict:** the main theorems (A, B(a), B(b)) are correct after fixing E1, and they can be built on. Remarks 1, 2 and 4 must be rewritten as indicated before the note is circulated.
