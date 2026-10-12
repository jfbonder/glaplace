# Mathematical referee report

**Paper:** Aragón & Fernández Bonder, *Convergence of the Glowinski decomposition–coordination method for the g-Laplacian…* (`/home/claude/paper_v3/new/main.tex`, version of 12 Oct, 1020 lines; I did not modify it).
**Primary sources used:** the PDFs in `/mnt/user-data/uploads/Repos/glaplace/biblio/`, read in full where relevant: Glowinski 1984 Ch. VI, GPT16 (cam16-10), Lieberman 1991, Cianchi–Maz'ya 2018, De Filippis–Mingione, Colombo–Mingione 2015a/b, BCM 2018, Esposito–Leonetti–Mingione, BDS, BGS, BOS, Carstensen–Ortner, Del Pezzo–Martínez, Eckstein–Bertsekas, Bùi–Combettes, Hästö–Ok.
Line numbers refer to `main.tex`.

---

## 0. Verdict

The proofs are correct: I found **no fatal error**. I checked every proof line by line (§2 below). The problems are in **how the paper places its results relative to Glowinski's book**, which is now available, and in a few literature statements.

1. **Glowinski's Theorem 5.1 applies much more widely than §2.6, the introduction and Table 3 say.** Hypothesis (iv) holds *literally* for p ≥ 2, through an explicit splitting. It also holds for e^t−1, e^t−1−t, the two-phase energies, t^{1.5}/1.5+t^3/3, t^3/log(1+t) and the manufactured problem, with F₀ = (ε/2)‖·‖². Under (S), the book therefore already gives **strong** convergence of (u_n, ν_n) in H¹₀ × L² for every test problem in Table 3 **except t^{1.5}/1.5**, and for all 0 < ρ < (1+√5)r/2. Theorem B(a) is essentially a special case. (MAJOR)
2. **Attribution.** Two things in the paper are presented as derived from later literature but come from Glowinski Ch. VI:
   - the Lyapunov functional E_n and the inequality version of (2.10) are Glowinski's (5.25) with ρ = r;
   - the Opial-type weak convergence of the whole multiplier sequence is stated in Ch. VI §7.

   In addition, the existence of a discrete saddle point (Theorem C(ii)) is stated in Ch. VI §4.2(1), with a reference to Rockafellar. The "new points" paragraph (line 90) must be revised. (MAJOR)
3. **Lieberman 1991, Thm 1.7 is an interior result.** The boundary statement is only a remark that points to Lieberman (1988). (MAJOR for the Table / Literature statement 3.3)
4. **"Convergence when (S) fails is open" (lines 93, 449, 843) is too weak.** Douglas–Rachford/ADMM theory in Hilbert space (Eckstein–Bertsekas Thm 7) gives: if (S) fails, the iterates η_n + rν_{n−1} are **unbounded**. Bùi–Combettes, transferred to ADMM, shows that weak convergence of η_n is sharp under (A1)–(A3). (MAJOR: a stated open problem is partly settled by known theory.)

The other issues are minor, mostly hypotheses missing from literature statements. Of the 54 TODO items, 27 are confirmed in full or in part, 3 are wrong (with corrections), and 24 could not be checked from the available PDFs. One partially confirmed item (#13) also records a misquote of Glowinski by Del Pezzo–Martínez. Section 4 lists all items.

---

## 1. Issues (most severe first)

### M1 (MAJOR): §2.6, line 79, line 84, line 399, lines 422–424, Table 3 "theory" column, §7 (line 841). What Glowinski's theorem actually covers.

**What the book says** (Glowinski 1984, Ch. VI, pp. 171–182, now available):
- Problem (P): min F(Bv) + G(v), with V, H Hilbert spaces and F, G convex, proper and l.s.c.
- Hypotheses of §4.1:
  - (4.1) dom(F∘B) ∩ dom G ≠ ∅;
  - (4.2) B injective with R(B) closed;
  - (4.3) lim_{|q|→∞} F(q)/|q| = +∞;
  - (4.4)–(4.5) F = F₀ + F₁, with F₀ and F₁ convex, proper and l.s.c., and F₀ Gâteaux-differentiable and "uniformly convex on the bounded sets of H". This is **defined** in (4.6) as: for every M there is δ_M continuous and strictly increasing with δ_M(0) = 0 such that (F₀′(q) − F₀′(p), q − p) ≥ δ_M(|q − p|) for |p|, |q| ≤ M.
- Thm 4.1 (ALG1) and **Thm 5.1 (ALG2)**: if 𝓛_r has a saddle point and 0 < ρ < (1+√5)r/2, then u^n → u **strongly** in V, p^n → p strongly in H, λ^{n+1} − λ^n → 0, (λ^n) is bounded, and weak cluster points of (λ^n) are multipliers.
- §7: by Opial, the whole sequence λ^n converges weakly.
- ALG2 (3.7)–(3.10) has the u-step first, the same order as in the paper. ({p⁰, λ¹} is arbitrary.)

These statements of the hypotheses and conclusions in §2.6 are **confirmed**: (i)–(iv), "strongly", the numbering and the golden-ratio bound.

**What is wrong in the paper:**

(a) **p > 2 (lines 79 and 422).** The paper says (iv) "must be understood in a weaker sense" because there is "no obvious splitting". There is one. Take h(t) = t^p/p for t ≤ 1 and h(t) = 1/p + (t²−1)/2 for t > 1. In general, take h = G on [0,1] and h(x,t) = G(x,1) + ½k(x,1)(t²−1) for t > 1. Then F₀ = ∫h(|q|) is C¹ on L² with Lipschitz gradient and uniformly monotone, F₁ = F − F₀ is convex, and (4.3) holds. The general statement is part (b) of the lemma below. So the book's Thm 5.1 applies literally to t^p/p with p ≥ 2 in H¹₀ × L². It also applies to p(x) ≥ 2 with p₊ < ∞, so the conclusion of [AFR23, Thm 4.3] is justified, at least through this splitting.

(b) **Uniformly convex energies.** Whenever G(x,·) − (ε/2)t² is convex for some ε > 0, take F₀ = (ε/2)‖q‖² and F₁ = F − F₀. Glowinski allows F₁ to be nonsmooth and to take the value +∞. This covers:
- e^t − 1 (g′ = e^t ≥ 1; ∂F₁ is multivalued);
- e^t − 1 − t;
- t²/2 + a(x)t⁴/4 for any measurable a ≥ 0 (two-phase and manufactured);
- t^{1.5}/1.5 + t³/3 (g′ = ½t^{−1/2} + 2t ≥ 3/2);
- t³/log(1+t) (g′ ≥ 2, checked numerically on [10⁻⁸, 10⁸]).

So these sentences are **wrong**:
- line 424: "Theorem A requires neither (iii) nor the differentiability in (iv), so it covers e^t−1" (the classical theorem covers it too, with *strong* convergence);
- Table 3: "A" for pq1.5_3 and for expdead.

(c) **Theorem B(a)** follows from the book's theorem whenever ess sup_x k(x,1) < ∞, by part (b) of the lemma below, and then for all ρ < (1+√5)r/2. The only cases of B(a) not covered by the classical theorem are those with ess sup_x k(x,1) = ∞ and inf_t g′(x,t) = 0, for example t^p + a(x)t^q with p > 2 and a unbounded. The new parts of the continuous theory are therefore Theorem B(b), and Theorem A (weak convergence) for G with no quadratic lower convexity at all (e.g. G with affine pieces, or t^p with p < 2).

(d) Since ℓ is linear here, the book's Remark 5.1 and §6.3 also record the larger range 0 < ρ < 2r, due to Gabay–Mercier. This is worth one sentence in place of "ρ ≠ r … (1+√5)/2" (line 424).

(e) Book §6.3 (p. 186) says explicitly that for the s-Laplacian (with V = W^{1,s}₀, H = L^s) Theorems 4.1/5.1 "cannot be applied" for s ≠ 2 because the spaces are not Hilbert, and that they apply after discretisation (Thms 4.2/5.2). This is useful context for the paper and should be cited.

**Fix.** Rewrite §2.6 from the book, delete "We have not been able to consult", and add:

```latex
\begin{lemma}[Glowinski's hypotheses (iii)--(iv) on $H=L^2$]\label{lem:split}
Let $G$ be a Young integrand and $F(q)=\int_\Omega G(x,\abs q)\dx$.
\begin{enumerate}[label=(\alph*)]
\item If for some $\varepsilon>0$ the function $t\mapsto G(x,t)-\frac\varepsilon2t^2$ is convex on $[0,\infty)$ for a.e.\ $x$,
then $F_0=\frac\varepsilon2\norm{\cdot}^2$, $F_1=F-F_0$ satisfy (iii)--(iv), with $\delta_M(s)=\varepsilon s^2$.
\item If $G(x,\cdot)\in C^1$, $g(x,0)=0$, $k(x,\cdot)=g(x,\cdot)/t$ is nondecreasing,
$k_*(t):=\operatorname*{ess\,inf}_xk(x,t)>0$ for $t>0$ and $K:=\operatorname*{ess\,sup}_xk(x,1)<\infty$, let
$h(x,t)=G(x,t)$ for $t\le1$ and $h(x,t)=G(x,1)+\frac{k(x,1)}2(t^2-1)$ for $t>1$.
Then $F_0(q)=\int_\Omega h(x,\abs q)\dx$ and $F_1=F-F_0$ satisfy (iii)--(iv): $F_0\in C^1(H)$, and
$(F_0'(q)-F_0'(p),q-p)\ge\frac38\,k_*\!\big(\min\{\tfrac{\norm{q-p}}{4\abs\Omega^{1/2}},1\}\big)\norm{q-p}^2$.
\end{enumerate}
\end{lemma}
```

Proof sketch. (a) φ₁ = G − εt²/2 is convex with φ₁′(0+) = g(x,0) ≥ 0, so it is nondecreasing; hence φ₁(x,|ξ|) is convex and nonnegative, and F₁ is l.s.c. by Fatou. F ≥ F₀ gives (iii).
(b) h′(x,t) = t·k(x, min{t,1}) ≤ Kt, so F₀′ is a Lipschitz Nemytskii map on L². G − h is convex and nondecreasing because, for t > s ≥ 1, g(t) − g(s) ≥ (t−s)k(x,s) ≥ (t−s)k(x,1). The monotonicity bound is Lemma 2.9(i) applied to k(x, min{t,1}), followed by the δ-splitting of the proof of Theorem B(a) with δ = ‖e‖/(2|Ω|^{1/2}). For (iii), G ≥ h ≥ k_*(1)(t²−1)/2 when t ≥ 1.

Then: change the Table 3 "theory" column to "G5.1" for every row except p1.5 (B(b)); state that (S) is still needed; and rephrase the abstract, the contributions (items 2 and 4) and §7 accordingly.

### M2 (MAJOR): line 90 ("new points"), Remark 2.13 [rem:A-literature], Theorem C(ii). Attribution.

- **Lyapunov functional.** Glowinski (5.25) with ρ = r reads
(|λ̄ⁿ|² + r²|p̄^{n−1}|²) − (|λ̄^{n+1}|² + r²|p̄ⁿ|²) ≥ 2r(F₀′(pⁿ) − F₀′(p), p̄ⁿ) + r²|Buⁿ − pⁿ|² + r²|pⁿ − p^{n−1}|².
Dividing by 2r, this is exactly (2.10) as an inequality, with E_n the paper's functional; the book derives it via (5.17)–(5.24), using the cross-term sign from (5.10) at n−1. Remark 2.13 credits Bartels–Milicevic and Boyd et al. It must credit Glowinski (5.25) / Fortin–Glowinski first.
- **Whole-sequence weak convergence via Opial.** This is stated in Ch. VI §7, p. 194 ("using some results due to Opial … weak convergence … of the whole sequence {λⁿ}"; the proof is referred to Glowinski–Lions–Trémolières, App. 2).
- **Theorem C(ii).** Ch. VI §4.2(1), p. 174: "Since the constraint Bv − q = 0 is linear, if (P) has a solution, then 𝓛 and 𝓛_r have a saddle point (see Rockafellar [1]…)". So the existence of discrete saddle points is classical finite-dimensional convex duality. The proof via the paper's Lemma 2.15 is fine, but it is not new.

**Fix.** In line 90, remove "the existence of discrete saddle points without any regularity (Theorem C(i)–(ii))" from the list of new points. Say that Theorems A and C follow Glowinski's argument (eq. (5.25)) with hypotheses (4.3)–(4.6) removed. The genuinely new items are Theorem B(b), Proposition 2.16, Remark 5.2 and Proposition 5.3 (see m11 for the last).

### M3 (MAJOR): Literature statement 3.3 (line 479) and Table 1 row 1 (line 522). Lieberman Thm 1.7.

Lieberman 1991, p. 320: Theorem 1.7 assumes (1.10a–d) for 0 < α ≤ 1, for x, y ∈ Ω, z, w ∈ [−M₀, M₀] and p ∈ ℝⁿ. Concretely:
- a^{ij}ξ_iξ_j ≥ (g(|p|)/|p|)|ξ|² and |a^{ij}| ≤ Λg(|p|)/|p|;
- |A(x,z,p) − A(y,w,p)| ≤ Λ₁(1 + g(|p|))[|x−y|^α + |z−w|^α];
- |B| ≤ Λ₁(1 + g(|p|)|p|).

The structure function g is C¹ with δ ≤ tg′/g ≤ g₀ (that is (1.1)′; Theorem 1.7 depends on δ). The conclusion is that every W^{1,G}(Ω) solution with |u| ≤ M₀ is in **C^{1,β}(Ω)**, with the estimate (1.11) on Ω′ ⋐ Ω: this is **interior**. The text after the theorem says only: "Again boundary Hölder gradient estimates for Dirichlet and conormal boundary conditions are proved similarly and we omit their statements and proofs except to point out that Theorems 1 and 2 of [15] remain valid with our conditions (1.10a–d)". Reference [15] is G. M. Lieberman, *Boundary regularity for solutions of degenerate elliptic equations*, Nonlinear Anal. 12 (1988).

**Fix.** "Lieberman [Thm 1.7] proves interior C^{1,β} regularity for bounded weak solutions; the corresponding global result under (1.10) is asserted in [Lieberman 1991, p. 320] by reference to [Lieberman 1988, Thms 1–2]." Add Lieberman 1988 to the bibliography and keep Antonini (Cor. 1.4) as the verified global statement. Also mention that Lieberman allows u-dependence and requires |u| ≤ M₀; for f ∈ L^∞ this is satisfied by the L^∞ bound.

### M4 (MAJOR): lines 93 (i), 449, 843 (open problem 1), Remark 2.12 [rem:weak-only]. What happens without (S), and sharpness of weak convergence.

(a) For ρ = r and the u-step first, the iteration is the ADMM of Eckstein–Bertsekas Thm 8 with f = −ℓ, M = B and g = F. The identification z_n := η_n + rν_{n−1}, η_n = J_{r∂F*}(z_n), see E–B (Y1)–(Y3), is algebraic and holds in Hilbert space. The relevant DR operators are A = ∂(ι_{\{B*η=ℓ\}}) = N_{\{B*η=ℓ\}}, whose subdifferential needs no qualification condition because R(B) is closed, and ∂F*. One checks directly:

> zer(A + ∂F*) ≠ ∅ ⇔ (S).

(If a ∈ ∂F*(η) and −a ∈ R(B), then a = Bu, η ∈ ∂F(Bu) and B*η = ℓ.) By E–B Thm 7, which is stated for a general Hilbert space, "If zer(A+B) = ∅, then {z^k} is unbounded".

Hence: **if (S) fails, (η_n + rν_{n−1}) is unbounded in L², so (η_n, ν_n) cannot converge or even stay bounded in L² × L².** Together with Theorem A: (S) ⇔ the iterates are bounded. This should replace "convergence … is open" by a precise statement, for example "(η_n, ν_n) is unbounded; the behaviour of u_n (cf. GPT16, Remark 23) is open". Suggested proposition:

```latex
\begin{proposition}
Assume \ref{A1}--\ref{A3} and $\rho=r$. Then (S) holds if and only if $(\eta_n+r\nu_{n-1})_n$ is bounded in $H$.
\end{proposition}
```

Proof: Theorem A for one direction; Gabay / Eckstein–Bertsekas (Y1)–(Y3) and E–B Thm 7 for the other. The authors should write out the five-line identification of variables.

(b) **Bùi–Combettes.** In BC, A = N_V with V = {e₀}^⊥ and B = (P_V P_K P_V)^{−1} − Id. Since P_V P_K P_V = ∇[(½|·|² − ½d_K²)∘P_V] is a firmly nonexpansive gradient of a convex function, B = ∂ψ for a convex ψ. So the example lies inside convex minimisation.

Through the correspondence in (a), it is an instance of (A1)–(A3) with V = ℝ, B: t ↦ te₀ (an isometry onto the closed range span{e₀} = V^⊥), ℓ = 0, F = ψ* − ψ*(0) ≥ 0 and r = 1. In that instance the BC "shadow" sequence x_n is the **multiplier** η_n, which converges weakly but not strongly. Meanwhile ν_n = (z_n − η_n)/r → 0 and u_n ∈ ℝ converge strongly.

So BC supports exactly this statement: *the weak convergence of η_n in Theorem A cannot be improved under (A1)–(A3)*. It says nothing about ν_n or u_n. Remark 2.12 should say this precisely instead of "only weak convergence of (ν_n, η_n) and of u_n … is the general situation".

### m1 (MINOR): Proposition 2.14 [prop:glowinski-pless2], lines 412–420.

"uniformly convex … with modulus δ(·)" is not Glowinski's definition, which is (4.6), a monotonicity bound on F₀′. Add the derivation. From (4.6),
F₀(q) − F₀(m) − (F₀′(m), q−m) = ∫₀¹ (F₀′(m + s(q−m)) − F₀′(m), q−m) ds ≥ ∫₀¹ δ₁(s|q−m|)/s ds.
Applied with m = q_M/2 to the points 0 and q_M, this gives ½(F₀(0) + F₀(q_M)) − F₀(q_M/2) ≥ ∫₀¹ δ₁(s/2)/s ds > 0. All three points lie in the unit ball, and F₁ is finite there because F is finite on L² for p < 2. The rest of the proof (and the proof that (iii) fails) is correct.

### m2 (MINOR): Remark 2.13 / line 330, GPT16 citations.

GPT16 §3.3 (pp. 22–23) identifies ALG2 (with ρ = r) with Douglas–Rachford in the **primal** special case V = H, B = I with F, G differentiable. It does not show "DR applied to the dual". Cite Gabay 1983 / Eckstein–Bertsekas §5 for the dual statement. Also, GPT16's sentence "convergence of the multiplier sequence is no better than weak in general" (p. 19–20) is stated for ALG1. For ALG2 use Glowinski 1984, §7.

### m3 (MINOR): Remark 2.11 [rem:init].

The statement is correct, but there is a cleaner identity. The three nonnegative terms combine:
(r/2)‖s_n‖² + (r/2)‖ν_n − ν_{n−1}‖² + r(s_n, ν_n − ν_{n−1}) = (r/2)‖Bu_n − ν_{n−1}‖².
So
E_n − E_{n−1} = −m_n − (r/2)‖Bu_n − ν_{n−1}‖² for every n ≥ 1, whatever η₁ is.
Monotonicity of E_n then never needs η₁ ∈ ∂F(ν₀); the cross-term sign (n ≥ 2) is used only to split the summability of ‖s_n‖² and ‖ν_n − ν_{n−1}‖².

### m4 (MINOR): Literature statement 3.4(b), lines 487–490, and Table 1 row 3.

**Confirmed** (Cianchi–Maz'ya 2018, Thm 2.3 and Rem. 2.8, §4):
- the theorem holds for **n ≥ 2**;
- the hypotheses are a ∈ C¹(0,∞), −1 < i_a ≤ s_a < ∞, Ω bounded convex and f ∈ L²;
- u is the "generalized" (approximable) solution, which is unique and "agrees with the weak solution whenever f is summable enough … for a weak solution to exist". For d = 2 and i_a > −1 one has W^{1,G}₀ ⊂ W^{1,2+i_a}₀ ↪ L², so f ∈ L² is in the dual, and the generalized solution is the weak/energy solution.

**Fix.** Drop "conditionally" and "(to be checked)". Add the hypothesis a ∈ C¹(0,∞), i.e. G ∈ C²(0,∞). Optionally note that the same argument gives (S) in d ≥ 3 for t^p/p when p ≥ 2 − 2/d, since |∇u|^{p−1} = |η| ∈ L^{2*}. This makes Theorem B(b) unconditional on convex planar domains: say so after Theorem B.

### m5 (MINOR): Literature statement 3.4(a), Table 1 row 2.

Cianchi–Maz'ya 2018, Rem. 2.6–2.7 states that the 2011/2014 papers give global gradient boundedness "for n ≥ 3 provided that f ∈ L^{n,1}", with ∂Ω ∈ W²L^{n−1,1}, and also on convex domains ((2.17) "continues to hold if Ω is any bounded convex domain"). So n ≥ 3 is **confirmed** secondhand; replace "possibly only for d ≥ 3" by "for d ≥ 3".

### m6 (MINOR): Literature statement 3.6 (line 498), "Energies not covered" (line 502), Table 1.

- **De Filippis–Mingione Thm 1.3:** (1.19) reads verbatim "q/p ≤ 1 + 1/n − 1/d, if n ≥ 2 and, when n = 2, also q/p < p": **confirmed**. Add the condition on f: f ∈ L(n,1) for n > 2, and f ∈ L²(log L)^a with a > 2 for n = 2 ((1.9)).
- **DFM Thms 1.4–1.5** concern exp(c₂(x)|Dw|^p) with **p > 1** ((1.12)) and iterated exponentials with p₀ ≥ p_m > 1 ((1.20)). So e^t − 1 − t and e^t − 1 (exponent 1) are **not covered** as stated. Replace "we have not checked" by this.
- **Hästö–Ok Thm 1.1:** C^{1,α}_loc requires (VA1) **and** ω(r) ≤ cr^β, plus (A0), (Inc)_{p−1} and (Dec)_{q−1}. The last one forces Δ₂, so exponentials are excluded. The theorem concerns minimisers of ∫φ(x,|∇u|) **without a right-hand side**.
- **Colombo–Mingione 2015a/b and BCM 2018** treat functionals with f ≡ 0 as well. Say "for f ≡ 0", since the paper's problem has f ≠ 0.

### m7 (MINOR): §3.3 and Table 1, BDS.

Thms 30 and 32 are **confirmed** in arXiv v2 numbering (Ω = (−1,1)^d, d ≥ 2; for Thm 32, p > 1 and q > p + α max{1, (p−1)/(d−1)}). The gap is shown for G = F + S° with S° a *linear continuous functional on W^{1,φ(·)}*, with zero boundary values, not for an L² load. Add "with a suitable linear term". This does not affect the logic, because a gap rules out Lipschitz minimisers.

### m8 (MINOR): §5, line 603, BOS and Carstensen–Ortner.

- **Carstensen–Ortner Thm 2.1** assumes W continuous on Ω̄ × ℝ^m × ℝ^{m×n} (so continuous in x), homogeneous Dirichlet data, exact energies on a polyhedral, exactly triangulated domain, and shape regularity. The conclusion is lim E_ℓ = inf E(A^∞), with A^∞ the Lipschitz class. State these hypotheses.
- **BOS** says conforming minima converge "to the so-called H-minimum", with **H = W ∩ W^{1,∞}** (Lipschitz, not "the closure of smooth functions"), citing CO10 for exact quadrature.
- **BOS Thm 1** (CR converges to min_W) requires the specialised one-point quadrature with (A1)–(A3) or (B), and the polynomial growth (1) with 1 < p₋ < p₊ < ∞. So exponential G are not covered.
- BOS §4.2 shows that with the "minimum" quadrature even the *conforming* scheme converges to the W-minimiser. This strengthens the paper's remark on barycentre sampling: in the gap regime, the outcome of the paper's scheme can depend on the quadrature. Say so.

### m9 (MINOR): Proposition 5.3(i), line 631. Attribution.

The interpolation argument (|∇I_hu| ≤ C_σL, the mean-value estimate along edges, dominated convergence) is the argument of CO10 (Lemma 2.1, proof of Thm 2.1), adapted to inscribed polygons. Cite it. The adaptation is correct (§2.12).

### m10 (MINOR): line 568.

"for G ≡ 0 … J_h is unbounded below": this needs ℓ_h ≠ 0.

### m11 (MINOR): line 779 and §6.6.

The ADMM dual residual is r∇*(ν_n − ν_{n−1}), a functional on V_h. ‖ν_n − ν_{n−1}‖ only bounds it. Fine as a recommendation, but say "an upper bound for the dual residual".

### m12 (MINOR): Del Pezzo–Martínez, Remark 4.2 (line 591).

Confirmed: (H4) fails for B = R_h + ∇, and their Remark 4.4 says "(H4) is only used to conclude uₙ → u". Two additions are worth making:
- their (H2) quotes Glowinski **without** the uniform convexity of F₀. This is harmless in their finite-dimensional setting (book Thms 4.2/5.2 need F₀ strictly convex and C¹), but it is a misquote;
- they work with p(x) ∈ [p₁, 2], i.e. sub-quadratic, which is consistent with the finite-dimensional theory.

The table of L² errors on five meshes, with no fitted order, is confirmed.

### m13 (MINOR, typos)

- Line 597: "A gap it is".
- Line 486: the parenthesis "(For [CM2011] we have only the abstract …)" can be removed after m5.

---

## 2. Things verified correct (line by line)

1. **Lemma 2.2** (l.s.c. via Fatou on an a.e.-convergent subsequence, plus Mazur). Correct.
2. **Proposition 2.3** (saddle points ⇔ (2.1)). Both directions are correct. The q-variation uses convexity correctly, and r-independence holds. This matches the book's Thm 2.1.
3. **Lemma 2.5** (pointwise ∂F). Correct. The sets Ω^z_M and the countable union over z ∈ ℚ^d handle measurability properly; boundedness of Ω is needed so that q′ ∈ L², and it is assumed.
4. **Lemma 2.6** (∂φ for radial φ), **(2.4)** (∂F single-valued or empty). Correct.
5. **Lemma 2.7** (well-posedness). The coercivity bound, strict convexity, the characterisation η_n + r(Bu_n − ν) ∈ ∂F(ν), and (2.8) are correct. η₁ ∈ ∂F(ν₀) is needed only for n = 0 in (2.7); see m3.
6. **Lemma 2.8 / Corollary 2.9** (pointwise prox). All four parts are correct, including the case ŵ·z < 0 in (iv).
7. **Theorem A.**
   - Step 1: the identity is checked; both polarisation identities are right.
   - Step 2: the signs follow from monotonicity of ∂F and η_n ∈ ∂F(ν_{n−1}).
   - Step 3: correct. The limit of (η_{n+1}, ν_n) is taken through the equation, which is the standard closedness argument for maximal monotone graphs, and weak l.s.c. finishes it.
   - Step 4: correct (Opial; |x_n − z|² = 2E_n(z); B injective).
   - Rellich holds on any bounded open set (extension by zero).
   - No hidden use of dom F = H, maximal monotonicity or growth.
8. **Lemma 2.10** (ψ bounds). (i): the identity and the bound ½k(|a−b|/2)|a−b|² are correct. (ii): the eigenvalues of DT, absolute continuity across the origin (|c_s| = |a−b||s−s₀|, p − 2 > −1), the lower bound (p−1)(|a|+|b|)^{p−2}|a−b|², and the Hölder bound 2^p|a−b|^{p−1} are all correct.
9. **Theorem B.** (a): the δ-splitting and the uniqueness argument (every minimiser carries η*) are correct. (b): the Hölder exponents 2/p and 2/(2−p), the convergence in L^s for s < 2 by interpolation, and 2(p−1) ≤ p for the η-estimate are correct.
10. **Proposition 2.14.** The failure of (iii) is correct (with ε < 2(2−p)); the failure of (iv) is correct up to m1.
11. **Lemma 2.15** (multiplier lemma, Hahn–Banach with Ψ′(q₀;·)). Correct, including the continuity bound.
12. **Proposition 2.16.**
    - (i) Correct. W^{1,p}₀ ↪ L² is used exactly where needed.
    - (ii) Correct. In the converse, on Lipschitz Ω, H¹₀ ∩ {∇v ∈ L^p} ⊂ W^{1,p}₀ (TODO #53). Argument: the zero extension ũ has compact support, ũ ∈ W^{1,1} and ∇ũ ∈ L^p, so ũ ∈ W^{1,p}(ℝ^d) by Poincaré on a ball; a W^{1,p}(ℝ^d) function supported in Ω̄ belongs to W^{1,p}₀(Ω) for Lipschitz (segment-property) Ω, by the standard translation argument in charts. A textbook reference can be added.
13. **Lemmas 3.1–3.2** (Lipschitz solutions, explicit flux). Correct.
14. **Theorem C.**
    - (i): ‖q‖_{L¹} ≥ a_min^{1/2}‖q‖_h, coercivity and the strict-convexity argument for φ are correct.
    - (ii): correct (but classical, M2).
    - (iii): correct.
    - Remark 4.1 (nonhomogeneous data): the change of variables is exact, and the initialisation is admissible because η₁ = 0 ∈ ∂F_h(0).
15. **Lemma 5.1.** Correct (W^{1,1}₀ density, finiteness of J(u) through the subgradient inequality at 0).
16. **Remark 5.2** (radial (S)). Correct: τ is measurable as the countable sup, bounded, the directions agree, the dead core is handled with Lemma 2.6(i), and the two-phase interface is a null set.
17. **Proposition 5.3.**
    - (i)–(ii): J_h = J on V_h (G autonomous; the barycentre rule is exact for affine functions; Ω_h ⊂ B₁).
    - Interpolation (TODO #51): ∇u is in fact continuous on the closed disk, since τ(0) = 0 and τ is continuous, and the edge-wise mean-value estimate is valid because u is C¹ off the origin and Lipschitz.
    - Compactness of Q, positivity of μ(ε) by strict convexity, and D_h = J(u_h) − J(u): correct.
    - (iii): ∫_{1/2}^1 gives the constant 3/16 ≥ 1/8. Correct.
18. **Claims checked for the test G:**
    - t^{1.5}/1.5 + t³/3: tg′/g ∈ [½, 2];
    - t²/2 + t⁴/4: tg′/g ∈ [1, 3];
    - t² log(1+t): tg′/g → 2 at 0 and → 1 at ∞, and k = 2 log(1+t) + t/(1+t) is nondecreasing;
    - t³/log(1+t): k is nondecreasing and g′ ≥ 2 (numerical check);
    - |∇u| ∼ |x|^{1/3} for p = 4.

    All correct.

---

## 3. Literature claims checked against the sources

| Claim (location) | Status |
|---|---|
| Glowinski Ch. VI: hypotheses (i)–(iv), ALG2 strong convergence of (u_n, ν_n), ρ < (1+√5)r/2, Thms 4.1/5.1 (§2.6) | **CONFIRMED** (pp. 171–182). But see M1: the hypotheses hold for p ≥ 2 and for all uniformly convex G |
| Glowinski: uniform convexity required (TODO 1) | **CONFIRMED**: (4.5)–(4.6). Del Pezzo–Martínez omit it |
| Glowinski "Rem. 4.4" | **NOT FOUND** in Ch. VI (Remarks 4.1–4.2 only); "Remark 4.4" is Del Pezzo–Martínez's own |
| GPT16: ν-step first; Remark 21 (ρ = r ⇒ ADMM); multipliers only weak | **CONFIRMED** (the last for ALG1) |
| GPT16 §3.3: ALG2 = DR | **CONFIRMED only in the primal case** V = H, B = I (m2) |
| E–B Thm 7 (weak convergence of z^k), Thm 8 (finite-dim, full column rank, KT pair) | **CONFIRMED** (Math. Prog. 55 (1992) 293–318) |
| Bùi–Combettes: DR converges only weakly | **CONFIRMED**. Within ADMM it concerns the multiplier only (M4b) |
| Lieberman Thm 1.7 C^{1,β}(Ω̄) | **WRONG**: interior; boundary via Lieberman 1988 (M3) |
| CM2018 Thm 2.3 for d = 2, notion of solution | **CONFIRMED** (n ≥ 2; generalized = weak) |
| CM2011: d ≥ 3 | **CONFIRMED** secondhand (CM2018 Rem. 2.6) |
| CM2015a Thm 1.2 strict q/p < 1+α/n; CM2015b q ≤ p+α for bounded minimisers; BCM Thms 1.1, 1.3, 1.4 | **CONFIRMED** (interior, f ≡ 0) |
| DFM Thm 1.3 incl. "q/p < p if n = 2" | **CONFIRMED** verbatim; add the f-condition |
| DFM Thms 1.4–1.5 vs e^t − 1 | **Not covered** (need p > 1) |
| Hästö–Ok C^{1,α}_loc under (VA1) | **Incomplete**: also needs ω(r) ≤ cr^β; f ≡ 0 |
| ELM: W^{1,q}_loc for q/p < (n+α)/n; counterexample 1 < p < n < n+α < q with a gap | **CONFIRMED** (Thms 2, 3) |
| BDS Thms 30, 32 | **CONFIRMED** (v2 numbering; gap with a linear functional S°) |
| BGS Thm 2.3 (bounded Lipschitz domain, q ≤ p + α max{1, p/d}, W^{1,q} data, C_c^∞ density) | **CONFIRMED**; sharp for p ≤ d vs BDS: consistent |
| CO10: conforming FEM fails ⇔ Lavrentiev | **CONFIRMED** with hypotheses (m8) |
| BOS: conforming → H-minimum; CR Thm 1; quadrature sensitivity | **CONFIRMED**, but "closure of smooth functions" should be "W ∩ W^{1,∞}"; CR hypotheses missing (m8) |
| Del Pezzo–Martínez: (H4) fails, used only for u_n; errors without fitted order | **CONFIRMED** |

---

## 4. TODO_VERIFY.md, item by item

Legend: C = CONFIRMED, W = WRONG (with correction), NF = NOT FOUND or not checkable from the provided PDFs.

**A. Glowinski / ADMM**
1. **C** for the hypotheses (4.1)–(4.6), the numbering (Ch. VI Thm 4.1 = ALG1, Thm 5.1 = ALG2), "strongly" ((5.2)–(5.3)), ρ < (1+√5)r/2, and uniform convexity being required. **W** for the p > 2 sub-item: (iv) holds literally via the splitting in M1, so the AFR23 application for p(x) ≥ 2 (with p₊ < ∞) is justified. Extra: whole-sequence weak convergence of λ, book §7. The s-Laplacian is excluded at the continuous level in W^{1,s} (book §6.3).
2. **C**: ν-step first, Remark 21, "no better than weak" (for ALG1). §3.3 DR identification: only in the primal case B = I (m2). Pages 19–94: **NF** (preprint pagination).
3. **C** (indirect): the golden ratio is in book Thm 5.1, and Ch. VI "follows Fortin and Glowinski [2]". Pages: **NF**.
4. **NF** (Gabay not in the PDFs). E–B §5 / Thm 8 give the ADMM–DR link and can be cited instead.
5. **C** secondhand: E–B Cor. 6.3 attributes weak convergence of DR to Lions–Mercier [31].
6. **C**: Thm 7 (weak convergence of z^k; unbounded if no zero) and Thm 8 (full column rank, Kuhn–Tucker pair, finite dimension). DOI: **NF**.
7. **C** secondhand: BC write "x_n ⇀ x … was first established in [7]" (Svaiter). Volume/pages: **NF**.
8. **C** (two maximally monotone operators, separable H). Pages: **NF**. Scope: see M4b.
9. **C** secondhand: BC cite "[1, Theorem 26.11(iii)]".
10. **NF**.
11. **NF**.
12. **NF** in the PDFs. From memory, Boyd App. A uses V = (1/ρ)‖y − y*‖² + ρ‖B(z − z*)‖², consistent with the claim; but cite Glowinski (5.25) first (M2).
13. **C**, plus the misquote of uniform convexity (m12). The "order 0.8" removal is correct: they fit no order.
14. **NF** (AFR23 not provided).

**B. Regularity**

15. **C**: g = g(t), x-dependence (1.10c) with Λ₁(1+g(|p|))[|x−y|^α + |z−w|^α], and u-dependence allowed. **W**: Thm 1.7 is interior (M3); the boundary version is via Lieberman, Nonlinear Anal. 12 (1988), Thms 1–2. Antonini's attribution is consistent with this.
16. **NF**.
17. **C** (n ≥ 3; convex domains also included per CM2018 (2.17)). DOI: **NF**.
18. **NF** (JEMS paper not provided).
19. **C**: (a) n ≥ 2; (b) generalized (approximable) solution, unique, equal to the weak solution when f is in the dual, which holds for d = 2; (c) a ∈ C¹(0,∞), −1 < i_a ≤ s_a < ∞, Ω convex bounded, f ∈ L². The conditional conclusion becomes unconditional (m4).
20. **NF**.
21. **NF**.
22. **C** (Thm 1.2 under (1.3) q/p < 1 + α/n; interior; f ≡ 0). ARMA 215 (2015) 443–496 **C** from the PDF.
23. **C** (Thm 1.1: u ∈ W^{1,p} ∩ L^∞, q ≤ p+α, Du ∈ C^{0,β}_loc). Journal: ARMA, doi 10.1007/s00205-015-0859-9 (from the PDF). Vol. 218/pages: **NF** in the text.
24. **C** (Thms 1.1, 1.3, 1.4 as stated; Thm 1.4 is local, on balls B ⋐ Ω). Journal data: **NF** (arXiv).
25. **C**: "q/p < p when n = 2" is in the original. Thms 1.4–1.5 do **not** cover e^t − 1 − t or e^t − 1 (p > 1 required) (m6). Volume/pages: **NF**.
26. **C** interior, but C^{1,α} also needs ω(r) ≤ cr^β, and f ≡ 0 (m6).
27. **NF** (no Missaoui PDF). Nothing in the provided sources contradicts the claim.

**C. Counterexamples / Lavrentiev**

28. **C** (Thms 2, 3 of ELM; scalar, isotropic, Lipschitz boundary datum).
29. **NF**.
30. **NF**. Secondhand: DFM attribute the double phase energy to Zhikov [77, 78]. The exact log-Hölder no-gap reference is not found.
31. **C** (v2 numbering; note the linear term S°, m7). Published numbering: **NF**.
32. **C** (Thm 2.3; bounded Lipschitz Ω; the v1 typo is confirmed — §3 has the correct q ≤ p + α max{1, p/d}).
33. **NF**.
34. **NF**.
35. **NF**.
36. **C**, with the hypotheses of m8.
37. **C** with corrections: H = W ∩ W^{1,∞}; Thm 1 hypotheses; quadrature sensitivity (§4.2) (m8). Volume/pages: **NF**.
38. **C** secondhand (BOS: [Ort11; OP11] for autonomous convex integrands).
39. **NF**.

**D. FE rates**

40–45. **NF** (none provided).
46. **NF**. The book §6.3 refers to Glowinski–Marrocco [5] for "a detailed analysis, including error estimates", which is consistent.
47. **NF**. Not contradicted by any source read; I also know of no general O(h²) L² proof for the p-Laplacian. Keep "to our knowledge".

**E. Code / data**

48–50, 52, 54: outside the mathematical scope (not checked).
51. **C** (Prop. 5.3 proof, including the interpolation step; see §2.17). Remark 5.2 is rigorous: keep it.
53. **C** (argument in §2.12).

**F. Bibliographic data checkable from the PDFs**

- Glowinski Ch. VI pp. 166–194 **C**.
- E–B 55 (1992) 293–318 **C**.
- CM2015a 215 (2015) 443–496 **C**.
- ELM 204 (2004) 5–55 **C**.
- CO10 CMAM 10 (2010) no. 2, 137–163 **C**.
- Lieberman CPDE 16 (1991) 311–361: consistent (pp. 314–320 seen).
- All other items in F: **NF**.
