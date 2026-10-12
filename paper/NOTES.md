# NOTES: editorial and mathematical decisions (main.tex)

## Structure (as requested in BRIEF.md)
1. §2 Continuous level under (S): abstract structure (A1)–(A3); Prop. 2.3 (saddle points); Lemma 2.4 (pointwise ∂F); Lemma 2.5 (∂Φ); algorithm; Lemma 2.6 (well-posedness); Lemma 2.7 (pointwise problem g(x,t)+rt=r|w|: existence, uniqueness, bracket [0,|w|], dead-core threshold); **Theorem A** (with the displayed Lyapunov identity (2.9)); Lemma 2.12; **Theorem B** (a)/(b); §2.6 relation with Glowinski's theorem, with Prop. 2.13 (why p ≥ 2); §2.7 Lemma 2.14 (multiplier lemma) and Prop. 2.15 (what (S) means for t^p/p).
2. §3 When (S) holds or fails: Lemmas 3.1/3.2; "Literature statements" 3.3–3.6; energies not covered; counterexamples; Table 1.
3. §4 Discrete problem, **Theorem C** with full proof (coercivity, then the multiplier via the Hahn–Banach lemma 2.14, then Theorem A verbatim); Remark 4.1 (nonhomogeneous data, discrete); Remark 4.2 (nonhomogeneous data, continuous); Remark 4.3 (literature, no rates).
4. §5 Lavrentiev: known gap and no-gap results; consequences for conforming P1; Lemma 5.1; Remark 5.2 (the test problems); Proposition 5.3 (convergence in energy of P1 minimisers, autonomous radial case).
5. §6 Numerics; §7 conclusions and open problems; software availability; acknowledgements TODO.

## Mathematical decisions
- **No Rockafellar citation.** Rockafellar's theorems (maximal monotonicity; subdifferentials of integral functionals) are not in the REFERENCES files. I avoided them:
  - Step 3 of Theorem A uses the subgradient inequality plus weak lower semicontinuity, not maximal monotonicity.
  - The pointwise characterisation of ∂F is proved directly (Lemma 2.4).
  - The prox optimality condition is proved directly (Lemma 2.6(ii)).
- **Multiplier lemma (Lemma 2.14).** This is a short Hahn–Banach argument. It gives (S) for the discrete problem (Thm C(ii)), and for 1<p<2 it gives "(S) ⟺ J attains its minimum on H¹₀ ⟺ u_p ∈ H¹₀" (Prop. 2.15), with complete proofs and no external convex-analysis theorem.
- **Assumption (A2).** F ≥ 0, convex, lsc and proper (weakened after refereeing from F(0)=0, so that the shifted functional for nonhomogeneous data is covered). F ≥ 0 avoids needing affine minorants.
- **REVIEW corrections incorporated.**
  - E1: the multiplier step of Prop. 2.3 (⇒).
  - E2: what (S) means for p<2 vs p>2.
  - E3: superlinearity in Thm C.
  - E4: no claim that strong convergence fails for (P).
  - G1: ∂F as an inclusion / pointwise.
  - G2: uniform-in-x bound in Lemma 3.1.
  - G3: initialisation (Remark 2.9).
  - G4: nonsmooth prox (Lemma 2.7).
  - G5: Thm B(b) strengthened. η_n → η̄ strongly and L^s convergence for s<2, using a self-proved Hölder bound |T(a)−T(b)| ≤ 2^p|a−b|^{p−1} (constant 2^p, weaker than 2^{2−p}, but with an elementary proof). Uniqueness of the minimiser is added in Thm B(a).
- **Remark 5.2 deliberately goes slightly beyond Correction 6 of the brief.** For the *explicit radial test data* (disk, f = c), the flux is η = −cx/2. Lemma 3.2 then gives (S) rigorously, also for e^t−1−t, e^t−1 and the indicator two-phase coefficient, and Lemma 5.1 rules out any Lavrentiev obstruction.
  - So Theorems A/B do apply to *those particular problems*; the column "theory" in Table 3 says so.
  - The text stresses that this is a statement about these explicit data only. For general data no (S) is known for e^t−1−t and e^t−1, and no claim is made (§3.2, "Energies not covered").
  - **If the authors prefer to follow the brief literally, delete Remark 5.2 and the "theory" column of Table 3.**
- **Cianchi–Maz'ya ARMA 2018, Thm 2.3.** This is used to argue (conditionally) that (S) holds on convex plane domains with f ∈ L² for every autonomous G of that class, including t^p/p with all p>1. This is the route suggested in REVIEW E2(iii). It is marked conditional on the notion of solution.
- **Ordering of ALG2.** Glowinski–Pan–Tai write the ν-step before the u-step. Our order (u, then ν) is the one of Aragón–Fernández Bonder–Rubio and of the code. This is flagged in TODO_VERIFY.

## Numerical decisions
- **Main torsion study.** All numbers in Tables 3–4 and Fig. 1 come from `materials/convergence_order.csv`. Orders were refitted by me from the CSV (3 finest meshes) and agree with the README to 2 decimals. Fig. 1 was regenerated from the CSV (`scripts/make_figures.py`; two-phase legend now reads a = 1_{|x|>1/2}) with the 4 cases with reduced sup-order highlighted. The other 6 cases are drawn in gray, because 10 colour-coded series are unreadable.
- **New runs (by me, glaplace from /home/claude/glaplace, unchanged).** Raw output is in `results/`.
  - `scripts/lyapunov_identity.py` gives `results/lyapunov.json` and `lyapunov_curves.npz`, used for Table 8 and Fig. 3. The saddle proxy is η* = η_{L+1} + r(ν_L − ν_{L−1}), which satisfies the discrete multiplier equation exactly (verified, relative residual ≤ 4e−13).
  - `materials/lyapunov_check.py` was also run (`results/lyapunov_check_original.txt`). It gives consistent conclusions: V nonincreasing, cross terms and m_n ≥ 0 up to 1e−26.
  - `scripts/extra_experiments.py sc` gives `results/sc.txt` (supercloseness, Table 2; tolerance 1e−10).
  - `scripts/extra_experiments.py pto1` gives `results/pto1.txt` (p = 1.2, 1.1, Table 5). **All 8 runs hit max_iter = 6000 without reaching 1e−9.**
  - `scripts/extra_experiments.py rman` gives `results/rman.txt`. This re-runs the manufactured order study and the r study; the r-iteration counts reproduce the README exactly.
  - `scripts/extra_experiments.py rdead` gives `results/rdead.txt`. This re-runs the dead-core r study and reproduces the README (r = 1: 6000*, 3: 2139, 10: 650, 30: 244). It also adds r = 100 (529 iterations), whose error differs in the 3rd digit; this is reported as a caveat on the stopping criterion.
  - `materials/examples/two_phase.py` (no --plot) gives `results/two_phase.txt` (§6.8 numbers).
- **Supercloseness.** True versus nodal: the factor is 7–25, and the nodal "order" is about 2.5.
- **"7 p≈1 runs" mentioned in the brief.** This refers to a different data set: `paper_v2/results/convergence_order_true.csv`, from the pxlaplace code, stopped with rel_u. I did not use it. Instead I re-ran p = 1.2 and 1.1 with glaplace (residual stop). In that study 8 of 8 runs did not converge. In the paper the limitation is now stated as "twelve runs did not reach the tolerance" (2 expdead + 8 p→1 + the two r = 1 runs of Table 7).
- **Timings.** These are not reported (shared environment).
- **Exact radial solutions.** These were checked against the closed forms for p = 1.5, 2, 4 (max relative deviation 9e−9 for p = 4).
- **two_phase.png and exp_dead_core.png.** These were copied from `paper_v3/figures/`. They were produced by the repo scripts and not regenerated by me; the parameters are stated in the captions.

## Compilation
- `pdflatex main.tex` was run twice (three times after cleaning): 29 pages, no errors, no overfull boxes, no undefined references.
- The bibliography is an inline `thebibliography`.
- The T1 fontenc and microtype font expansion are disabled because the TeX installation lacks scalable EC fonts.

## Changes after refereeing (REFEREE_REPORT.md)

All must-fix errors E1–E11 and the gaps G1–G17 were addressed; minor points 1–11 were addressed (12 = TODO placeholders remain). The pre-referee version is kept as `main_before_referee.tex`.

**Must-fix errors**
- **E1 (nonhomogeneous data).**
  - (A2) is weakened to "F convex, lsc, proper, F ≥ 0".
  - Lemma 2.6(ii) now uses dom F ≠ ∅. The initialisation sentence in §2.3 is revised (ν₀ = η₁ = 0 is admissible when F(0) = 0; by Rem. 2.9 the condition is not needed).
  - Remark 4.1 now explains the shifted discrete functional F̃_h (F̃_h(0) ≠ 0) and that the code's initialisation corresponds to q₀ = −∇ũ_h, η₁ = 0 ∈ ∂F_h(0).
  - New Remark 4.2 gives the continuous lift u = u_D + v, F̃(q) = F(q+∇u_D), so Prop. 2.3 and Thms A–B apply and Lemma 3.2 holds with u − u_D ∈ H¹₀.
  - §6.4 now cites Remarks 4.1/4.2.
- **E2 (point values).** §4.1 now assumes G(x,·) is a Young function for every x and f has values at the barycentres (or a fixed representative / exact integration). Theorem C refers to this.
- **E3 (superlinearity).** For t^p log^a(1+t) the condition is now "p > 1, or p = 1 and a > 0".
- **E4 (Lieberman ranges).** These are now (1,2] for t² log(1+t) and [1,2) for t³/log(1+t), with the remark that tg'/g ≥ 1 (g/t nondecreasing).
- **E5.** The false "isotropic ⇒ uniformly elliptic" sentence is replaced by the referee's wording. Marcellini 1991 §6 is cited specifically for the counterexample.
- **E6.** The double-phase (S) condition now reads |∇u|^{p−1} + a|∇u|^{q−1} ∈ L².
- **E7.** In Table 1, the Cianchi–Maz'ya 2011 row is downgraded to "yes, conditionally (range of d unverified; d ≥ 3 for JEMS); not established for d = 2 from the sources read". Fact 3.4(a) adds "(possibly only for d ≥ 3)".
- **E8.** The abstract and Contribution (2) now state Theorem B(a) with all of its hypotheses: (S), C¹, g(0) = 0, g/t nondecreasing and bounded below by a positive function of t uniformly in x.
- **E9.** The Breit–Diening–Schwarzacher attribution for discontinuous coefficients is deleted ("not covered by any estimate we know"). The "≈ 0.8" for Del Pezzo–Martínez is replaced by "report L² errors on a sequence of meshes without fitting an order".
- **E10.** Twelve non-converged runs are now counted. The expdead row of Table 3 is marked † with a caption note. The p→1 fits are said to use only unconverged runs, and the expdead observation notes the same.
- **E11.** The formula is now u(x) = ∫_{|x|}^1 …, and the figure legend reads 1_{|x|>1/2}.

**Gaps and overstatements**
- **G1 (convergence of u_h → u).**
  - New Proposition 5.3, with full proof, for the autonomous radial problems: J(u_h) → J(u), the Bregman remainder D_h → 0, ∇u_h → ∇u in measure, and u_h → u in H¹₀ when g(0) = 0 and g/t is nondecreasing and positive.
  - The proof uses J_h = J on V_h (exact barycentre rule for constant f, autonomous G), Lemma 5.1 for the lower bound, and interpolation (shape regularity, u ∈ C¹(B̄₁∖{0})) for the upper bound.
  - Proposition 5.3 does not cover the indicator two-phase problem; the text says so. The text states exactly what Lemma 5.1 gives: m_W = m_Lip = J(u) and J(v_h) ≥ J(u).
  - The gap is now defined with smooth functions, with m_Lip introduced separately. "What is not proved" gains items (iv) (u_h → u) and (v) (d = 2). Contribution (5) and the conclusions are revised to match.
  - The integrability note Φ(x,∇u) ≤ η·∇u ∈ L^∞ is added to the proof of Lemma 5.1.
- **G2 (why p ≥ 2).**
  - The introduction and §2.6 now claim only that the hypotheses as reproduced in GPT16 fail for p < 2 (Prop. 2.13).
  - New paragraph: for p > 2, F = +∞ off L^p, so F is not Gâteaux differentiable on H. F − (c/2)‖·‖² is not convex near 0, so there is no obvious splitting as in (iv). Hypothesis (iv) must therefore be read in a weaker sense for the theorem to apply as in AFR23. This is flagged in TODO_VERIFY item 1.
  - The subsection is retitled "the restriction p ≥ 2".
- **G3.** The word "strongly" is now marked "as we read [GPT16]".
- **G4.** Boyd et al. App. A is hedged ("as far as we recall").
- **G5.** Remark 5.2 adds the hypothesis inf_s g₀(s,t) → ∞, and proofs that τ is measurable (as a monotone inverse) and bounded.
- **G6.** Remark 2.11 is softened to "the general theory gives only weak convergence".
- **G7 (p → 1 runs).**
  - New run (`extra_experiments.py pto1_hist`, `results/pto1_hist.txt`) recording the L² error at iterations 1500/3000/4500/6000.
  - For p = 1.2 the error is stable to about 3×10⁻⁵ (relative). For p = 1.1 it still drifts by about 3% between iterations 3000 and 6000; this is reported.
  - "Errors are nevertheless meaningful" and "consistent with Thm B(b)" are removed.
- **G8.** The text now recommends monitoring the dual residual as well; no rerun was done.
- **G9 (Lyapunov check).**
  - `lyapunov_identity.py` now also reports statistics over n ≤ n_a, the iterations with E_{n−1}/E₀ > 10⁻²⁰.
  - Table 8 is rebuilt on that range. There E_n is strictly decreasing and both sign terms are strictly positive (min 8.6e−25 … 4.7e−21).
  - The text now says that the identity check tests the implementation, not the theorem; that the full-range values are round-off (≤ 3e−28); and that the decay rate near n ≈ L is affected by the proxy.
- **G10.** Iteration growth is now described as sublinear (≈ h^{−1/2} on the finest meshes for p = 1.5), and the list of growing cases adds t² log(1+t) and e^t−1.
- **G11.** e^t−1 (sup order 1.11) is added to observation (2).
- **G12.** The heuristic in §6.8 now excludes the origin (for p ≠ 2) and the dead core.
- **G13.** "What is not proved" is extended with items (iv) and (v).
- **G14.** The Table 3 caption says the theory column concerns the continuous iteration and that Theorem C covers the discrete iteration in all cases.
- **G15.** The ess-inf sentence in the proof of Thm B(a) is removed; the proof uses the hypothesis at t = δ/2 only.
- **G16.** The condition "q/p < p if d = 2" is added to the De Filippis–Mingione statement.
- **G17.** The only DOI (Antonini) is removed, so no entry has a DOI. The MEM bibliographic data are listed individually in TODO_VERIFY section F.

**Minor points**
- **Notation.**
  - The dimension is d throughout §3 (the integrability exponent is renamed σ); the iteration index stays n.
  - g/t is renamed k(x,t).
  - Triangles are named K (the stiffness matrix is A; the coercivity constant is M).
  - The Diening–Kreuzer function is renamed 𝒱.
  - (q,q')_h replaces (p,q)_h. ψ in Lemma 2.7 is renamed χ. Orders are written κ.
  - The double-phase energy is written explicitly, not as H(x,t).
- **Supercloseness wording** is corrected ("I_h u is superclose to u_h"; spurious orders 2.1–2.5).
- **Theorem A** last line now reads "ū is the unique minimiser".
- **Prop. 2.15(ii)** gains a converse on Lipschitz domains, with proof.
- **Table 7** caption: 4.413–4.414×10⁻⁵.
- **Actual mesh sizes** are added: h_max = 0.288, 0.136, 0.0625, 0.0313, 0.0161 (`results/meshinfo.txt`).
- **Exact-solution check** is now backed by `results/exactcheck.txt`.
- **Section titles** with math use `\texorpdfstring`; there are no hyperref warnings.

**Compilation.** pdflatex, 3 passes: 29 pages, no errors, no overfull boxes, no undefined references. Pages were checked visually.

**Remains open.** Everything in TODO_VERIFY. In particular:
- the Glowinski hypotheses for p > 2 (item 1);
- Cianchi–Maz'ya in d = 2;
- the MEM bibliographic data;
- whether to keep Remark 5.2 / Prop. 5.3.
