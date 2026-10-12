# TODO_VERIFY — items the authors must check before submission

Labels as in REFERENCES_A/B:
- FT = read in full text via WebFetch (model-mediated)
- ABS = abstract only
- SEC = known from citations in other papers
- MEM = from memory

Items 1–52 date from the first draft (some are updated after refereeing); items 53–54 and section F were added after refereeing. Each item gives the exact claim in main.tex, where it is made, and the source.

## A. Glowinski's theorem and ADMM theory

1. **Glowinski 1984, Ch. VI, Thms 4.1/5.1 (ABS: chapter metadata only).**
   - Claim (§2.6): the hypotheses are (i)–(iv) as listed: saddle point; B injective with closed range; F(q)/‖q‖ → ∞; F = F0 + F1 with F0 Gâteaux-differentiable and uniformly convex on bounded sets.
   - Claim (§2.6): ALG2 with 0 < ρ < (1+√5)r/2 gives (u_n, ν_n) → (u, Bu) strongly, with the multiplier converging only weakly.
   - Claim (§2.6): the numbering is Ch. VI, Thms 4.1 and 5.1.
   - Source: Glowinski–Pan–Tai 2016 §3.2 (FT-WF); the numbering is inferred from Del Pezzo–Martínez and Aragón–FB–Rubio.
   - **Check against the book.** In particular check whether uniform convexity is required: Del Pezzo–Martínez omit it.
   - **NEW (after refereeing): the case p > 2.** For F(q) = ∫|q|^p/p, F = +∞ on L² \ L^p, so F is not Gâteaux differentiable on H. Also, F − (c/2)‖·‖² is not convex near 0, so (iv) as reproduced by GPT16 does not hold literally. §2.6 now says that (iv) must be read in a weaker sense (e.g. F₀ differentiable on dom F) for the theorem to apply as in Aragón–FB–Rubio Thm 4.3. **Check the exact hypothesis in the book, and whether the application for p(x) ≥ 2 in [AFR23] is justified.**
   - "Strongly" in the conclusion is our reading of GPT16 (REFERENCES_B: "by context"). Check it.
2. **Glowinski–Pan–Tai 2016 (FT-WF).**
   - Claim (§2.6): ALG2 is written there with the ν-step before the u-step.
   - Claim (§2.6): "multipliers converge in general only weakly".
   - Claim (§1): Remark 21 identifies ρ = r with ADMM.
   - Claim (Rem. 2.10): §3.3 identifies ALG2 with Douglas–Rachford.
   - Page range 19–94 is MEM.
3. **Fortin–Glowinski 1983, Ch. III (MEM).**
   - Claim (§2.6): the golden-ratio bound ρ < (1+√5)r/2 is due to them; this is cited together with GPT16.
   - Pages are not verified.
4. **Gabay 1983, Ch. IX (MEM).**
   - Claim (Rem. 2.10): ADMM is Douglas–Rachford applied to the dual.
5. **Lions–Mercier 1979 (secondhand).**
   - Claim (Rem. 2.10): weak convergence of Douglas–Rachford; cited generically.
6. **Eckstein–Bertsekas 1992.**
   - Thm 7: weak convergence of Douglas–Rachford (Rem. 2.10).
   - Thm 8: finite-dimensional ADMM, requiring a Kuhn–Tucker pair and M of full column rank (Rem. 4.2).
   - Status FT-WF with OCR noise. The DOI was not included.
7. **Svaiter 2011, Thm 1 (FT-WF).**
   - Claim: weak convergence of Douglas–Rachford. Volume/pages are MEM.
8. **Bùi–Combettes 2020 (FT-WF).**
   - Claim (Rem. 2.11): an example in which Douglas–Rachford converges weakly but not strongly (two maximally monotone operators in a separable Hilbert space).
   - Pages are MEM.
9. **Bauschke–Combettes 2017.**
   - Claim (Rem. 2.10): cited generically for weak convergence of Douglas–Rachford (Thm 26.11 in the REFERENCES; MEM).
   - No theorem number is used in the text.
10. **Moradifam–Nachman arXiv:1112.1960 (ABS).**
    - Claim (Rem. 2.10): weak convergence of ADMM/split Bregman in infinite-dimensional Hilbert spaces via Douglas–Rachford.
11. **Bartels–Milicevic arXiv:1704.06069 (FT-WF).**
    - Claim (Rem. 2.10): the Lyapunov quantity D_j² of their Thm 3.7 equals 2rE_n for τ = r.
    - Claim (Rem. 4.2): h-dependent constants are observed in their numerics.
    - The journal version has not been identified.
12. **Boyd et al. 2011.** (Rem. 2.10 is now hedged "as far as we recall".) §6.6 also cites Boyd et al. for "monitoring primal and dual residuals is standard practice" (§3.3; MEM).
    - §3.2 (FT-WF), cited in Rem. 4.2: residual and dual convergence in finite dimension.
    - App. A (MEM), cited in Rem. 2.10: the same Lyapunov functional. **Check App. A.**
13. **Del Pezzo–Martínez arXiv:1401.2442 (FT-WF).**
    - Claim (Rem. 4.2): convergence of ALG1/ALG2 for DG via Glowinski's theorems in finite dimension.
    - Claim (Rem. 4.2): the injectivity hypothesis (H4) fails for their B and is used only for u_n → u.
    - The claimed order of about 0.8 has been REMOVED after refereeing. The text now says only that they report L² errors on a sequence of meshes without fitting an order. Check that this is accurate.
14. **Aragón–Fernández Bonder–Rubio 2023 (FT-WF).**
    - Claim (§2.6): Thm 4.3 applies Glowinski's Thm 5.1 for p(x) ≥ 2, giving strong H¹₀ convergence.
    - Claim (§6.10): the numerics are on a single mesh.

## B. Regularity, Assumption (S)

15. **Lieberman 1991 (structure condition FT from the ANU proceedings; Thm 1.7 SEC/MEM).** Claims in Literature statement 3.3 and Table 1:
    - g depends on t only;
    - the x-dependence is only a Hölder perturbation Λ(1+g(|ξ|))|x−y|^θ (MEM);
    - Thm 1.7 gives C^{1,β}(Ω̄) for bounded solutions, ∂Ω ∈ C^{1,α}, C^{1,α} data (SEC/MEM);
    - Antonini attributes boundary regularity to another Lieberman paper.
    - **Check the boundary part of Thm 1.7 in the Comm. PDE text.**
16. **Antonini 2026, Cor. 1.4 (FT; arXiv v3).**
    - Claim: bounded C^{1,α} domain, f ∈ L^q with q > n, C^{1,α} Dirichlet data, −1 < i_a ≤ s_a < ∞, x-Hölder dependence with a single b, gives u ∈ C^{1,β}(Ω̄).
    - In the REFERENCES the integrability is written "f ∈ L^d, d > n"; main.tex writes it as L^q, q > n. Check the exponent.
    - Check the journal data (doi 10.1016/j.na.2026.114266).
17. **Cianchi–Maz'ya, Comm. PDE 2011 (ABS + MEM).** In Table 1 this is now "yes, conditionally; not established for d = 2 from the sources read".
    - Claim: global Lipschitz bound on convex domains or with ∂Ω ∈ W²L^{n−1,1}, f ∈ L^{n,1}.
    - **The range of n (the REFERENCES say n ≥ 3, MEM) and the DOI must be checked.**
18. **Cianchi–Maz'ya, JEMS 2014, Thm 3.5 (FT).**
    - Claim: p-type a with 2 ≤ p < n, n ≥ 3; the paper shows that boundedness of the gradient can fail without the domain assumptions.
    - Journal data MEM.
19. **Cianchi–Maz'ya, ARMA 2018, Thm 2.3 (FT).**
    - Claim (after Literature statement 3.4): on convex Ω with f ∈ L², a(|∇u|)∇u ∈ W^{1,2}(Ω).
    - **Check:** (a) that the theorem holds in dimension n = 2; (b) the notion of solution (weak solution in W^{1,G}_0, or generalized/approximable); (c) the hypotheses on a.
    - The conditional conclusion "(S) holds on convex plane domains for every autonomous G with −1 < i_a ≤ s_a < ∞, including all p > 1" depends on this.
20. **Fan 2007 (SEC).**
    - Claim (Literature statement 3.5, Table 1): p ∈ C^{0,μ}(Ω̄), ∂Ω ∈ C^{1,β}, bounded weak solutions are C^{1,α}(Ω̄).
    - The hypotheses were not read in the original; the DOI is MEM.
21. **Coscia–Mingione 1999, Acerbi–Mingione 2001 (MEM).**
    - Claim: interior C^{1,β} for Hölder exponents; log-Hölder exponents give Hölder continuity of u.
    - Acerbi–Mingione 2002 (ARMA 164, MEM) is cited only as a motivation for electrorheological fluids.
22. **Colombo–Mingione 2015a, Thm 1.2 (FT preprint).**
    - Claim: q/p < 1+α/n (strict) gives Du ∈ C^{0,β}_loc.
23. **Colombo–Mingione 2015b (ABS).**
    - Claim: bounded minimisers with q ≤ p+α are C^{1,β} locally.
24. **Baroni–Colombo–Mingione 2018 (FT).**
    - Thm 1.1: borderline q/p ≤ 1+α/n.
    - Thm 1.3: u ∈ C^{0,γ} and q < p + α/(1−γ).
    - Thm 1.4: no gap, local.
    - Journal data MEM.
25. **De Filippis–Mingione 2021 (FT).** The text now includes the extra condition "q/p < p if d = 2" exactly as recorded in REFERENCES_A. It looks odd (perhaps a garbled transcription); check it.
    - Thm 1.3: a ∈ W^{1,d}, d > n, q/p ≤ 1+1/n−1/d gives local Lipschitz bounds.
    - Thms 1.4–1.5: exponential integrands.
    - **Whether Thms 1.4–1.5 cover e^t−1−t or e^t−1 was NOT checked** (stated as such in §3.2).
    - Volume/pages MEM.
26. **Hästö–Ok 2022, JEMS (FT, main theorem).**
    - Claim: C^{1,α}_loc under (VA1); interior only.
    - Volume/pages MEM.
27. **Claim (§3.2, Table 1): no refereed global (up to the boundary) gradient bound is known for x-dependent double phase.**
    - This is based on the librarian's search.
    - **Not cited in the paper:** Missaoui, arXiv:2601.11495 (unrefereed single-author preprint) claims global C^{1,α} for Musielak–Orlicz equations, with t^p + a t^q (a ∈ C^{0,α}, q ≤ p+α) listed as an example. Check it; if it is correct and refereed, Table 1 and §3.2 can be updated.

## C. Counterexamples and Lavrentiev

28. **Esposito–Leonetti–Mingione 2004 (ABS + SEC).**
    - Claim (§3.3): regularity means higher integrability when q/p < (n+α)/n, and this is sharp.
    - Claim (§3.3): a counterexample with a gap and minimiser ∉ W^{1,q}_loc in the range p < n < n+α < q, extending Zhikov.
29. **Fonseca–Malý–Mingione 2004 (ABS).**
    - Claim: for 1 < p < N < N+α < q there is a W^{1,p} minimiser with singular set of Hausdorff dimension > N−p−ε.
    - DOI MEM (not included).
30. **Zhikov 1986, 1995 (MEM/SEC).**
    - Claim: checkerboard examples for a variable exponent (<2 in two opposite quadrants, >2 in the others) and for double phase with Lipschitz a and p < 2 < 3 < q (1995, Ex. 3.1).
    - Claim: log-Hölder continuity excludes the gap ("attributed to Zhikov, cf. [Zhikov1995]"). **Find the exact reference.**
    - Claim (§1): the double phase energy goes back to Zhikov 1986.
    - Citation data MEM.
31. **Balci–Diening–Surnachev 2020 (FT, arXiv v2 numbering).**
    - Thm 32: gap for q > p + α max{1, (p−1)/(d−1)}.
    - Thm 30: gap for an exponent with modulus 1/log(e+1/t)^κ, κ < 1.
    - **Theorem numbers may be shifted in the published version** (BGS cite it as [2, Thm 34]).
    - Article number MEM.
32. **Bulíček–Gwiazda–Skrzeczkowski 2022, Thm 2.3 (FT).**
    - Claim: bounded Lipschitz domain, a ∈ C^α(Ω̄), q ≤ p + α max{1, p/d}, data in W^{1,q}, no boundedness; then no gap; sharp for p ≤ d.
    - Pages/DOI MEM. The v1 statement has p and q swapped (typo).
33. **Borowski–Chlebicka–De Filippis–Miasojedow 2024 (ABS).**
    - Claim: a sharp scale of weights excluding the gap. Vol/art MEM.
34. **Giaquinta 1987; Marcellini 1989/1991; Hong 1992 (SEC via Hirsch–Schäffner).**
    - Claim: unbounded minimisers for anisotropic autonomous (p,q) integrands, e.g. p = 2, q = 4, n ≥ 6 (Giaquinta).
    - Hong's content was not verified at all.
35. **Ball–Mizel 1985 (MEM/SEC).**
    - Claim: a 1D strictly convex superlinear integrand whose singular minimiser does not satisfy the integrated EL equation, with a gap.
36. **Carstensen–Ortner 2010 (ABS).**
    - Claim (§5): conforming FEM fail exactly when the Lavrentiev phenomenon is present.
37. **Balci–Ortner–Storn 2022 (FT).**
    - Claim (§5): conforming discrete minima converge to the infimum over the closure of smooth functions, with exact or accurate quadrature.
    - Thm 1: Crouzeix–Raviart converges to min_W.
    - Claim: 2D experiments show an energy gap; the conforming scheme is sensitive to quadrature.
    - Vol/pages MEM.
38. **Ortner 2011 (ABS).**
    - Claim: Crouzeix–Raviart converges to the correct minimiser for convex problems.
39. **Balci–Kaltenbach arXiv:2303.10687 (ABS).**
    - Cited only as "see also". The journal is unknown.

## D. FE rates (§6.9)

40. **Barrett–Liu 1993 (ABS).**
    - Claim: order h in W^{1,p} for p < 2 and in W^{1,q}, q ∈ [1,2], for p > 2, for regular solutions.
41. **Liu–Barrett 1994 (FT-WF, exponents partly garbled).**
    - Claim: quasi-norm estimates.
    - Claim: the L² bound of Thm 4.3 is obtained via Hölder.
42. **Ebmeyer–Liu 2005 (ABS).**
    - Claim: quasi-norm interpolation and optimal estimates.
43. **Diening–Růžička 2007, Diening–Kreuzer 2008.**
    - Claim: order h for ‖F(∇u) − F(∇u_h)‖₂ when F(∇u) ∈ W^{1,2}, for Δ2 ∩ ∇2 N-functions.
    - Source: eq. (1.4) of Diening–Kreuzer (FT-WF), attributed there to Diening–Růžička.
44. **Del Pezzo–Martínez 2015 (FT-WF).**
    - Claim: order h^{p₁/2}, and order h under extra regularity, in W^{1,p(·)}, for p ≤ 2.
45. **Breit–Diening–Schwarzacher 2015 (FT-WF).**
    - Claim: order h^α in the quasi-norm for p ∈ C^{0,α}. DOI MEM (not included).
    - The earlier attribution to them for discontinuous coefficients has been REMOVED.
46. **Glowinski–Marrocco 1975 (FT-WF, OCR garbled).**
    - Claim: convergence without rate (Thm 4.1), and W^{1,p} rates for u ∈ W^{2,p} (Thm 6.1).
47. **Claim (§6.9): "We have found no proof of an O(h²) estimate in L² for the p-Laplacian or g-Laplacian".**
    - This is based on the librarian's search only. Please double-check with your own knowledge of the literature.

## E. Items about the code / data / the brief

48. **"7 p≈1 runs that did not reach tolerance" (brief, Part 5).**
    - This was not found in materials/README.md. The 7 runs are in paper_v2/results/convergence_order_true.csv (pxlaplace data), which is not used here.
    - main.tex instead reports a new glaplace run for p = 1.2 and 1.1, in which 8/8 runs did not converge in 6000 iterations.
    - After refereeing, the count in §6.9 is "twelve runs": 8 p→1, 2 expdead, and the two r = 1 runs of Table 7.
    - Confirm this is what you intended.
49. **Exact radial solutions "to ~1e−10".**
    - The paper only states my own check against the closed forms for p = 1.5, 2, 4.
    - The accuracy for the other G (in particular the dead core and the two-phase jump) was not independently checked.
50. **Software availability.**
    - Check that github.com/jfbonder/glaplace is public, that the license is MIT as in the README, and that the example scripts used here are in the repo.
    - The scripts in new/scripts/ (Lyapunov identity, extra runs, figures) are not in the repo; decide whether to add them.
51. **Remark 5.2 and Proposition 5.3 (deliberate deviation from Correction 6 of the brief).**
    - Remark 5.2 claims rigorously that (S) holds for the explicit radial data, including e^t−1−t and e^t−1, so Theorems A/B apply to those particular problems.
    - The proof is in the paper (Lemmas 2.5, 3.2, 5.1). **Decide whether to keep it.**
    - Prop. 5.3 (new) proves J(u_h) → J(u), ∇u_h → ∇u in measure, and H¹₀ convergence when g/t is nondecreasing, for the autonomous radial problems with shape-regular meshes. Its proof is self-contained. Check the interpolation step: the mean-value estimate |∇I_h u|_K − ∇u(x)| ≤ C_σ osc_K ∇u and the boundedness of ∇I_h u.
52. **Placeholders.** Author addresses and the funding acknowledgement are TODO placeholders in main.tex.

53. **Prop. 2.15(ii), converse (new).** The converse on Lipschitz domains uses the standard fact that H¹₀ ∩ {∇v ∈ L^p} ⊂ W^{1,p}₀ (zero trace plus W^{1,p} implies W^{1,p}₀ on Lipschitz domains). No reference is cited; add one if wanted.
54. **Stopping-criterion recommendation (new, §6.6).** We recommend monitoring the dual residual but did not rerun the experiments with it.

## F. Bibliographic data flagged MEM in REFERENCES_A/B but printed in the bibliography (G17 of the referee)

Verify each against the publisher record before submission:
- GPT16, pp. 19–94.
- Svaiter 2011, vol. 49 and pp. 280–287.
- Bùi–Combettes 2020, pp. 1118–1120.
- Balci–Ortner–Storn 2022, vol. 151, pp. 779–805.
- BGS 2022, pp. 209–240.
- BDS 2020, Art. 180.
- BCM 2018, Art. 62.
- De Filippis–Mingione 2021, vol. 242, pp. 973–1057.
- Hästö–Ok 2022, vol. 24(4), pp. 1285–1334.
- Cianchi–Maz'ya JEMS 16 (2014) 571–595.
- BCDM 2024, Art. 35.
- Gabay 1983, pp. 299–331.
- Acerbi–Mingione 2002, ARMA 164, 213–259.
- Zhikov 1986 and 1995 (all data).
- Fan 2007 (vol. 235(2), 397–417).
- Lavrentiev (not cited).
- Gabay–Mercier (not cited).

No DOIs are printed. The only one, Antonini's, was unverified and has been removed.
