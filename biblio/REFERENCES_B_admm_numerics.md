# References B: ADMM / augmented Lagrangian theory and FE numerics for p-, p(x)-, Orlicz Laplacians

Prepared 2026-10-12 for the paper on convergence of the Glowinski decomposition-coordination (DCM, ALG2, ADMM) method for the g-Laplacian.

## 0. Download status (read first)

**No PDF could be downloaded.** The sandbox egress policy refuses curl CONNECT (HTTP 403) to arxiv.org, export.arxiv.org, hal.science, web.stanford.edu, ucla.edu, optimization-online.org and every other academic host I tried. The proxy README says not to retry or route around such denials, so I did not. The WebFetch tool, which reads a page and summarises it through a model, could reach arXiv and several open-access hosts, so many texts below were **read in full text through WebFetch**. These are model-mediated readings and not my own line-by-line reading. Treat exact constants and equation numbers as "verify against the PDF before citing".

Other notes:
- `/home/claude/biblio/` was **empty** when I started. The Del Pezzo–Martínez arXiv:1401.2442 PDF that was said to be there was not.
- pub.uni-bielefeld.de returned a bot-protection page, and two api.crossref.org DOI lookups returned HTTP 429 (rate limited). I took those citations from search results and from L. Diening's publication page.

Verification labels used below:
- **FT-WF**: full text read via WebFetch.
- **ABS**: abstract or metadata page only.
- **MEM**: from memory, not verified in this session.

To fetch them yourself, every arXiv id below is at `https://arxiv.org/pdf/<id>`.

---

## TOPIC C: ADMM and augmented Lagrangian theory

### C1. Glowinski 1984, Chapter VI
```bibtex
@incollection{Glowinski1984ChVI,
  author    = {Glowinski, Roland},
  title     = {Decomposition-Coordination Methods by Augmented {L}agrangian: Applications},
  booktitle = {Numerical Methods for Nonlinear Variational Problems},
  series    = {Springer Series in Computational Physics},
  publisher = {Springer}, address = {Berlin, Heidelberg}, year = {1984},
  chapter   = {VI}, pages = {166--194}, doi = {10.1007/978-3-662-12613-4_6},
  note      = {Reprinted 2008, Scientific Computation series}
}
```
- Link: https://link.springer.com/chapter/10.1007/978-3-662-12613-4_6 (paywalled; preview disabled). **NOT DOWNLOADED.**
- Verified on the Springer page: pages 166–194, DOI as above. The abstract says the chapter treats min F(Bv)+G(v), with B ∈ L(V,H) and F, G convex, proper and l.s.c. It also says the chapter "follows Fortin and Glowinski".
- **Theorem content: NOT verified from the book itself.** See Section (a) below for what accessible sources reproduce.
- Status: ABS (chapter metadata only).

### C2. Fortin–Glowinski 1983 (eds.), Chapter III, and Gabay 1983, Chapter IX
```bibtex
@book{FortinGlowinski1983,
  editor    = {Fortin, Michel and Glowinski, Roland},
  title     = {Augmented {L}agrangian Methods: Applications to the Numerical Solution of Boundary-Value Problems},
  series    = {Studies in Mathematics and its Applications}, volume = {15},
  publisher = {North-Holland}, address = {Amsterdam}, year = {1983}
}
@incollection{FortinGlowinski1983ChIII,
  author    = {Fortin, Michel and Glowinski, Roland},
  title     = {On Decomposition-Coordination Methods Using an Augmented {L}agrangian},
  booktitle = {Augmented {L}agrangian Methods},
  series    = {Studies in Mathematics and its Applications}, volume = {15},
  publisher = {North-Holland}, year = {1983}, chapter = {III},
  note      = {ScienceDirect pii S0168-2024(08)70028-6; pages not verified}
}
@incollection{Gabay1983,
  author    = {Gabay, Daniel},
  title     = {Applications of the method of multipliers to variational inequalities},
  booktitle = {Augmented {L}agrangian Methods},
  series    = {Studies in Mathematics and its Applications}, volume = {15},
  publisher = {North-Holland}, year = {1983}, chapter = {IX}, pages = {299--331}
}
```
- Link: https://www.sciencedirect.com/bookseries/studies-in-mathematics-and-its-applications/vol/15/suppl/C (robots-blocked to WebFetch). **NOT DOWNLOADED.**
- Content:
  - Chapter III defines ALG1 and ALG2 and proves convergence for 0 < ρ < (1+√5)/2·r (MEM).
  - Tao–Yuan (arXiv:2006.08309, read FT-WF) confirm that the golden-ratio bound is attributed to Fortin–Glowinski [1983]. The only locators they give are "Eqs. (5.26)–(5.27)", which are the monotone decrease of a Lyapunov quantity.
  - Gabay's Chapter IX identifies ADMM with Douglas–Rachford on the dual. This is MEM, and the identification is also stated in Boyd et al. §3.5 (FT-WF).
- Status: MEM, plus indirect confirmation.

### C3. Glowinski–Le Tallec 1989
```bibtex
@book{GlowinskiLeTallec1989,
  author    = {Glowinski, Roland and Le Tallec, Patrick},
  title     = {Augmented {L}agrangian and Operator-Splitting Methods in Nonlinear Mechanics},
  series    = {SIAM Studies in Applied Mathematics}, volume = {9},
  publisher = {SIAM}, address = {Philadelphia}, year = {1989},
  doi       = {10.1137/1.9781611970838}
}
```
- **NOT DOWNLOADED** (paywalled). The DOI is from memory.
- Content (MEM): Chapter 3 restates the ALG1, ALG2 and ALG3 convergence theorems in Hilbert spaces. They give strong convergence of (u^n, p^n) and only weak convergence of λ^n, under assumptions of the type listed in (a).
- Status: MEM.

### C4. Glowinski–Pan–Tai 2016 (best accessible reproduction of the Glowinski hypotheses)
```bibtex
@incollection{GlowinskiPanTai2016,
  author    = {Glowinski, Roland and Pan, Tsorng-Whay and Tai, Xue-Cheng},
  title     = {Some Facts About Operator-Splitting and Alternating Direction Methods},
  booktitle = {Splitting Methods in Communication, Imaging, Science, and Engineering},
  editor    = {Glowinski, R. and Osher, S. J. and Yin, W.},
  publisher = {Springer}, year = {2016}, pages = {19--94},
  doi       = {10.1007/978-3-319-41589-5_2},
  note      = {Preprint: UCLA CAM Report 16-10}
}
```
- Link: https://ww3.math.ucla.edu/camreport/cam16-10.pdf (open access; curl blocked by policy). **NOT DOWNLOADED.** The page range 19–94 is MEM.
- Content (§3.2):
  - Setting: V and H are real Hilbert spaces, B ∈ L(V,H), and F, G are proper, l.s.c. and convex. Problem (39) is u = argmin_v [F(Bv)+G(v)].
  - Assumptions for ALG1, which ALG2 inherits:
    - (i) L has a saddle point.
    - (ii) B is injective and R(B) is closed in H.
    - (iii) F(q)/|q| → +∞ as |q| → ∞.
    - (iv) F = F₀+F₁, with F₀ and F₁ proper, l.s.c. and convex, and F₀ Gâteaux-differentiable and **uniformly convex on bounded sets** of H.
  - ALG1 with 0 < ρ < 2r gives (u^n, p^n) → (u, Bu) in V×H (eq. (46)).
  - ALG2 (eqs. (47)–(50)) does the p-step first, then the u-step, then λ^{n+1} = λ^n + ρ(Bu^n − p^n). If 0 < ρ < (1+√5)/2·r, the same conclusion holds. This is stated in the text, with no theorem number.
  - The multiplier sequence converges "no better than weakly in general".
  - Remark 21: ρ = r places the algorithms in the ADMM family.
  - §3.3: with V = H, B = I and ρ = r, ALG2 is a disguised Douglas–Rachford scheme with r = 1/Δt.
- Status: FT-WF.

### C5. Glowinski–Marrocco 1975
```bibtex
@article{GlowinskiMarrocco1975,
  author  = {Glowinski, R. and Marrocco, A.},
  title   = {Sur l'approximation, par {\'e}l{\'e}ments finis d'ordre un, et la r{\'e}solution, par p{\'e}nalisation-dualit{\'e}, d'une classe de probl{\`e}mes de {D}irichlet non lin{\'e}aires},
  journal = {RAIRO Analyse Num{\'e}rique}, volume = {9}, number = {R2}, year = {1975}, pages = {41--76},
  doi     = {10.1051/m2an/197509R200411}
}
```
- Link: https://www.esaim-m2an.org/articles/m2an/pdf/1975/02/m2an197509R200411.pdf (free; curl blocked). **NOT DOWNLOADED.**
- Content (the OCR was partly garbled):
  - Thm 4.1: P1 solutions converge strongly in W₀^{1,p} without a rate, under a minimum-angle condition.
  - Thm 6.1: W^{1,p}-norm rates when u ∈ W^{2,p}. Eq. (6.1) covers p > 2 and eq. (6.2) covers 1 < p ≤ 2.
  - Thm 8.1: fractional regularity W^{s,p}, 1 ≤ s ≤ 2, by nonlinear interpolation.
  - The exponents were illegible. MEM: the p ≥ 2 rate is usually quoted as h^{1/(p−1)}. The p < 2 exponent is uncertain.
  - §9, Thm 9.1: penalty-duality (ALG1-type) iteration **in the discrete finite-dimensional setting**. It gives u_h^n → u_h and y_h^n → ∇u_h for ρ_n in a bounded interval; the upper bound was garbled.
  - §10: a continuous dual iteration whose convergence is only claimed ("on démontrerait…").
  - §12: numerics on the disk for p from 1.1 to 50, with timings.
- Status: FT-WF (partly garbled).

### C6. Gabay–Mercier 1976
```bibtex
@article{GabayMercier1976,
  author  = {Gabay, Daniel and Mercier, Bertrand},
  title   = {A dual algorithm for the solution of nonlinear variational problems via finite element approximation},
  journal = {Computers \& Mathematics with Applications}, volume = {2}, number = {1}, year = {1976}, pages = {17--40},
  doi     = {10.1016/0898-1221(76)90003-1}
}
```
- **NOT DOWNLOADED** (Elsevier; open-archive status unknown). The DOI and pages are MEM.
- Content (MEM): ADMM convergence for 0 < ρ < 2r under strong-monotonicity and injectivity-type assumptions, in Hilbert space. It also covers the FE setting.
- Status: MEM.

### C7. Lions–Mercier 1979
```bibtex
@article{LionsMercier1979,
  author  = {Lions, P.-L. and Mercier, B.},
  title   = {Splitting algorithms for the sum of two nonlinear operators},
  journal = {SIAM J. Numer. Anal.}, volume = {16}, number = {6}, year = {1979}, pages = {964--979},
  doi     = {10.1137/0716071}
}
```
- **NOT DOWNLOADED.**
- Content, as reported by Svaiter and by Bùi–Combettes (FT-WF):
  - For A and B maximal monotone, the DR governing sequence z_k = x_k + λb_k converges **weakly**.
  - Lions–Mercier did not prove weak convergence of the shadow sequence x_k = J_{λB} z_k.
  - Strong convergence requires extra conditions, e.g. one operator uniformly or strongly monotone (MEM).
- Status: secondhand FT-WF.

### C8. Eckstein–Bertsekas 1992
```bibtex
@article{EcksteinBertsekas1992,
  author  = {Eckstein, Jonathan and Bertsekas, Dimitri P.},
  title   = {On the {D}ouglas--{R}achford splitting method and the proximal point algorithm for maximal monotone operators},
  journal = {Math. Programming}, volume = {55}, year = {1992}, pages = {293--318},
  doi     = {10.1007/BF01581204}
}
```
- Link: https://www.mit.edu/~dimitrib/Eckstein_Bertsekas_D-R.pdf (author page; curl blocked). Also LIDS-P-1919: https://dspace.mit.edu/bitstream/1721.1/3160/1/P-1919-20783099.pdf. **NOT DOWNLOADED.** The DOI is MEM.
- Content (FT-WF, with OCR noise in symbols):
  - **Thm 3** (generalized proximal point, Hilbert H): z^{k+1} = (1−ρ_k)z^k + ρ_k w^k, with ‖w^k − J_{c_kT} z^k‖ ≤ ε_k, Σε_k < ∞, 0 < inf ρ_k ≤ sup ρ_k < 2 and inf c_k > 0. If zer T ≠ ∅, z^k converges **weakly**. Otherwise z^k is unbounded.
  - **Thm 5**: characterises the fixed points of the DR operator.
  - **Thm 7** (generalized DR, Hilbert H, relaxation in (0,2), summable errors): z^k converges **weakly** to a point of Z* = {u+λb : b ∈ Bu, −b ∈ Au}.
  - **Thm 8** (generalized ADMM): **finite-dimensional**, R^n and R^m, with M of full column rank. With a Kuhn–Tucker pair, x^k → x*, w^k → Mx* and p^k → a dual solution. The penalty λ equals the dual step, with relaxation ρ_k ∈ (0,2).
- Status: FT-WF.

### C9. Boyd–Parikh–Chu–Peleato–Eckstein 2011
```bibtex
@article{Boyd2011ADMM,
  author  = {Boyd, Stephen and Parikh, Neal and Chu, Eric and Peleato, Borja and Eckstein, Jonathan},
  title   = {Distributed Optimization and Statistical Learning via the Alternating Direction Method of Multipliers},
  journal = {Foundations and Trends in Machine Learning}, volume = {3}, number = {1}, year = {2011}, pages = {1--122},
  doi     = {10.1561/2200000016}
}
```
- Link: https://web.stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf (author page; curl blocked). **NOT DOWNLOADED.**
- §3.2, verified FT-WF:
  - Assumption 1: f and g are closed, proper and convex.
  - Assumption 2: the unaugmented Lagrangian L₀ has a saddle point. No rank conditions are needed.
  - Conclusions: the residual r^k → 0, the objective tends to p*, and the dual iterates y^k → y*.
  - Primal iterates x^k, z^k are **not** claimed to converge without additional assumptions.
  - The dual step equals ρ, the penalty parameter.
  - The setting is **finite-dimensional**.
- Appendix A body (MEM; WebFetch truncated the text before it):
  - V^k = (1/ρ)‖y^k−y*‖² + ρ‖B(z^k−z*)‖².
  - (A.1): V^{k+1} ≤ V^k − ρ‖r^{k+1}‖² − ρ‖B(z^{k+1}−z^k)‖².
  - (A.2) and (A.3) are two-sided bounds on p^{k+1} − p*.
  - The proof uses only telescoping. It needs neither Opial nor compactness, but compactness is implicit because the setting is finite-dimensional.
- Status: §3.2 FT-WF; Appendix A MEM.

### C10. Svaiter 2011
```bibtex
@article{Svaiter2011,
  author  = {Svaiter, B. F.},
  title   = {On weak convergence of the {D}ouglas--{R}achford method},
  journal = {SIAM J. Control Optim.}, volume = {49}, number = {1}, year = {2011}, pages = {280--287},
  doi     = {10.1137/100788100}
}
```
- arXiv:1007.2173. Also https://optimization-online.org/wp-content/uploads/2010/07/2676.pdf. **NOT DOWNLOADED.** The volume and pages are MEM; the DOI was confirmed by search.
- Content, **Theorem 1** (real Hilbert H; A and B maximal monotone; λ > 0; summable errors α_k, β_k; zer(A+B) ≠ ∅): the iterates (x_k, b_k) and (y_k, a_k) converge **weakly** to (x̄, b̄) and (x̄, −b̄), with b̄ ∈ Bx̄ and −b̄ ∈ Ax̄. This upgrades Lions–Mercier, who had weak convergence only of x_k + λb_k.
- Status: FT-WF.

### C11. Bùi–Combettes 2020
```bibtex
@article{BuiCombettes2020,
  author  = {B{\`u}i, Minh N. and Combettes, Patrick L.},
  title   = {The {D}ouglas--{R}achford Algorithm Converges Only Weakly},
  journal = {SIAM J. Control Optim.}, volume = {58}, number = {2}, year = {2020}, pages = {1118--1120},
  doi     = {10.1137/19M1308451}
}
```
- arXiv:1912.09564. **NOT DOWNLOADED.** The volume and pages are MEM; the DOI was confirmed by search.
- Content:
  - **Theorem 1** recalls DR: y_n ⇀ y (Lions–Mercier) and x_n = J_B y_n ⇀ x ∈ zer(A+B) (Svaiter; also Bauschke–Combettes, Thm 26.11(iii)).
  - **Counterexample 2** (separable infinite-dimensional H): A = N_V for a closed subspace V, and B = (P_V P_K P_V)^{-1} − Id for a closed convex cone K built from Hundal's curve. Then x_n ⇀ 0, but x_n does **not** converge in norm.
  - **Counterexample 4** gives the same failure for Spingarn's partial inverses.
- My remark, not stated in the paper: P_V P_K P_V is the gradient of a convex function, so B is a subdifferential and A = ∂ι_V. The counterexample therefore lives within convex minimization.
- Status: FT-WF.

### C12. Bauschke–Combettes book, and Bauschke–Moursi 2017
```bibtex
@book{BauschkeCombettes2017,
  author    = {Bauschke, Heinz H. and Combettes, Patrick L.},
  title     = {Convex Analysis and Monotone Operator Theory in {H}ilbert Spaces},
  edition   = {2}, series = {CMS Books in Mathematics}, publisher = {Springer}, year = {2017},
  doi       = {10.1007/978-3-319-48311-5}
}
@article{BauschkeMoursi2017,
  author  = {Bauschke, Heinz H. and Moursi, Walaa M.},
  title   = {On the {D}ouglas--{R}achford algorithm},
  journal = {Math. Program.}, volume = {164}, year = {2017}, pages = {263--284},
  doi     = {10.1007/s10107-016-1086-3}
}
```
- Bauschke–Combettes: https://pcombet.math.ncsu.edu/livre1.pdf exposes only the 1st-edition front matter. In the 1st edition, DR is §25.2 (p. 366) and §27.2 (p. 400). **NOT DOWNLOADED.**
- Bauschke–Combettes, 2nd edition, **Thm 26.11** (MEM):
  - Hypotheses: A and B maximally monotone, zer(A+B) ≠ ∅, γ > 0, λ_n ∈ [0,2] with Σλ_n(2−λ_n) = ∞.
  - (i)–(ii): x_n ⇀ a fixed point of the DR operator.
  - (iii): the shadow J_{γB} x_n ⇀ a zero of A+B. Bùi–Combettes cite exactly 26.11(iii).
  - Later items: **strong** convergence of the shadow if A or B is uniformly monotone on bounded sets (or demiregular).
  - The Fejér/Opial tool is Lemma 2.47 and Thm 5.5 (Fejér monotone with weak cluster points in C gives weak convergence). These numbers are MEM and uncertain.
- Bauschke–Moursi, arXiv:1604.04603: weak convergence of the shadow sequence in the inconsistent convex-feasibility case, plus a general sufficient condition. The journal data and DOI are MEM.
- Status: ABS / MEM.

### C13. Moradifam–Nachman (split Bregman in Hilbert space)
```bibtex
@misc{MoradifamNachman2011,
  author = {Moradifam, Amir and Nachman, Adrian},
  title  = {Convergence of the alternating split {B}regman algorithm in infinite-dimensional {H}ilbert spaces},
  year   = {2011}, eprint = {1112.1960}, archivePrefix = {arXiv}
}
```
- **NOT DOWNLOADED.** The journal version is not identified.
- Content: **weak** convergence of alternating split Bregman (ADMM) in infinite-dimensional Hilbert spaces, including an inexact version. The proof goes through the dual, Svaiter's DR theorem and Setzer's ADMM–DR equivalence.
- Status: ABS.

### C14. Bartels–Milicevic (ADMM with variable step sizes)
```bibtex
@misc{BartelsMilicevic2017,
  author = {Bartels, S{\"o}ren and Milicevic, Marijo},
  title  = {Alternating direction method of multipliers with variable step sizes},
  year   = {2017}, eprint = {1704.06069}, archivePrefix = {arXiv}
}
```
- **NOT DOWNLOADED.** Journal version not confirmed.
- Setting: Hilbert spaces X and Y, B ∈ L(X,Y), F and G convex, proper and l.s.c., coercive. Assumption 3.2: a saddle point exists.
- **Thm 3.7**: for non-increasing steps τ_j bounded below, the energy D_j² = ‖δ_λ^j‖² + τ_j²‖Bδ_u^j‖² decreases. This is the user's V_n up to scaling, with ρ = r = τ. The residuals R_j are square-summable.
- Cor 3.9: strong convexity of G gives an error bound by the residual.
- Prop 3.11 and Cor 3.12: R_J² = O(1/J).
- Thm 3.14 and Prop 3.15: linear convergence when G is strongly convex with Lipschitz gradient and B′ has a bounded left inverse.
- Numerics: ROF (TV) with P1 and P0 elements. The constants are **h-dependent**: α_B ≈ h/2 and the optimal τ ~ h^{-1}. The multiplier scales like h^{-d/2}.
- The p-Laplacian is mentioned only as motivation, and p < 2 is not discussed.
- Status: FT-WF.

### C15. Other C3 items
- **Bartels**, *Numerical Methods for Nonlinear Partial Differential Equations*, Springer SSCM 47, 2015, DOI 10.1007/978-3-319-13797-1. Status: MEM. It contains a Hilbert-space ADMM convergence proof with an energy argument (exact section and theorem numbers not verified). **NOT DOWNLOADED.**
- **Bartels–Diening–Nochetto**, "Unconditional stability of semi-implicit discretizations of singular flows", SIAM J. Numer. Anal. 56 (2018) 1896–1914, DOI 10.1137/17M1159166 (DOI from search; volume and pages MEM), arXiv:1712.03260.
  - Proves unconditional energy stability of a semi-implicit scheme for singular p-Laplace flows (p < 2), with a regularised singular factor treated explicitly.
  - Error estimates degrade as the regularisation parameter shrinks.
  - It is not an ADMM analysis.
  - Status: ABS.
- **Hintermüller–Stadler**, SIAM J. Sci. Comput. 28 (2006) 1–23, DOI 10.1137/040613263. This is a semismooth Newton primal-dual method, not ADMM. Open copy: https://math.nyu.edu/~stadler/papers/HintermullerStadler06.pdf. Status: ABS.
- **Chambolle–Pock**, J. Math. Imaging Vis. 40 (2011) 120–145, DOI 10.1007/s10851-010-0251-1. A primal-dual method with finite-dimensional theory. Status: MEM.
- **Ito–Kunisch**, *Lagrange Multiplier Approach to Variational Problems and Applications*, SIAM 2008, DOI 10.1137/1.9780898718614. It has augmented Lagrangian chapters in Hilbert space. Status: MEM.
- **Diening–Fornasier–Tomasi–Wank** and **Balci–Diening–Storn**: see D-section entries D11 and D12 (Kačanov iterations, an alternative to ADMM).
- **No source found** that proves convergence of ADMM, DCM or split Bregman for the p-Laplacian with p < 2 in the continuous W^{1,p} setting. The closest are:
  - Glowinski–Marrocco §10: a claimed L^p–L^{p'} dual iteration.
  - Aragón–Fernández Bonder–Rubio (D2): they restrict to p ≥ 2 precisely so that L^p ⊂ L² and W^{1,p}_0 ⊂ H^1_0.
  - Del Pezzo–Martínez (D1): finite-dimensional DG.
  This is a genuine gap the paper can address.

---

## TOPIC D: FE error analysis and numerics

### D1. Del Pezzo–Martínez, DCM for p(x) (arXiv:1401.2442)
```bibtex
@misc{DelPezzoMartinez2014DCM,
  author = {Del Pezzo, Leandro M. and Mart{\'i}nez, Sandra},
  title  = {The decomposition-coordination method for the $p(x)$-{L}aplacian},
  year   = {2014}, eprint = {1401.2442}, archivePrefix = {arXiv}
}
```
- **NOT DOWNLOADED.** No journal version found.
- Setting:
  - V = S^k(T_h) (DG) and H = S^l × S^l with l ≥ k−1.
  - Bv = R_h(v) + ∇v, where R_h is the lifting operator.
  - F(q) = ∫|q|^{p(x)} (later ∫|q|^{p(x)}/p(x)). G collects the quadratic jump and boundary penalties, weighted by h^{−2/p′}.
- Hypotheses, reproduced as those of [Glowinski, Ch. VI]:
  - (H1) F and G are l.s.c., and dom(F∘B) ∩ dom G ≠ ∅.
  - (H2) F is convex and Gâteaux-differentiable, with F(q)/|q| → ∞.
  - (H3) R(B) is closed.
  - (H4) B is injective.
  - The authors note that H4 **fails** for their B, because of its kernel, and say H4 is used only to get u^n → u.
- Lemma 4.1 corresponds to [Gl, VI, Thm 2.1]: a saddle point gives a solution with Bu = η.
- **Thm 1.3** (ALG1): saddle point and 0 < α₀ ≤ ρ_n ≤ α₁ < 2r.
- **Thm 1.4** (ALG2): saddle point and 0 < ρ < r(1+√5)/2.
- Conclusions of both: u_n → u in V, η^n → Bu, λ^{n+1}−λ^n → 0, and λ^n bounded.
- **Remark 4.4** cites [Gl, VI, Thms 4.1 and 5.1] for ‖Bu^n−η^n‖ → 0, η^n → η and λ^n bounded under H1–H3. Strong convergence of u^n then comes from convexity and differentiability of G (eq. 4.13).
- Typos noted: the proof of Thm 1.4 is headed "Proof of Thm 1.3", and there is a u vs u_h confusion.
- Numerics:
  - Setup: [−1,1]², P0 DG, r = ρ = 1, p(x) = 1 + (b(x₁+x₂)/2 + 1 + b)^{-1}.
  - The L² errors decrease by about 4.1× from m = 100 to m = 2916 for b = 0. My estimate is an order of about 0.8, not stated by the authors.
  - Iterations grow with b, from 3 to 63.
- **Everything here is finite-dimensional.**
- Status: FT-WF.

### D2. Aragón–Fernández Bonder–Rubio 2023
```bibtex
@article{AragonFBRubio2023,
  author  = {Arag{\'o}n, Adriana and Fern{\'a}ndez Bonder, Juli{\'a}n and Rubio, Diana},
  title   = {Effective numerical computation of $p(x)$-{L}aplace equations in {2D}},
  journal = {Int. J. Comput. Math.}, volume = {100}, number = {11}, year = {2023}, pages = {2111--2123},
  doi     = {10.1080/00207160.2023.2263103}, eprint = {2204.08054}, archivePrefix = {arXiv}
}
```
- **NOT DOWNLOADED.** The journal data were verified via Crossref.
- Scheme (2.1)–(2.3):
  - −Δu_n = div(η_n − ν_{n−1}) + f.
  - |ν_n|^{p−2}ν_n + ν_n = η_n + ∇u_n, solved by scalar bisection on t^{p−1}+t = |η+∇u|.
  - η_{n+1} = η_n + ∇u_n − ν_n.
  - Discretisation: P1 for u, P0 for ∇u, η, ν and p, with one fixed stiffness matrix.
- **Thm 4.2** is Glowinski's Ch. VI Thm 2.1 (saddle point characterisation).
- **Thm 4.3**: for **p(x) ≥ 2**, with the proof "exactly Thm 5.1 of [Gl, Ch. VI]":
  - u_n → u* strongly in H¹₀.
  - ν_n → ν* strongly in L².
  - η_{n+1} − η_n → 0.
  - η_n bounded in L².
- **Thm 4.4**: r = 1 recovers (2.1)–(2.3).
- **Remark 4.5** calls the restriction p ≥ 2 "of technical nature".
- Numerics: one 100×100 mesh, sup and L^p errors, with no rates reported.
  - p = 20: error ≈ 2e−3 / 1e−3.
  - p = 1.1: error ≈ 2e−4 / 3e−4.
  - p(x) ∈ [4/3, 2]: error ≈ 1.5e−5 / 2.7e−4.
  - Qualitative torsion examples, including a discontinuous exponent (1.2 | 4).
- Minor sign and denominator typos were flagged.
- Status: FT-WF.

### D3. Del Pezzo–Martínez 2015, order of convergence for p(x)
```bibtex
@article{DelPezzoMartinez2015,
  author  = {Del Pezzo, Leandro M. and Mart{\'i}nez, Sandra},
  title   = {Order of convergence of the finite element method for the $p(x)$-{L}aplacian},
  journal = {IMA J. Numer. Anal.}, volume = {35}, number = {4}, year = {2015}, pages = {1864--1887},
  doi     = {10.1093/imanum/dru050}, eprint = {1311.3230}, archivePrefix = {arXiv}
}
```
- **NOT DOWNLOADED.** The journal data were verified via Crossref. The task's guess of "Bustos" is wrong: the authors are only Del Pezzo and Martínez.
- Setting: Ω ⊂ R² convex; P1 elements on polygonal Ω^h ⊂ Ω; 1 < p₁ ≤ p(x) ≤ p₂ ≤ 2, log-Hölder; f ∈ L^{q(·)} with q₁ > 2; g ∈ H².
- **Thm 1.1**: ‖u−u^h‖_{W^{1,p(·)}(Ω^h)} ≤ C h^{p₁/2}. This uses H² regularity (their Thm 2.9, which needs convex Ω and Lipschitz p).
- **Thm 1.2**: order **h** in W^{1,p(·)}, under two conditions:
  - (1.4) ∫|∇u|^{p−2}|D²u|² < ∞.
  - (1.5) u ∈ C^{2,α⁺}(τ) on each triangle, with α⁺ = (2−p⁺)/p⁺.
  - These are verifiable in radial examples.
- Numerics:
  - Setup: [−1,1]², the same p(x) family as D1, DCM with r = ρ = 1.
  - The fitted W^{1,p(·)} order is ≈ 0.95–1.00 for all b.
  - The constant grows like 1/(p₁−1).
  - **No L² rates are given.**
- Status: FT-WF.

### D4. Del Pezzo–Lombardi–Martínez 2012 (DG for p(x))
```bibtex
@article{DelPezzoLombardiMartinez2012,
  author  = {Del Pezzo, Leandro M. and Lombardi, Ariel L. and Mart{\'i}nez, Sandra},
  title   = {Interior penalty discontinuous {G}alerkin {FEM} for the $p(x)$-{L}aplacian},
  journal = {SIAM J. Numer. Anal.}, volume = {50}, number = {5}, year = {2012}, pages = {2497--2521},
  doi     = {10.1137/110820324}
}
```
- **NOT DOWNLOADED.** Data verified via Crossref.
- Content: convergence of IP-DG (strong convergence of the discrete gradients) with **no rate**. This is according to Breit–Diening–Schwarzacher (FT-WF) and D1.
- Status: secondhand.

### D5. Breit–Diening–Schwarzacher 2015 (P1 for p(·))
```bibtex
@article{BreitDieningSchwarzacher2015,
  author  = {Breit, Dominic and Diening, Lars and Schwarzacher, Sebastian},
  title   = {Finite element approximation of the $p(\cdot)$-{L}aplacian},
  journal = {SIAM J. Numer. Anal.}, volume = {53}, number = {1}, year = {2015}, pages = {551--572},
  doi     = {10.1137/130946046}, eprint = {1311.5121}, archivePrefix = {arXiv}
}
```
- **NOT DOWNLOADED.** The DOI is MEM; the rest was verified on Diening's page.
- **Thm 4.4**: p ∈ C^{0,α}(Ω̄) with p⁻ > 1, and F(·,∇v) ∈ W^{1,2}, where F(x,ξ) = (κ+|ξ|)^{(p(x)−2)/2}ξ. Then ‖F(·,∇v) − F(·,∇v_h)‖₂ ≤ c h^α.
- **Thm 4.8**: the same rate for a piecewise-frozen exponent p_T.
- Remark 4.5: for α < 1, F ∈ W^{1,2} is not expected. A Nikolskii N^{α,2} regularity is conjectured.
- No L² or W^{1,p} rates for u − u_h, and no numerics.
- Status: FT-WF.

### D6. Barrett–Liu 1993
```bibtex
@article{BarrettLiu1993,
  author  = {Barrett, John W. and Liu, W. B.},
  title   = {Finite element approximation of the $p$-{L}aplacian},
  journal = {Math. Comp.}, volume = {61}, number = {204}, year = {1993}, pages = {523--537},
  doi     = {10.1090/S0025-5718-1993-1192966-4}
}
```
- Link: https://www.ams.org/journals/mcom/1993-61-204/S0025-5718-1993-1192966-4/ (free access; the PDF fetch failed). **NOT DOWNLOADED.**
- Content (abstract):
  - P1 elements on a polygonal approximation, planar domain, p ∈ (1,∞).
  - For sufficiently regular solutions, which can be achieved for a subclass of data, the bounds are "optimal error bounds".
  - The error is measured in W^{1,q}, with **q = p when p < 2** and **q ∈ [1,2] when p > 2**. This means order h in W^{1,p} for p < 2, and in W^{1,q} (q ≤ 2) for p > 2.
  - Numerical results are included.
- Status: ABS.

### D7. Liu–Barrett 1994 (quasi-norm)
```bibtex
@article{LiuBarrett1994,
  author  = {Liu, W. B. and Barrett, John W.},
  title   = {Quasi-norm error bounds for the finite element approximation of some degenerate quasilinear elliptic equations and variational inequalities},
  journal = {RAIRO Mod{\'e}l. Math. Anal. Num{\'e}r.}, volume = {28}, number = {6}, year = {1994}, pages = {725--744}
}
```
- Link: https://www.numdam.org/item/M2AN_1994__28_6_725_0.pdf (free). **NOT DOWNLOADED.**
- Content:
  - Assumptions (A1)–(A2) on k(t)t, with parameters α₁ and α₂; the p-Laplacian has α₁ = α₂ = 1.
  - The quasi-norm satisfies "all the properties of a norm except the homogeneity property".
  - **Thm 3.1**: a quasi-norm best-approximation bound, with an L² term and an h² term.
  - **Thm 4.1 / Cor 4.1** (p ≤ 2): optimal W^{1,p} rate for u ∈ W^{2,2+α₁(p−2)}.
  - **Thm 4.2** (p ≥ 2): optimal H¹ rate.
  - **Thm 4.3**: local W^{1,r}(D) rate h^s, plus ‖u−u_h‖_{L²} ≤ C h^s via Hölder.
  - Thm 4.4: obstacle problem.
  - No numerics.
- Status: FT-WF (exponents partly garbled).

### D8. Ebmeyer–Liu 2005; Liu–Yan 2001; Hirn 2013
```bibtex
@article{EbmeyerLiu2005,
  author  = {Ebmeyer, Carsten and Liu, WenBin},
  title   = {Quasi-norm interpolation error estimates for the piecewise linear finite element approximation of $p$-{L}aplacian problems},
  journal = {Numer. Math.}, volume = {100}, number = {2}, year = {2005}, pages = {233--258},
  doi     = {10.1007/s00211-005-0594-5}
}
@article{LiuYan2001,
  author  = {Liu, W. B. and Yan, N. N.},
  title   = {Quasi-norm local error estimators for $p$-{L}aplacian},
  journal = {SIAM J. Numer. Anal.}, volume = {39}, number = {1}, year = {2001}, pages = {100--127}
}
@article{Hirn2013,
  author  = {Hirn, Adrian},
  title   = {Finite element approximation of singular power-law systems},
  journal = {Math. Comp.}, volume = {82}, year = {2013}, pages = {1247--1268}
}
```
- Ebmeyer–Liu (data verified via the Kent repository; **NOT DOWNLOADED**):
  - Quasi-norm interpolation estimates for several interpolators.
  - Optimal a priori P1 bounds under weaker regularity, for the 2D and 3D p-Laplacian.
  - Status: ABS.
- Liu–Yan: a posteriori quasi-norm estimators. Status: MEM.
- Hirn: P1 for p < 2 systems, with an h^{?} rate in the F-quasi-norm. Status: MEM, details not verified.
- All three **NOT DOWNLOADED.**

### D9. Diening–Růžička 2007; Diening–Kreuzer 2008; Belenki–Diening–Kreuzer 2012; Belenki–Berselli–Diening–Růžička 2012; Diening–Kreuzer–Süli 2013; Diening–Ettwein 2008
```bibtex
@article{DieningRuzicka2007,
  author  = {Diening, L. and R{\r u}{\v z}i{\v c}ka, M.},
  title   = {Interpolation operators in {O}rlicz--{S}obolev spaces},
  journal = {Numer. Math.}, volume = {107}, number = {1}, year = {2007}, pages = {107--129},
  doi     = {10.1007/s00211-007-0079-9}
}
@article{DieningKreuzer2008,
  author  = {Diening, Lars and Kreuzer, Christian},
  title   = {Linear convergence of an adaptive finite element method for the $p$-{L}aplacian equation},
  journal = {SIAM J. Numer. Anal.}, volume = {46}, number = {2}, year = {2008}, pages = {614--638},
  doi     = {10.1137/070681508}
}
@article{BelenkiDieningKreuzer2012,
  author  = {Belenki, L. and Diening, L. and Kreuzer, Ch.},
  title   = {Optimality of an adaptive finite element method for the $p$-{L}aplacian equation},
  journal = {IMA J. Numer. Anal.}, volume = {32}, number = {2}, year = {2012}, pages = {484--510},
  doi     = {10.1093/imanum/drr016}
}
@article{BBDR2012,
  author  = {Belenki, L. and Berselli, L. C. and Diening, L. and R{\r u}{\v z}i{\v c}ka, M.},
  title   = {On the finite element approximation of $p$-{S}tokes systems},
  journal = {SIAM J. Numer. Anal.}, volume = {50}, number = {2}, year = {2012}, pages = {373--397},
  doi     = {10.1137/10080436X}
}
@article{DieningKreuzerSuli2013,
  author  = {Diening, L. and Kreuzer, Ch. and S{\"u}li, E.},
  title   = {Finite element approximation of steady flows of incompressible fluids with implicit power-law-like rheology},
  journal = {SIAM J. Numer. Anal.}, volume = {51}, number = {2}, year = {2013}, pages = {984--1015}
}
@article{DieningEttwein2008,
  author  = {Diening, L. and Ettwein, F.},
  title   = {Fractional estimates for non-differentiable elliptic systems with general growth},
  journal = {Forum Math.}, volume = {20}, number = {3}, year = {2008}, pages = {523--556}
}
```
- Data verified via Crossref (Diening–Růžička, BBDR) and Diening's publication page (the others). The Diening–Kreuzer DOI is MEM.
- Links: Diening–Kreuzer preprint https://www.math.lmu.de/~diening/archive/diening_kreuzer-AFEM2007.pdf (FT-WF). **None downloaded.**
- Diening–Kreuzer content (FT-WF):
  - The function is F(a) = (φ′(|a|)/|a|)^{1/2} a, with F(a) = (κ+|a|)^{(p−2)/2}a for the p-Laplacian. Shifted N-functions φ_a are defined by φ_a′(t)/t = φ′(a+t)/(a+t) (eq. 3.2).
  - **Lemma 16**: J(u₁) − J(u₂) ~ ‖F(∇u₁) − F(∇u₂)‖₂².
  - **Eq. (1.4)**, attributed to [Diening–Růžička 2007]: the quasi-best-approximation ‖F(∇u) − F(∇u_h)‖₂² ≤ c min_{v_h} ‖F(∇u) − F(∇v_h)‖₂² ≤ c h² ‖∇F(∇u)‖₂². This gives **order h in the F-quasi-norm when F(∇u) ∈ W^{1,2}**, for all 1 < p < ∞ and general Δ₂ ∩ ∇₂ N-functions φ.
  - **Thm 20 / Cor 21**: linear convergence of AFEM with Dörfler marking.
- BBDR 2012 (MEM): the same F-quasi-norm order h for p-Stokes velocity, and reduced order for the pressure, h^{2/p′}.
- Diening–Ettwein (MEM): Nikolskii/fractional regularity, F(∇u) ∈ N^{s,2}, for general-growth systems.
- Status: ABS / FT-WF / MEM as indicated.

### D10. Diening–Stroffolini–Verde; Balci–Diening–Surnachev
- These are regularity papers, not FE papers. Diening's page lists:
  - Manuscripta Math. 129 (2009), DOI 10.1007/s00229-009-0277-0.
  - J. Differential Equations 253 (2012) 1943–1958, DOI 10.1016/j.jde.2012.06.010.
  - ESAIM COCV 17 (2011) 178–189.
- Balci–Diening–Surnachev, "New examples on Lavrentiev gap using fractals", Calc. Var. PDE 59 (2020): Lavrentiev examples for double phase. Status: MEM.

### D11. Diening–Fornasier–Tomasi–Wank 2020 (relaxed Kačanov, p ≤ 2)
```bibtex
@article{DieningFornasierTomasiWank2020,
  author  = {Diening, Lars and Fornasier, Massimo and Tomasi, Roberto and Wank, Maximilian},
  title   = {A relaxed {K}a{\v c}anov iteration for the $p$-{P}oisson problem},
  journal = {Numer. Math.}, volume = {145}, year = {2020}, pages = {1--34},
  doi     = {10.1007/s00211-020-01107-1}, eprint = {1702.03844}, archivePrefix = {arXiv}
}
```
- The DOI was confirmed via search. The volume, pages and the fourth author are MEM; the arXiv v2 page lists only Diening, Fornasier and Wank.
- Content: an iteratively reweighted, relaxed Kačanov scheme for 1 < p ≤ 2, with at least an **algebraic** convergence rate.
- Status: ABS.

### D12. Balci–Diening–Storn 2023 (relaxed Kačanov, large p)
```bibtex
@article{BalciDieningStorn2023,
  author  = {Balci, Anna Kh. and Diening, Lars and Storn, Johannes},
  title   = {Relaxed {K}a{\v c}anov scheme for the $p$-{L}aplacian with large exponent},
  journal = {SIAM J. Numer. Anal.}, volume = {61}, number = {6}, year = {2023}, pages = {2775--2794},
  doi     = {10.1137/22M1528550}, eprint = {2210.06402}, archivePrefix = {arXiv}
}
```
- Data verified via Crossref.
- Content: globally convergent for 2 ≤ p < ∞. Each step solves a weighted Poisson problem. The rate is "independent of the underlying mesh".
- Status: ABS.

### D13. Double phase / Musielak–Orlicz FE
```bibtex
@article{BalciOrtnerStorn2022,
  author  = {Balci, Anna Kh. and Ortner, Christoph and Storn, Johannes},
  title   = {Crouzeix--{R}aviart finite element method for non-autonomous variational problems with {L}avrentiev gap},
  journal = {Numer. Math.}, volume = {151}, year = {2022}, pages = {779--805},
  doi     = {10.1007/s00211-022-01303-1}, eprint = {2106.06837}, archivePrefix = {arXiv}
}
```
- The DOI was confirmed via search; the volume and pages are MEM.
- Content: "conforming schemes fail due to the Lavrentiev gap". Crouzeix–Raviart converges to the global minimiser for non-autonomous non-standard growth integrands, the double phase case included. No rates appear in the abstract.
- Status: ABS.
- **Highly relevant warning for the double-phase experiments.** When a Lavrentiev gap exists, conforming P1 (which our DCM uses) may converge to the W^{1,∞}-closure minimiser rather than the true one. Our smooth test cases presumably have no gap; say so explicitly.
- Also relevant: Heid, arXiv:2304.01566, "A damped Kačanov scheme for … relaxed p(x)-Poisson". It proves convergence of the relaxed and discrete solutions, with no rates. Journal not identified. Status: ABS.
- I found **no paper with a priori FE rates specifically for double phase**. arXiv:2410.14995, "Discarding Lavrentiev's gap…", was only seen as a search hit and was not read.

### D14. Chen–Levine–Rao 2006
```bibtex
@article{ChenLevineRao2006,
  author  = {Chen, Yunmei and Levine, Stacey and Rao, Murali},
  title   = {Variable exponent, linear growth functionals in image restoration},
  journal = {SIAM J. Appl. Math.}, volume = {66}, number = {4}, year = {2006}, pages = {1383--1406},
  doi     = {10.1137/050624522}
}
```
- Status: MEM. Not verified, not downloaded. To my recollection it analyses a variable-exponent (TV to p(x)) model with a gradient-flow-type numerical scheme, **not** the DCM. Cite only as a motivating application.

---

## Answers

### (a) Hypotheses of Glowinski's Ch. VI Thm 5.1 (ALG2) as reproduced by accessible sources

No accessible source quotes Thm 5.1 verbatim. Three sources reproduce it, and they partly disagree.

1. **Glowinski–Pan–Tai 2016 (C4)**, the most complete and written by Glowinski himself:
   - Assumptions:
     - V and H are Hilbert spaces and B ∈ L(V,H).
     - F and G are proper, l.s.c. and convex.
     - (i) A saddle point of L_r exists.
     - (ii) B is injective and R(B) is closed.
     - (iii) F(q)/|q| → ∞.
     - (iv) F = F₀+F₁, with F₀ Gâteaux-differentiable and uniformly convex on bounded sets.
   - ALG2 condition: 0 < ρ < (1+√5)/2·r.
   - Conclusion: (u^n, p^n) → (u, Bu) in V×H (strongly, by context). λ^n converges only weakly in general.
   - The p-step comes before the u-step.
2. **Del Pezzo–Martínez (D1)** lists (H1)–(H4):
   - (H1) l.s.c. and the domain condition.
   - (H2) F convex and Gâteaux-differentiable, with F(q)/|q| → ∞.
   - (H3) R(B) closed.
   - (H4) B injective.
   - Under H1–H3 they get Bu^n − η^n → 0, η^n → η and λ^n bounded (citing Thms 4.1 and 5.1). H4 is used only for u^n → u.
   - They do **not** list uniform convexity. This is probably a simplification, since uniform convexity is what gives strong convergence of p^n in infinite dimensions. In their finite-dimensional DG setting it is irrelevant.
3. **Aragón–FB–Rubio (D2)** applies Thm 5.1 with V = H¹₀ and H = L², for p ≥ 2. The conclusions are u_n → u* in H¹₀ strongly, ν_n → ν* in L² strongly, η_{n+1} − η_n → 0 and η_n bounded. **Weak** convergence of η_n is not stated.

**Recommendation:** cite the hypotheses as in Glowinski–Pan–Tai §3.2.2. Obtain the book (Springer, DOI 10.1007/978-3-662-12613-4_6) to confirm the exact numbering: Thm 4.1 for ALG1, Thm 5.1 for ALG2, and Remark 4.4. The Del Pezzo–Martínez citation pattern supports this numbering, but I could not verify it.

### (b) Is "saddle point + closed range ⇒ strong convergence" true for ADMM in infinite dimensions?

**No, not in general. It needs an additional uniform-convexity or uniform-monotonicity hypothesis.**

- **General case: weak convergence only.**
  - Eckstein–Bertsekas Thms 3 and 7, Svaiter Thm 1, Bauschke–Combettes Thm 26.11 and Moradifam–Nachman all give **weak** convergence of the DR and ADMM iterates in Hilbert space.
  - Bùi–Combettes 2020 give a counterexample where DR fails to converge strongly. It is built from subdifferentials of convex functions (A = N_V; B a subdifferential by my remark in C11).
  - Through Gabay's equivalence (ADMM = DR on the dual, ρ = r), ADMM cannot converge strongly in general either.
  - Strong convergence of the (u^n, p^n) iterates is obtained when F₀ is **uniformly convex on bounded sets**. That is Glowinski's (iv), and in operator language it corresponds to uniform monotonicity in Bauschke–Combettes Thm 26.11. Uniform convexity gives p^n → p strongly. B injective with closed range, i.e. bounded below, then transfers this to u^n → u.
  - The multiplier converges only weakly, even then.
  - Closed range plus injectivity by itself only converts convergence of Bu^n into convergence of u^n of the same type.
- **Finite-dimensional sources.** Boyd et al. and Eckstein–Bertsekas Thm 8 are finite-dimensional, so their strong-convergence statements do not transfer.
- **The Lyapunov functional.** The user's V_n = (1/2r)‖η_{n+1} − η*‖² + (r/2)‖ν_n − ν*‖² with ρ = r is, up to a factor 1/2, Boyd et al. Appendix A's V^k = (1/ρ)‖y^k − y*‖² + ρ‖B(z^k − z*)‖² (MEM). There the splitting variable z plays the role of ν and B = −I in their constraint Ax + Bz = c. Bartels–Milicevic Thm 3.7 uses the same quantity, D_j² = ‖δ_λ‖² + τ²‖Bδ_u‖², which is the user's quantity times 2r when τ = r = ρ.
- **Opial.** Opial's lemma does **not** appear in Boyd, which is finite-dimensional. In Hilbert space it is the standard step: Fejér monotonicity of (η_n, ν_n) plus all weak cluster points being saddle points gives weak convergence (Bauschke–Combettes Lemma 2.47 / Thm 5.5; numbering MEM).
- **Glowinski's ρ ≠ r proof.** Glowinski's proof for ρ ≠ r uses a different functional, involving ‖λ^n − λ‖², rρ‖p^n − p‖² and a ‖Bu^n − p^n‖² term. This is MEM; Tao–Yuan cite Fortin–Glowinski eqs. (5.26)–(5.27) as the monotone decrease.

### (c) Which works prove or observe L² and H¹ (or W^{1,p}) rates, with exponents

| Work | Proved | Observed |
|---|---|---|
| Glowinski–Marrocco 1975 | W^{1,p}, p-dependent; p ≥ 2 commonly quoted as h^{1/(p−1)} (MEM); p < 2 exponent unverified | disk tests, errors vs p |
| Barrett–Liu 1993 | **order h** in W^{1,p} (p < 2) and in W^{1,q}, q ∈ [1,2] (p > 2), for regular u | numerics confirm |
| Liu–Barrett 1994 | quasi-norm optimal; L² bound ‖u−u_h‖_{L²} ≤ Ch^s only via Hölder (Thm 4.3) | none |
| Ebmeyer–Liu 2005 | optimal quasi-norm / W^{1,p} under weaker regularity | (not read) |
| Diening–Růžička 2007, Diening–Kreuzer 2008, BBDR 2012 | ‖F(∇u)−F(∇u_h)‖₂ ≤ c h ‖∇F(∇u)‖₂ for all 1 < p < ∞ and Orlicz φ (Δ₂ ∩ ∇₂) | AFEM numerics (not read) |
| Breit–Diening–Schwarzacher 2015 | F-quasi-norm h^α for p ∈ C^{0,α} | none |
| Del Pezzo–Martínez 2015 | W^{1,p(·)}: h^{p₁/2} (H² regularity); h under (1.4)–(1.5); p ≤ 2 | ≈ order 1 in W^{1,p(·)} |
| Del Pezzo–Martínez DCM (arXiv:1401.2442) | none | L² ≈ 0.8 (my estimate from their table; DG-P0) |
| Aragón–FB–Rubio 2023 | none | single mesh, no rates |

Main points:
- **No source found proves an O(h²) L² rate for the p-Laplacian.** An Aubin–Nitsche duality argument is not available in general for the degenerate or singular case. Our observed order 2 in L² and order 1 in ‖∇(u−u_h)‖_{L²} for smooth G and smooth solutions therefore goes beyond what is proved. It is consistent with linear-FE heuristics when g(t)/t is bounded above and below near the range of |∇u|.
- The reduced orders we see for large p and for a discontinuous coefficient are consistent with:
  - the regularity-limited quasi-norm theory (F(∇u) ∉ W^{1,2});
  - h^α for C^{0,α} exponents in Breit–Diening–Schwarzacher;
  - the p-dependent constants seen by Del Pezzo–Martínez.

## Doubtful items to check before citing

1. The numbering of Glowinski Ch. VI, Thms 4.1/5.1 and Remark 4.4. I inferred it from Del Pezzo–Martínez and Aragón–FB–Rubio; the book was not accessed.
2. Whether Glowinski's Thm 5.1 requires uniform convexity of F₀. Glowinski–Pan–Tai say yes; Del Pezzo–Martínez omit it.
3. The Glowinski–Marrocco rate exponents: the OCR was garbled.
4. The Boyd et al. Appendix A formulas (MEM).
5. Bauschke–Combettes theorem and lemma numbers other than 26.11(iii), which is confirmed via Bùi–Combettes.
6. Journal data marked MEM:
   - Gabay–Mercier, Lions–Mercier pages, Svaiter volume and pages, Bùi–Combettes pages, Bauschke–Moursi.
   - Breit–Diening–Schwarzacher DOI, Diening–Kreuzer DOI.
   - Diening–Fornasier–Tomasi–Wank volume, pages and fourth author.
   - Balci–Ortner–Storn volume and pages.
   - Liu–Yan, Hirn, Chen–Levine–Rao.
7. The Bartels–Milicevic and Moradifam–Nachman journal versions were not identified.
