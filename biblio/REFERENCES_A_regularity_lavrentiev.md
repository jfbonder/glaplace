# References A: gradient regularity (Assumption (S)), counterexamples, Lavrentiev phenomenon

Prepared 2026-10-12 for the paper on the convergence of Glowinski ADMM for the g-Laplacian.

**Download status: NOTHING WAS DOWNLOADED.** The session's egress proxy denies `arxiv.org` (organization policy: `connect_rejected`, HTTP 403 on CONNECT). The other open-access hosts tested (cvgmt.sns.it, air.unipr.it, flore.unifi.it, degruyterbrill.com) were unreachable from curl, and wrap.warwick.ac.uk returned 403. Per instructions, I did not try to route around these blocks. Every entry below gives a link so you can fetch it yourself. The texts were read through the WebFetch tool, which returns a summary produced by a model, not a verbatim copy. "Read in full text" therefore means the theorem statements were taken from the full paper through that tool. Check theorem numbers against the PDF before citing.

Verification labels:
- **FT**: theorem statements read from the full text (arXiv PDF/HTML or an open repository copy).
- **ABS**: abstract or bibliographic page only.
- **SEC**: known only from how other papers cite it.
- **MEM**: from memory, not checked against any source in this session. Treat as uncertain.

Notation: n (or d) is the dimension. H(x,t) = t^p + a(x) t^q, with 1 < p ≤ q and 0 ≤ a ∈ C^{0,α}.

---

## TOPIC A — Regularity

### A1. Lieberman 1991
```bibtex
@article{Lieberman1991,
  author  = {Lieberman, Gary M.},
  title   = {The natural generalization of the natural conditions of {L}adyzhenskaya and {U}ral'tseva for elliptic equations},
  journal = {Comm. Partial Differential Equations},
  volume  = {16}, number = {2-3}, pages = {311--361}, year = {1991},
  doi     = {10.1080/03605309108820761}}
```
- Link: https://doi.org/10.1080/03605309108820761 (paywalled). A related shorter version, Proc. Centre Math. Appl. ANU vol. 24 (interior results only), is at https://maths.anu.edu.au/files/CMAProcVol24-Lieberman.pdf
- File: NOT DOWNLOADED.
- Content:
  - **Structure condition (from the ANU proceedings, FT, OCR-garbled).** g = G' satisfies δ ≤ t g'(t)/g(t) ≤ g0 for all t > 0, with δ, g0 > 0. Here **g depends only on t**.
  - **x-dependence (MEM).** It enters only through the vector field A(x,z,p), with the ellipticity and growth of A_p controlled by g(|p|)/|p|. The Hölder condition is of the form |A(x,z,p) − A(y,w,p)| ≤ Λ(1+g(|p|))(|x−y|^θ + |z−w|^θ), and the lower-order term satisfies |B| ≤ Λ(1+g(|p|)|p|).
  - **Theorem 1.7 (SEC).** Other papers cite it for C^{1,α}(Ω̄) of weak solutions with bounded right-hand side; example: arXiv 1310.5903, Remark 2.2, says "by Lieberman [12, theorem 1.7] … u ∈ C^{1,α}(Ω̄)". Martínez–Wolanski (arXiv math/0602388) cite the interior C^{1,α} estimate from the proof of Theorem 1.7 (p. 346), and Theorem 1.2 for a sup bound.
  - **Hypotheses of Theorem 1.7 (MEM).** u bounded (|u| ≤ M0), ∂Ω ∈ C^{1,α}, Dirichlet data φ ∈ C^{1,α}. Conclusion: u ∈ C^{1,β}(Ω̄), with β and the norm depending on δ, g0, n, α, Λ, M0, Ω.
  - **Caveat.** Antonini (Nonlinear Anal. 2026, arXiv 2601.07140) credits Lieberman's Orlicz paper only with *interior* gradient Hölder continuity. Antonini attributes boundary gradient regularity to another Lieberman paper (probably the 1988 p-growth paper, Nonlinear Anal. 12) and calls his own global Orlicz result the natural generalization of that paper. The boundary statement of Lieberman 1991 Thm 1.7 must therefore be checked in the PDF.
- Status: structure condition FT (ANU version). Theorem 1.7 SEC/MEM. **Not verified in the Comm. PDE text.**

### A1'. Antonini 2026 (modern global C^{1,β} for Orlicz-type operators, verified)
```bibtex
@article{Antonini2026,
  author  = {Antonini, Carlo Alberto},
  title   = {Local and global {$C^{1,\beta}$}-regularity for uniformly elliptic quasilinear equations of $p$-{L}aplace and {O}rlicz-{L}aplace type},
  journal = {Nonlinear Anal.}, year = {2026}, doi = {10.1016/j.na.2026.114266},
  eprint  = {2601.07140}, archivePrefix = {arXiv}}
```
- arXiv: 2601.07140 (v3, 27 Aug 2026). File: NOT DOWNLOADED.
- Hypotheses:
  - −div A(x,∇u) = f, with A_ξ controlled by a(|ξ|) = b(|ξ|)/|ξ| and −1 < i_a ≤ s_a < ∞.
  - x-Hölder dependence |A(x,ξ) − A(y,ξ)| ≤ Λ(1+b(|ξ|))|x−y|^α, still with a **single** Young function b, not x-dependent.
- Results:
  - **Corollary 1.4 (Dirichlet).** Bounded C^{1,α} domain, f ∈ L^d with d > n, Dirichlet data in C^{1,α}. Then u ∈ C^{1,β}(Ω̄).
  - **Corollary 1.5.** The Neumann analogue.
  - **Theorems 1.2–1.3.** Local-at-the-boundary versions.
- Status: FT (HTML; the bibliography was truncated).

### A2. Cianchi–Maz'ya
```bibtex
@article{CianchiMazya2011,
  author  = {Cianchi, Andrea and Maz'ya, Vladimir G.},
  title   = {Global {L}ipschitz regularity for a class of quasilinear elliptic equations},
  journal = {Comm. Partial Differential Equations}, volume = {36}, number = {1}, pages = {100--133}, year = {2011},
  doi     = {10.1080/03605301003657843}}  % DOI from memory: verify
@article{CianchiMazya2014JEMS,
  author  = {Cianchi, Andrea and Maz'ya, Vladimir G.},
  title   = {Gradient regularity via rearrangements for $p$-{L}aplacian type elliptic boundary value problems},
  journal = {J. Eur. Math. Soc.}, volume = {16}, pages = {571--595}, year = {2014}}  % journal data from memory
@article{CianchiMazya2014ARMA,
  author  = {Cianchi, Andrea and Maz'ya, Vladimir G.},
  title   = {Global boundedness of the gradient for a class of nonlinear elliptic systems},
  journal = {Arch. Ration. Mech. Anal.}, volume = {212}, pages = {129--177}, year = {2014}}  % MEM
@article{CianchiMazya2018,
  author  = {Cianchi, Andrea and Maz'ya, Vladimir G.},
  title   = {Second-order two-sided estimates in nonlinear elliptic problems},
  journal = {Arch. Ration. Mech. Anal.}, volume = {229}, pages = {569--599}, year = {2018},
  doi     = {10.1007/s00205-018-1223-7}, eprint = {1703.07446}}
```
**(a) Comm. PDE 2011.**
- Links: record at https://flore.unifi.it/handle/2158/398534 (no open PDF found).
- Content: Lipschitz continuity of solutions of Dirichlet and Neumann problems for a class including the p-Laplacian, under minimal integrability of the data and of the boundary curvature, and also on arbitrary bounded convex domains (ABS).
- Hypotheses (MEM, consistent with the 2014 and 2018 papers):
  - −div(a(|∇u|)∇u) = f with −1 < i_a ≤ s_a < ∞, where i_a = inf t a'(t)/a(t) and s_a = sup t a'(t)/a(t).
  - f ∈ L^{n,1}(Ω), n ≥ 3.
  - Ω convex, or ∂Ω ∈ W²L^{n−1,1}.
- Conclusion: ∇u ∈ L^∞(Ω), with ‖∇u‖_∞ ≤ C b^{-1}(C‖f‖_{L^{n,1}}), where b(t) = t a(t).
- Status: ABS + MEM.

**(b) "Gradient regularity via rearrangements" (JEMS 2014).**
- Open copy: https://users.mai.liu.se/vlama82/pdf/cm_rearr-gradient-12072011.pdf. Read FT.
- Hypotheses: a(t) with t a'(t)/a(t) ≥ p−2 and t a(t) ≤ C(t^{p−1}+1), p ∈ [2,n), n ≥ 3. Dirichlet or Neumann data. ∂Ω ∈ W²L^{n−1,1} or Ω convex.
- Results:
  - **Theorems 1.1 and 1.3.** A rearrangement estimate of |∇u|* by f**.
  - **Theorem 3.5.** f ∈ L^{n,1} gives ‖∇u‖_∞ ≤ C‖f‖_{L^{n,1}}^{1/(p−1)}.
  - The paper shows by example that boundedness of the gradient can fail without these domain assumptions.
- Status: FT.

**(c) ARMA 2018.**
- arXiv: 1703.07446.
- Results:
  - **Theorem 2.1 (local).** f ∈ L²_loc gives a(|∇u|)∇u ∈ W^{1,2}_loc.
  - **Theorem 2.3 (convex Ω, no boundary regularity).** Dirichlet or Neumann data and f ∈ L² give a(|∇u|)∇u ∈ W^{1,2}(Ω), with the two-sided bound C1‖f‖_2 ≤ ‖a(|∇u|)∇u‖_{W^{1,2}} ≤ C2‖f‖_2.
  - **Theorem 2.4.** The same conclusion for Lipschitz domains with ∂Ω ∈ W²L^{n−1,∞} (n ≥ 3), or W²L^{1,∞}log L (n = 2), plus a smallness condition on the curvature.
  - **Remarks 2.6–2.7.** With f ∈ L^{n,1} and ∂Ω ∈ W²L^{n−1,1}, the gradient is globally bounded (citing CiMa1–2). If moreover inf_{[0,M]} a > 0, then u ∈ W^{2,2}(Ω), and on convex Ω regardless of ∂Ω.
- Status: FT.
- All of (a)–(c) have **autonomous** a(|ξ|). x-dependent G is not covered.

### A3. Double phase

**A3.1 Colombo–Mingione, Regularity for double phase variational problems**
```bibtex
@article{ColomboMingione2015a,
  author = {Colombo, Maria and Mingione, Giuseppe},
  title  = {Regularity for double phase variational problems},
  journal= {Arch. Ration. Mech. Anal.}, volume = {215}, number = {2}, pages = {443--496}, year = {2015},
  doi    = {10.1007/s00205-014-0785-2}}
```
- Open preprint: https://cvgmt.sns.it/media/doc/paper/2503/cmarma3.pdf. File: NOT DOWNLOADED.
- Hypotheses: a ∈ C^{0,α} (eq. 1.2) and the **strict** inequality q/p < 1 + α/n (eq. 1.3).
- Results:
  - **Theorem 1.1.** Higher integrability, reverse Hölder, local boundedness.
  - **Theorem 1.2.** Du ∈ C^{0,β}_loc for local minimisers of ∫(|Dw|^p + a|Dw|^q) and of more general G(x,Dw).
  - **Theorem 1.3.** Morrey decay, giving u ∈ C^{0,γ}_loc for all γ < 1.
  - **Theorem 1.4.** The vectorial case.
  - **Theorem 1.5.** Fractional differentiability.
- All results are **interior**.
- When q/p > 1+α/n the paper records that minimisers may fail to be W^{1,q}_loc, the Lavrentiev gap appears (its Theorem 4.1 relates the absence of the gap to (1.3)), and minimisers may be discontinuous with singular set of dimension close to n−p.
- Status: FT (preprint version).

**A3.2 Colombo–Mingione, Bounded minimisers**
```bibtex
@article{ColomboMingione2015b,
  author = {Colombo, Maria and Mingione, Giuseppe},
  title  = {Bounded minimisers of double phase variational integrals},
  journal= {Arch. Ration. Mech. Anal.}, volume = {218}, number = {1}, pages = {219--273}, year = {2015},
  doi    = {10.1007/s00205-015-0859-9}}
```
- Record: https://air.unipr.it/handle/11381/2814708 (no file). File: NOT DOWNLOADED.
- Result (ABS): **bounded** minimisers are C^{1,β} (locally) under the sharp condition **q ≤ p + α**.
- Status: ABS. BCM 2018 restates the hypothesis as (1.3): u ∈ L^∞, a ∈ C^{0,α}, q ≤ p+α.

**A3.3 Baroni–Colombo–Mingione 2018**
```bibtex
@article{BaroniColomboMingione2018,
  author = {Baroni, Paolo and Colombo, Maria and Mingione, Giuseppe},
  title  = {Regularity for general functionals with double phase},
  journal= {Calc. Var. Partial Differential Equations}, volume = {57}, number = {2}, pages = {Art. 62, 48 pp.}, year = {2018},
  doi    = {10.1007/s00526-018-1332-z}, eprint = {1708.09147}}  % journal data MEM
```
- arXiv: 1708.09147. File: NOT DOWNLOADED.
- Standing assumptions: **either** (1.2) a ∈ C^{0,α} and q/p ≤ 1+α/n, now **including the borderline equality, here for the first time**; **or** (1.3) u ∈ L^∞, a ∈ C^{0,α} and q ≤ p+α.
- Results:
  - **Theorem 1.1.** Du ∈ C^{0,β0}_loc for local minimisers of general non-autonomous double phase integrands.
  - **Theorem 1.2.** Under a continuous modulus: u ∈ C^{0,θ}_loc for all θ < 1, with Morrey decay.
  - **Theorem 1.3.** If u ∈ C^{0,γ}(Ω) and **q < p + α/(1−γ)**, then Du ∈ C^{0,β1}_loc.
  - **Theorem 1.4 (no Lavrentiev).** Under (1.3), or under u ∈ C^{0,γ} with q ≤ p + α/(1−γ), on every ball B ⋐ Ω there are u_k ∈ W^{1,∞}(B) with u_k → u in W^{1,p} and F(u_k,B) → F(u,B).
- The paper says the counterexamples of its refs. [26, 27] (ELM 2004, FMM 2004) show (1.2) and (1.3) are sharp. It also notes that (1.2) and (1.3) coincide for p = n.
- All results are **interior**.
- Status: FT.

**A3.4 De Filippis–Mingione 2021**
```bibtex
@article{DeFilippisMingione2021,
  author = {De Filippis, Cristiana and Mingione, Giuseppe},
  title  = {Lipschitz bounds and nonautonomous integrals},
  journal= {Arch. Ration. Mech. Anal.}, volume = {242}, pages = {973--1057}, year = {2021},
  doi    = {10.1007/s00205-021-01698-5}, eprint = {2007.07469}}  % volume/pages MEM
```
- arXiv: 2007.07469. File: NOT DOWNLOADED.
- Results (all for minimisers of ∫F(x,Dw) − f·w, with Du ∈ L^∞_loc):
  - **Theorem 1.1.** (p,q) growth with |∂_x F| controlled by h ∈ L^d, d > n, and q/p < 1 + min{1/n − 1/d, m_p}. For p ≥ 2 or f ≡ 0 this is q/p < 1+1/n−1/d.
  - **Theorem 1.3 (double phase).** 0 ≤ a ∈ W^{1,d}, d > n, **q/p ≤ 1 + 1/n − 1/d**, and also q/p < p if n = 2.
  - **Theorems 1.4–1.5.** Exponential and iterated-exponential integrands.
  - **Theorem 1.6.** Uniformly elliptic A(x,|Dw|) with |f|, |Dc| in a Lorentz-type space X.
- All results are **interior/local**.
- The paper describes q/p ≤ 1+α/n with a ∈ C^{0,α} as necessary (citing ELM, FMM) and sufficient (citing BCM).
- Status: FT.

**A3.5 Hästö–Ok (generalized Orlicz, maximal regularity)**
```bibtex
@article{HastoOk2022,
  author = {H{\"a}st{\"o}, Peter and Ok, Jihoon},
  title  = {Maximal regularity for local minimizers of non-autonomous functionals},
  journal= {J. Eur. Math. Soc.}, volume = {24}, number = {4}, pages = {1285--1334}, year = {2022},
  doi    = {10.4171/JEMS/1118}, eprint = {1902.00261}}  % vol/pages MEM
@article{HastoOk2022ARMA,
  author = {H{\"a}st{\"o}, Peter and Ok, Jihoon},
  title  = {Regularity theory for non-autonomous partial differential equations without {U}hlenbeck structure},
  journal= {Arch. Ration. Mech. Anal.}, volume = {245}, pages = {1401--1436}, year = {2022}}  % MEM
```
- arXiv: 1902.00261. File: NOT DOWNLOADED.
- Hypotheses (Theorem 1.1): φ(x,t) with ∂_tφ satisfying (A0), (Inc)_{p−1} and (Dec)_{q−1}; local minimiser of ∫φ(x,|∇u|).
- Results:
  - **(1)** (VA1), i.e. φ⁺_{B_r}(t) ≤ (1+ω(r))φ⁻_{B_r}(t) on the relevant range, gives u ∈ C^{α}_loc for every α < 1.
  - **(2)** (VA1) with ω(r) ≤ c r^β gives u ∈ C^{1,α}_loc.
  - **Theorems 7.2 and 7.4.** The same under the weaker (wVA1).
- Special cases:
  - Double phase: (VA1) holds iff q/p < 1+β/n. The borderline q/p = 1+β/n follows via Theorem 7.4 and (wVA1).
  - Variable exponent: ω_p(r) log(1/r) → 0 gives C^α for all α; ω_p(r) ≲ r^β gives C^{1,α}.
  - Perturbed Orlicz a(x)ψ(t).
- **Interior only.**
- Status: FT (main theorem); the corollary details are summarised by the paper itself.

**A3.6 Byun–Oh**
```bibtex
@article{ByunOh2020,
  author = {Byun, Sun-Sig and Oh, Jehan},
  title  = {Regularity results for generalized double phase functionals},
  journal= {Anal. PDE}, volume = {13}, number = {5}, pages = {1269--1300}, year = {2020},
  doi    = {10.2140/apde.2020.13.1269}}
@article{ByunOh2017,
  author = {Byun, Sun-Sig and Oh, Jehan},
  title  = {Global gradient estimates for non-uniformly elliptic equations},
  journal= {Calc. Var. Partial Differential Equations}, volume = {56}, pages = {Art. 46}, year = {2017}}  % MEM
```
- 2020 paper (ABS): optimal condition on the modulating coefficient for Hölder continuity and Harnack for quasiminimisers with (G,H)-growth. Interior.
- 2017 paper (MEM): global Calderón–Zygmund (L^γ) estimates for H(x,Du) in Reifenberg-flat domains. **This is not L^∞.**
- Status: 2020 ABS; 2017 MEM.

**A3.7 Boundary or global results for double phase / Musielak–Orlicz (searched)**
- **Missaoui, "Global C^{1,α}-regularity for Musielak–Orlicz equations in divergence form", arXiv 2601.11495 (2026, single-author preprint, not refereed as far as I can see).**
  - Equation: −div A(x,u,Du) = B with G(x,t) satisfying 1 < g⁻ ≤ t g/G ≤ g⁺ < n and a uniform log-type condition (G2). Bounded solutions (|u| ≤ M).
  - **Theorem 1.5 (Dirichlet, ∂Ω ∈ C^{1,β3}, data C^{1,β3}).** u ∈ C^{1,α}(Ω̄). **Theorem 1.6.** The Neumann analogue.
  - Example 1.1(7) lists t^p + a(x)t^q with a ∈ C^{0,α} and q ≤ p+α as satisfying the hypotheses, but the paper does not verify this in detail.
  - Status: FT through WebFetch. **Treat as unverified: needs careful checking before you rely on it.**
- **Bögelein–Duzaar–Marcellini–Scheven, "Boundary regularity for elliptic systems with p,q-growth", J. Math. Pures Appl. 159 (2022) 250–293 (MEM).** Global gradient bound for autonomous (p,q) systems with f ∈ L^d, d > n (SEC via De Filippis–Piccinini).
- **De Filippis–Piccinini, "Borderline global regularity for nonuniformly elliptic systems", arXiv 2206.15330 (IMRN 2023, MEM).**
  - **Theorem 1.1.** Autonomous Uhlenbeck-structure (p,q) systems, Ω convex, n ≥ 3, 0 ≤ q−p < 2p/(n−1) (as stated in the introduction), f ∈ L(n,1), zero Dirichlet data. Gives Du ∈ L^∞ near every boundary point, hence globally.
  - **Theorem 1.2.** W^{2,2} up to the boundary.
  - Not double phase (FT).
- I found **no** verified global (up-to-the-boundary) Lipschitz or C^{1,α} theorem for **x-dependent** double phase in a refereed source. A paper titled "Tachikawa, boundary regularity of double phase" could not be located.

### A4. Variable exponent p(x)
```bibtex
@article{CosciaMingione1999,
  author = {Coscia, Alessandra and Mingione, Giuseppe},
  title  = {H{\"o}lder continuity of the gradient of $p(x)$-harmonic mappings},
  journal= {C. R. Acad. Sci. Paris S{\'e}r. I Math.}, volume = {328}, number = {4}, pages = {363--368}, year = {1999}}
@article{AcerbiMingione2001,
  author = {Acerbi, Emilio and Mingione, Giuseppe},
  title  = {Regularity results for a class of functionals with non-standard growth},
  journal= {Arch. Ration. Mech. Anal.}, volume = {156}, pages = {121--140}, year = {2001}}
@article{AcerbiMingione2002CR,
  author = {Acerbi, Emilio and Mingione, Giuseppe},
  title  = {Regularity results for electrorheological fluids: the stationary case},
  journal= {C. R. Math. Acad. Sci. Paris}, volume = {334}, number = {9}, pages = {817--822}, year = {2002},
  doi    = {10.1016/S1631-073X(02)02337-3}}
@article{AcerbiMingione2002,
  author = {Acerbi, Emilio and Mingione, Giuseppe},
  title  = {Regularity results for stationary electro-rheological fluids},
  journal= {Arch. Ration. Mech. Anal.}, volume = {164}, pages = {213--259}, year = {2002},
  doi    = {10.1007/s00205-002-0208-7}}  % MEM
@article{Fan2007,
  author = {Fan, Xianling},
  title  = {Global {$C^{1,\alpha}$} regularity for variable exponent elliptic equations in divergence form},
  journal= {J. Differential Equations}, volume = {235}, number = {2}, pages = {397--417}, year = {2007},
  doi    = {10.1016/j.jde.2007.01.008}}  % DOI MEM
```
- **Acerbi–Mingione 2002 (CR note; open at numdam, https://www.numdam.org/item/CRMATH_2002__334_9_817_0).**
  - Hypotheses: stationary electrorheological **systems** (symmetric gradient, pressure), p Hölder continuous (|p(x)−p(x0)| ≤ L|x−x0|^α), 9/5 < γ1 ≤ p ≤ γ2 (n = 3).
  - Conclusion: **partial** regularity. Du is Hölder on an open Ω0 with |Ω∖Ω0| = 0. Local.
  - The ARMA 164 paper is the full version (MEM).
  - Status: FT for the note. **Not** the scalar p(x)-Laplacian result.
- **Coscia–Mingione 1999 (MEM).** Minimisers of ∫|Du|^{p(x)} with p ∈ C^{0,α} have C^{1,β}_loc gradients (Hölder Du). Interior.
- **Acerbi–Mingione 2001 (MEM).** Functionals with p(x)-growth. Log-Hölder p gives Hölder continuity of minimisers; Hölder p gives C^{1,β}_loc. Interior.
- **Fan 2007 (SEC, via arXiv 2010.15098 and many others).**
  - Hypotheses: Theorems 1.1–1.3, p ∈ C^{0,μ}(Ω̄), **bounded** weak solutions, ∂Ω ∈ C^{1,β}, Dirichlet data. Also Neumann.
  - Conclusion: u ∈ C^{1,α}(Ω̄).
  - Status: SEC. **Hypotheses not read in the original.**
- **Log-Hölder** (limsup ω_p(r) log(1/r) < ∞) is what gives density of smooth functions, i.e. no Lavrentiev gap (Zhikov). By itself it gives only Hölder continuity of u, not Lipschitz.

### A5. Generalized Orlicz framework
```bibtex
@book{HarjulehtoHasto2019,
  author = {Harjulehto, Petteri and H{\"a}st{\"o}, Peter},
  title  = {Orlicz Spaces and Generalized {O}rlicz Spaces},
  series = {Lecture Notes in Mathematics}, volume = {2236}, publisher = {Springer, Cham}, year = {2019},
  doi    = {10.1007/978-3-030-15100-3}}
```
- Content: (A0), (A1), (aInc), (aDec), and density of smooth functions under (A0)–(A2) (MEM).
- Interior regularity for this framework is in Hästö–Ok (A3.5).

---

## TOPIC B — Counterexamples and the Lavrentiev phenomenon

### B1. Esposito–Leonetti–Mingione 2004 (and 1999)
```bibtex
@article{ELM2004,
  author = {Esposito, Luca and Leonetti, Francesco and Mingione, Giuseppe},
  title  = {Sharp regularity for functionals with $(p,q)$ growth},
  journal= {J. Differential Equations}, volume = {204}, number = {1}, pages = {5--55}, year = {2004},
  doi    = {10.1016/j.jde.2003.11.007}}
@article{ELM1999,
  author = {Esposito, Luca and Leonetti, Francesco and Mingione, Giuseppe},
  title  = {Higher integrability for minimizers of integral functionals with $(p,q)$ growth},
  journal= {J. Differential Equations}, volume = {157}, pages = {414--438}, year = {1999}}  % MEM
```
- Status: ABS (ResearchGate metadata) plus SEC (BGS 2022, BDS 2020, Balci–Diening–Surnachev 2023). File: NOT DOWNLOADED.
- **Positive result.** Integrands Hölder in x (exponent α): regularity of minimisers when q/p < (n+α)/n = 1+α/n. The abstract also says this bound is sharp. BGS cite this paper for the absence of the gap when q ≤ p + pα/d.
- Here "regularity" means higher integrability (Du ∈ L^q_loc) and fractional differentiability, obtained by analysing the Lavrentiev phenomenon (per BCM 2018). It does **not** mean C^{1,β}, which is CM 2015.
- **Counterexample (SEC).** Ω = (−1,1)^d with p < d < d+α < q. A weight a ∈ C^{0,α} and boundary data exist for which a Lavrentiev gap occurs and the minimiser is not W^{1,q}_loc. This extends Zhikov's 2D construction (with α = 1, p < 2 < 3 < q) to dimension d and α ∈ (0,1].
- The paper also answers a question of Marcellini by constructing a minimiser with an isolated singularity (ABS).

### B2. Fonseca–Malý–Mingione 2004
```bibtex
@article{FMM2004,
  author = {Fonseca, Irene and Mal{\'y}, Jan and Mingione, Giuseppe},
  title  = {Scalar minimizers with fractal singular sets},
  journal= {Arch. Ration. Mech. Anal.}, volume = {172}, number = {2}, pages = {295--307}, year = {2004},
  doi    = {10.1007/s00205-003-0301-6}}  % DOI MEM
```
- Link: abstract at https://math.cmu.edu/users/nw0z/publications/03-CNA-004/004abs. File: NOT DOWNLOADED.
- Result: for 1 < p < N < N+α < q and any ε > 0, there is a ∈ C^{0,α} such that ∫(|Du|^p + a|Du|^q) has a local minimiser in W^{1,p} whose set of non-Lebesgue points is closed with Hausdorff dimension > N−p−ε. So the minimiser is discontinuous, and in particular Du ∉ L^∞.
- Status: ABS.

### B3. Classical (p,q) counterexamples: Giaquinta 1987, Marcellini 1989/1991, Hong 1992
```bibtex
@article{Giaquinta1987,
  author = {Giaquinta, Mariano}, title = {Growth conditions and regularity, a counterexample},
  journal= {Manuscripta Math.}, volume = {59}, number = {2}, pages = {245--248}, year = {1987}}
@article{Marcellini1989,
  author = {Marcellini, Paolo}, title = {Regularity of minimizers of integrals of the calculus of variations with nonstandard growth conditions},
  journal= {Arch. Rational Mech. Anal.}, volume = {105}, number = {3}, pages = {267--284}, year = {1989}}
@article{Marcellini1991,
  author = {Marcellini, Paolo}, title = {Regularity and existence of solutions of elliptic equations with $p,q$-growth conditions},
  journal= {J. Differential Equations}, volume = {90}, number = {1}, pages = {1--30}, year = {1991}}
@article{Hong1992,
  author = {Hong, Min-Chun}, title = {Some remarks on the minimizers of variational integrals with nonstandard growth conditions},
  journal= {Boll. Un. Mat. Ital. A (7)}, volume = {6}, number = {1}, pages = {91--101}, year = {1992}}
```
- Source: SEC, via Hirsch–Schäffner, arXiv 1911.12822 (CPAA/Comm. Contemp. Math.).
- Giaquinta: an **autonomous** convex anisotropic integrand with p = 2, q = 4 that has an **unbounded** minimiser in dimension n ≥ 6.
- Marcellini 1991, §6: local boundedness fails in general for q > (n−1)p/(n−1−p).
- These examples are anisotropic, of the form Σ|z_i|^{p_i}. They are **not** of the radial g-Laplacian form, which is uniformly elliptic in the Lieberman sense, so they do not affect the autonomous case of our problem.

### B4. Zhikov
```bibtex
@article{Zhikov1986,
  author = {Zhikov, Vasili{\u\i} V.}, title = {Averaging of functionals of the calculus of variations and elasticity theory},
  journal= {Izv. Akad. Nauk SSSR Ser. Mat.}, volume = {50}, number = {4}, pages = {675--710}, year = {1986},
  note   = {English transl.: Math. USSR-Izv. 29 (1987) 33--66}}
@article{Zhikov1995,
  author = {Zhikov, Vasili{\u\i} V.}, title = {On {L}avrentiev's phenomenon},
  journal= {Russian J. Math. Phys.}, volume = {3}, number = {2}, pages = {249--269}, year = {1995}}
@article{Zhikov1997,
  author = {Zhikov, Vasili{\u\i} V.}, title = {On some variational problems},
  journal= {Russian J. Math. Phys.}, volume = {5}, number = {1}, pages = {105--116}, year = {1997}}
```
- Citation data: MEM. Content: SEC through Balci–Diening–Surnachev 2023 (arXiv 2312.05604) and BDS 2020.
- **Zhikov 1995, Example 3.1:** checkerboard geometry, Ω = (−1,1)², double phase with Lipschitz a (α = 1) and **p < 2 < 3 < q**, giving a Lavrentiev gap.
- **Zhikov 1986 (MEM):** variable exponent equal to p1 < 2 in two opposite quadrants and p2 > 2 in the other two, in 2D, giving W ≠ H and a gap.
- **Zhikov (MEM):** log-Hölder continuity of p excludes the gap.

### B5. Balci–Diening–Surnachev 2020 (fractal examples)
```bibtex
@article{BDS2020,
  author = {Balci, Anna Kh. and Diening, Lars and Surnachev, Mikhail},
  title  = {New examples on {L}avrentiev gap using fractals},
  journal= {Calc. Var. Partial Differential Equations}, volume = {59}, number = {5}, pages = {Art. 180}, year = {2020},
  doi    = {10.1007/s00526-020-01818-1}, eprint = {1906.04639}}  % article number MEM
```
- arXiv: 1906.04639 (theorem numbers below are from the arXiv v2; BGS cite the published version as [2, Theorem 34], so numbering may shift by one or two). File: NOT DOWNLOADED.
- Construction: Ω = (−1,1)^d, d ≥ 2. Fractal Cantor-type contact sets replace Zhikov's point.
- **Theorem 32 (double phase).** For d ≥ 2, p > 1, α > 0 and **q > p + α·max{1, (p−1)/(d−1)}**, there is a ∈ C^{α} with H ≠ W, a Lavrentiev gap, and W-minimiser ≠ H-minimiser for suitable boundary data. Remark 33: this removes the dimensional threshold p < d < d+α < q of earlier examples.
- **Theorem 29.** A discontinuous exponent with values in [p⁻, p⁺].
- **Theorem 30.** A **uniformly continuous** exponent with modulus 1/log(e+1/t)^κ, κ < 1 (just below log-Hölder), with the same three conclusions.
- **Section 4.3.** Weighted p-energy.
- In the gap case the true (W-)minimiser cannot be approximated in energy by smooth or Lipschitz functions. In particular it is not in W^{1,∞}, and neither is anything that is an energy limit of conforming FE functions.
- Status: FT.

### B6. Bulíček–Gwiazda–Skrzeczkowski 2022
```bibtex
@article{BGS2022,
  author = {Bul{\'\i}{\v{c}}ek, Miroslav and Gwiazda, Piotr and Skrzeczkowski, Jakub},
  title  = {On a range of exponents for absence of {L}avrentiev phenomenon for double phase functionals},
  journal= {Arch. Ration. Mech. Anal.}, volume = {246}, pages = {209--240}, year = {2022},
  doi    = {10.1007/s00205-022-01816-x}, eprint = {2110.13945}}  % pages/DOI MEM
```
- arXiv: 2110.13945. File: NOT DOWNLOADED.
- **Theorem 2.3 (double phase).** Bounded Lipschitz Ω ⊂ R^d, a ∈ C^{α}(Ω̄), α ∈ (0,1], **q ≤ p + α·max(1, p/d)**, boundary data u0 ∈ W^{1,q}, **no boundedness assumed**.
  - Conclusion: the infima over u0+W_0^{1,p}, u0+W_0^{1,q} and u0+C_c^∞ coincide, and C_c^∞ is dense in W_0^{1,ψ}.
  - The theorem statement in v1 contains a typo with p and q swapped.
- Sharpness: sharp for p ≤ d, by BDS [2, Thm 34].
- Status: FT.

### B7. Borowski–Chlebicka–De Filippis–Miasojedow 2024
```bibtex
@article{BCDM2024,
  author = {Borowski, Micha{\l} and Chlebicka, Iwona and De Filippis, Filomena and Miasojedow, B{\l}a{\.z}ej},
  title  = {Absence and presence of {L}avrentiev's phenomenon for double phase functionals upon every choice of exponents},
  journal= {Calc. Var. Partial Differential Equations}, volume = {63}, pages = {Art. 35}, year = {2024},
  doi    = {10.1007/s00526-023-02640-1}, eprint = {2303.05877}}  % vol/art MEM
```
- Result: a new sharp scale of weights (beyond C^{0,α}) excluding the gap for every choice of exponents p, q > 1, with a counterexample at the edge of the scale.
- Status: ABS.

### B8. Ball–Mizel 1985 and the 1D origin
```bibtex
@article{BallMizel1985,
  author = {Ball, John M. and Mizel, Victor J.},
  title  = {One-dimensional variational problems whose minimizers do not satisfy the {E}uler--{L}agrange equation},
  journal= {Arch. Rational Mech. Anal.}, volume = {90}, number = {4}, pages = {325--388}, year = {1985}}
@article{Lavrentiev1926,
  author = {Lavrentiev, M.}, title = {Sur quelques probl{\`e}mes du calcul des variations},
  journal= {Ann. Mat. Pura Appl.}, volume = {4}, pages = {7--28}, year = {1927}}  % volume/pages uncertain
@article{Mania1934,
  author = {Mani{\`a}, B.}, title = {Sopra un esempio di {L}avrentieff},
  journal= {Boll. Un. Mat. Ital.}, volume = {13}, pages = {147--153}, year = {1934}}
```
- Content (MEM):
  - Integrand (x² − u³)² |u'|^{14} + ε|u'|² on [−1,1] with suitable boundary values.
  - Minimiser ≈ k|x|^{2/3}, singular at 0.
  - inf over W^{1,1} < inf over W^{1,∞}.
  - The minimiser does not satisfy the (integrated) Euler–Lagrange equation, and the integrand is strictly convex and superlinear in u'. SEC (arXiv 1503.06748) confirms the last point.
- Lavrentiev 1926/27 volume and pages are uncertain.
- Status: MEM/SEC.

### B9. Numerical consequences of the Lavrentiev gap
```bibtex
@article{BallKnowles1987,
  author = {Ball, John M. and Knowles, Greg}, title = {A numerical method for detecting singular minimizers},
  journal= {Numer. Math.}, volume = {51}, number = {2}, pages = {181--197}, year = {1987}}  % vol/pages MEM
@article{Li1995,
  author = {Li, Zhiping}, title = {A numerical method for computing singular minimizers},
  journal= {Numer. Math.}, volume = {71}, pages = {317--330}, year = {1995}}  % MEM
@article{CarstensenOrtner2010,
  author = {Carstensen, Carsten and Ortner, Christoph},
  title  = {Analysis of a class of penalty methods for computing singular minimizers},
  journal= {Comput. Methods Appl. Math.}, volume = {10}, number = {2}, pages = {137--163}, year = {2010},
  doi    = {10.2478/cmam-2010-0008}}
@article{Ortner2011,
  author = {Ortner, Christoph}, title = {Nonconforming finite-element discretization of convex variational problems},
  journal= {IMA J. Numer. Anal.}, volume = {31}, number = {3}, pages = {847--864}, year = {2011},
  doi    = {10.1093/imanum/drq004}}
@article{BalciOrtnerStorn2022,
  author = {Balci, Anna Kh. and Ortner, Christoph and Storn, Johannes},
  title  = {Crouzeix--{R}aviart finite element method for non-autonomous variational problems with {L}avrentiev gap},
  journal= {Numer. Math.}, volume = {151}, pages = {779--805}, year = {2022},
  doi    = {10.1007/s00211-022-01303-1}, eprint = {2106.06837}}  % vol/pages MEM
@article{BalciKaltenbach2023,
  author = {Balci, Anna Kh. and Kaltenbach, Alex},
  title  = {Error analysis for a {C}rouzeix--{R}aviart approximation of the variable exponent {D}irichlet problem},
  eprint = {2303.10687}, archivePrefix = {arXiv}, year = {2023}}  % journal: IMA J. Numer. Anal.? verify
```
- **Carstensen–Ortner 2010 (ABS; repository copy https://wrap.warwick.ac.uk/43755/, marked restricted).**
  - A **conforming FEM fails exactly when the Lavrentiev phenomenon is present**.
  - They propose and analyse a splitting/penalisation method that converges under mild conditions on the stored energy, with guidance on choosing the penalty parameter.
  - Status: ABS.
- **Ortner 2011 (ABS; author manuscript https://wrap.warwick.ac.uk/id/eprint/43773/1/WRAP_Ortner_8-lavrentiev_nc.pdf, curl got 403).** For **convex** problems, nonconforming Crouzeix–Raviart FEM converges to the correct (W-)minimiser, with convergence of the energy. Status: ABS.
- **Balci–Ortner–Storn 2022 (FT, arXiv 2106.06837).**
  - Conforming schemes satisfy W_h ⊂ H, so their discrete minima converge to **inf_H F** (citing CO10), given exact or accurate quadrature.
  - **Theorem 1.** For a convex Carathéodory φ(x,ξ) with p⁻–p⁺ growth and suitable quadrature (options A or B), CR-FEM satisfies F_h(u_h) → min_W F, with weak convergence of broken gradients (strong if strictly convex).
  - 2D experiments on (−1,1)²: Zhikov/BDS-type p(x) and double phase examples show a visible energy gap between P1 and CR. The conforming scheme is very sensitive to quadrature.
- **Ball–Knowles 1987 and Li 1995 (MEM).** Element-removal and truncation methods for detecting singular minimisers in 1D Ball–Mizel-type problems.
- **Balci–Kaltenbach (ABS).** CR for p(x)-Dirichlet. A medius estimate for uniformly continuous exponents; a priori estimates for Hölder exponents, optimal for Lipschitz exponents.

---

## What this implies for Assumption (S): ∇u ∈ L^∞(Ω) up to the boundary
Setting: homogeneous Dirichlet data and f ∈ L^∞ (or L^d, d > n).

| Setting | ∇u ∈ L^∞ up to ∂Ω? | Reference (status) |
|---|---|---|
| Autonomous G, δ ≤ t g'(t)/g(t) ≤ g0 (i.e. −1 < i_a ≤ s_a < ∞), ∂Ω ∈ C^{1,α}, f ∈ L^∞ | **yes** (even C^{1,β}(Ω̄)) | Lieberman 1991 Thm 1.7 (SEC/MEM; boundary part to verify); Antonini 2026 Cor. 1.4 (FT, f ∈ L^d, d > n) |
| Autonomous G, Ω **convex** (no boundary smoothness) or ∂Ω ∈ W²L^{n−1,1}, f ∈ L^{n,1} | **yes** (Lipschitz) | Cianchi–Maz'ya Comm. PDE 2011 (ABS+MEM); JEMS 2014 Thm 3.5 (FT, p-type a, n ≥ 3) |
| Same, flux regularity: a(\|∇u\|)∇u ∈ W^{1,2}(Ω), f ∈ L², convex Ω | (second-order result, not L^∞) | Cianchi–Maz'ya ARMA 2018 Thm 2.3 (FT) |
| G(x,t) = a(x)ψ(t)-type x-perturbations bounded by (1+g(\|ξ\|))\|x−y\|^α | **yes** (C^{1,β}(Ω̄)) | Lieberman 1991 / Antonini 2026 (structure allows only this x-dependence) |
| p(x), p ∈ C^{0,μ}(Ω̄), ∂Ω ∈ C^{1,β}, bounded solution | **yes** (C^{1,α}(Ω̄)) | Fan 2007 (SEC) |
| p(x) log-Hölder only | **unknown / not in general**: u Hölder, no Lipschitz theorem | Acerbi–Mingione 2001 (MEM) |
| p(x) merely continuous (not log-Hölder) | **no** in general: Lavrentiev gap possible | BDS 2020 Thm 30 (FT) |
| Double phase, a ∈ C^{0,α}, q/p ≤ 1+α/n | **only interior** (C^{1,β}_loc) | CM 2015a (FT, strict); BCM 2018 Thm 1.1 (FT, ≤) |
| Double phase, bounded u (automatic for f ∈ L^∞), q ≤ p+α | **only interior** (C^{1,β}_loc); global claimed in an unrefereed preprint | CM 2015b (ABS); BCM 2018 (FT); Missaoui arXiv 2601.11495 (preprint, unverified) |
| Double phase, u ∈ C^{0,γ}, q < p+α/(1−γ) | **only interior** | BCM 2018 Thm 1.3 (FT) |
| Double phase, a ∈ W^{1,d}, d > n, q/p ≤ 1+1/n−1/d | **only interior** (Du ∈ L^∞_loc) | De Filippis–Mingione 2021 Thm 1.3 (FT) |
| General Musielak–Orlicz φ(x,t) with (VA1)/(wVA1) | **only interior** (C^{1,α}_loc) | Hästö–Ok JEMS 2022 (FT) |
| Autonomous non-uniformly elliptic (p,q), convex Ω, f ∈ L(n,1) | **yes** (systems, Uhlenbeck structure) | De Filippis–Piccinini (FT); Bögelein–Duzaar–Marcellini–Scheven 2022 (MEM) |
| Double phase with q > p + α·max(1,(p−1)/(n−1)), or p < n < n+α < q | **no**: gap, discontinuous minimisers possible | BDS 2020; ELM 2004; FMM 2004 |

Practical conclusion: (S) is rigorously available for autonomous Lieberman-type G (C^{1,α} or convex domains) and for Hölder p(x) (Fan 2007). For double phase and general Musielak–Orlicz G, only **interior** Lipschitz/C^{1,β} is established in refereed work. Boundary versions should be stated as an assumption, or cited to the 2026 preprint with caution.

## Counterexamples (summary)

| Example | Condition | Type | Minimiser in H¹_0 / bounded gradient? |
|---|---|---|---|
| Ball–Mizel 1985 | 1D, polynomial strictly convex superlinear integrand | singular minimiser ~ \|x\|^{2/3}; inf_{W^{1,1}} < inf_{W^{1,∞}} | Du unbounded; not Lipschitz |
| Zhikov 1986 | p(x) = p1 < 2 < p2 on checkerboard quadrants, d = 2 | H ≠ W, gap | W-minimiser not approximable by smooth functions |
| Zhikov 1995 Ex. 3.1 | double phase, a Lipschitz, p < 2 < 3 < q, d = 2 | checkerboard, gap | as above |
| ELM 2004 | double phase, a ∈ C^{0,α}, p < d < d+α < q | gap; minimiser ∉ W^{1,q}_loc | Du ∉ L^q_loc, so not Lipschitz |
| FMM 2004 | same range p < N < N+α < q | W^{1,p} minimiser with closed singular set of dim > N−p−ε | discontinuous |
| Giaquinta 1987 | autonomous anisotropic, p = 2, q = 4, n ≥ 6 | unbounded minimiser | not bounded (anisotropic, not of g-Laplacian type) |
| Marcellini 1989/91 | q > (n−1)p/(n−1−p), anisotropic | unbounded minimiser | not bounded |
| BDS 2020 Thm 32 | double phase, d ≥ 2, q > p + α·max{1,(p−1)/(d−1)} | fractal (Cantor) contact sets, gap, h_W ≠ h_H | W-minimiser not in H; not Lipschitz |
| BDS 2020 Thm 30 | p(x) uniformly continuous, modulus 1/log(e+1/t)^κ, κ < 1 | gap | as above |
| BGS 2022 (positive) | q ≤ p + α·max(1,p/d), no boundedness | no gap | (sharp for p ≤ d) |
| BCM 2018 Thm 1.4 (positive) | bounded u and q ≤ p+α; or u ∈ C^{0,γ} and q ≤ p+α/(1−γ) | no gap (local) | — |

Numerical consequence: whenever a gap is present, conforming P1 FE energies converge to inf_H, not to min_W (Carstensen–Ortner 2010; Balci–Ortner–Storn 2022). An ADMM/Glowinski scheme on a conforming P1 space can at best approximate the H-minimiser. Assumption (S) excludes this situation.
