# Merge log: MATH_REFEREE_REPORT.md -> final/main.tex

Base: `edit/main.tex` (copied to `final/`; `edit/` and `new/` untouched). Each fix applied was checked against the PDFs
in `biblio/`: Glowinski 1984 Ch. VI pp. 166-194, Eckstein-Bertsekas 1992, Bui-Combettes, Lieberman 1991, Cianchi-Maz'ya 2018,
De Filippis-Mingione, Carstensen-Ortner 2010, Balci-Ortner-Storn, and GPT16 (cam16-10).
Result: **23 pages**, pdflatex run twice, no errors, no undefined references, no overfull boxes. Abstract: 130 words.

## Applied
- **M1 (scope of Glowinski's Thm 5.1).** Section 2.6 is rewritten from the book: (4.1)-(4.6), Thm 5.1 conclusions, §7 (Opial), Rem. 5.1 (rho<2r, Gabay-Mercier), §6.3 (s-Laplacian only after discretisation), the same step order. The sentence "not able to consult" and the "weaker sense for p>2" paragraph are deleted.
  New **Lemma 2.13** (splittings (a) F0 = eps/2 ||.||^2 and (b) h-truncation). Proof re-checked. One small addition: delta_M must be continuous, so k_* is replaced by its continuous minorant (1/tau) int_0^tau k_*.
  The coverage list was checked: g' >= 3/2 for t^1.5/1.5+t^3/3 and g' >= 2 for t^3/log(1+t) (numerically, on [1e-8, 1e8]).
  Table 2 "theory" column now reads "G" for every row except t^p/p with p<2; those rows keep B(b).
  The text after Theorem B says that (a) is Glowinski's whenever ess sup k(x,1) < infinity.
  The "covers e^t-1" claim is removed. AFR23 (p(x) >= 2) is now justified.
- **M2 (attribution).** The new "What is classical / What is new" paragraphs replace the contributions list. They credit Glowinski (5.25), §7 and §4.2, and E-B Thm 8.
  Remark 2.10 now credits Glowinski (5.25) first. Remark 4.3 adds that C(ii) is classical and that the book's Thm 5.2 needs a strictly convex C^1 part.
  - *Deviation:* the report lists as new only B(b), Prop 2.16, Rem 5.2 and Prop 5.3. Following the task, the text also counts as new:
    - Theorem A under (S) alone (without (iii)-(iv));
    - the B(a) cases with ess sup k(x,1) = infinity;
    - Lemma 2.13 and Prop 2.14;
    - Prop 2.11;
    - the pointwise solver.

    Each of these is stated precisely.
- **M3 (Lieberman).** Thm 1.7 is now described as interior; the boundary part is attributed to the p. 320 remark and Lieberman 1988 (added to the bibliography). The text mentions u-dependence and |u| <= M0. Table 1 is updated. Confirmed on p. 320 of the PDF.
- **M4.**
  - New **Proposition 2.11**: (S) holds iff (eta_n + r nu_{n-1}) is bounded. This uses E-B Thm 7 (general Hilbert space) and the identification (Y1)-(Y3) of Thm 8, with zer(A + dF*) nonempty iff (S); checked.
  - Bui-Combettes in Remark 2.10 now says only that weak convergence of eta_n is sharp under (A1)-(A3), and nothing about nu_n or u_n. The identification of their B as a subdifferential was checked.
  - "What is not proved" (i), the end of §2.7 and open problem (1) are updated, citing GPT16 Rem. 23.
- **m1** The proof of Prop 2.14 now uses Glowinski's (4.6) and the integral formula. The statement says "no splitting as in (iv)".
- **m2** The DR-dual identification is cited to Gabay / E-B §5, not to GPT16 §3.3, which covers only the primal case B = I.
- **m3** Remark 2.9 is rewritten with the identity E_n + m_n + (r/2)||Bu_n - nu_{n-1}||^2 = E_{n-1}.
- **m4** Cianchi-Maz'ya 2018: "conditionally" and "to be checked" are removed, and a in C^1 is added. The text explains generalized = weak solution for d=2. Theorem B(b) is now said to be unconditional on convex plane domains.
- **m5** CM2011 is stated for d >= 3, with the source "abstract only" removed.
- **m6**
  - De Filippis-Mingione: added the f-condition and the fact that Thms 1.4-1.5 need exponent p > 1, so the exponentials are not covered.
  - Hasto-Ok: added omega(r) <= c r^beta, Delta_2 and no right-hand side.
  - Colombo-Mingione / BCM: added f = 0.
- **m7** BDS: "with a suitable linear term".
- **m8** Carstensen-Ortner: hypotheses stated. BOS: the H-minimum over W ∩ W^{1,inf}, the Thm 1 hypotheses and the §4.2 quadrature sensitivity are stated.
- **m9** CO10 Lemma 2.1 is cited in the proof of Prop 5.3.
- **m10** "l_h != 0" is added.
- **m11** "bounds the ADMM dual residual" is added.
- **m12** Del Pezzo-Martinez: their (H2) omits uniform convexity (harmless), and p(x) is in [p1, 2].
- **m13** The typo "A gap it is" no longer existed in edit/. The CM2011 parenthesis is removed.
- **TODO_VERIFY.md** is updated: items are marked RESOLVED / PARTLY RESOLVED, and section G covers new items 55-58.
  - Two items were marked W in the report: #1 (p > 2) and #15 (Lieberman).
  - The report says three were WRONG but marks only these two. I took the third to be #2, the GPT16 §3.3 "DR for the dual" citation (m2). Please confirm.

## Not applied / doubtful
- The optional m4 remark (S for t^p/p, d >= 3, p >= 2 - 2/d) was not added, for lack of space. It is not wrong.
- Gabay-Mercier (rho < 2r) is cited only through Glowinski's Rem. 5.1; no bibliography entry was added.
- The issue number of Lieberman 1988 was not printed (it is not verified).
- No report fix was found to be wrong. The only refinement is the delta_M continuity point in Lemma 2.13(b) above.

## Space (to stay at 23 pages)
Several passages were condensed without changing their content or caveats:
- the contributions list;
- the Lyapunov-check subsection (§6.6), which keeps its key numbers;
- the counterexamples paragraph;
- the rates/limitations section;
- the observations, p -> 1 and r sections;
- Remarks 4.2, 4.3 and 5.2;
- the conclusions.

Figures were scaled: Fig. 1 to 0.9\textwidth and Fig. 2 to 0.7\textwidth. Theorems and proofs were not altered, except:
- Prop 2.14 (statement and proof, m1);
- Remark 2.9;
- the new Lemma 2.13 and Prop 2.11.

## For the authors
- Check Lemma 2.13 and Prop 2.11 (new mathematics).
- Check the third "WRONG" item.
- Check the bibliographic data still marked MEM in TODO_VERIFY F.
- The TODO placeholders (addresses, funding) are unchanged.
