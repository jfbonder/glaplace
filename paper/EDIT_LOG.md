# Edit log (editorial pass; working copy `paper_v3/edit`, original `paper_v3/new` untouched)

**Page count:** 29 -> 23 (body ends on p. 21; references about 2.3 pages). Same class, margins and font size.
Tables: 8 -> 2. Figures: 4 -> 2. Bibliography: 53 -> 52 entries (only `Marcellini1989` dropped, no longer cited);
every remaining entry is cited and every \cite key resolves. Compiles cleanly with latexmk + pdflatex (no errors, no
undefined references, no overfull boxes; one harmless underfull vbox).

Mathematics: all theorem-like environments, proofs and displayed equations are unchanged (checked mechanically).
See `EDITOR_FLAGS.md` for points the authors should check.

## Front matter
- Title shortened to "Convergence of Glowinski's decomposition--coordination method for the g-Laplacian" (the old subtitle
  listed the contents). **Authors decide**; the running head is unchanged.
- Abstract rewritten: 143 words, no itemised results, no numbers.
- `\tableofcontents` removed (unusual in a ~20-page journal article). **Authors decide.**

## Introduction
- Paragraphs on the method and on the classical theory merged and shortened; the p>2 subtlety is now only in Section 2.6.
- "Contributions" reduced from six long items to five short ones; the "what is new vs standard" paragraph kept.
- "What is not proved" kept with all five caveats (no rates; only rho = r; (S) is an assumption on the solution;
  u_h -> u not proved in general, two-phase excluded; d = 2), slightly reworded.

## Sections 2-5 (theory)
- Remarks 2.10 (related results) and 2.11 (weak, not strong, convergence) merged into one remark (both labels kept).
  The Bui-Combettes caveat and "we do not know whether strong convergence can fail for (P)" kept.
- Section 2.6: meta-commentary compressed into one sentence ("we have not been able to consult Glowinski 1984");
  the caveat that hypothesis (iv) must be read on dom F for p>2 kept. `G_2` -> `\Gamma` (notation clash with G).
- Section 3: the six-item counterexample list became one paragraph; literature statements shortened but every
  "not verified / conditionally" caveat kept; Table 1 kept (the "when does (S) hold" summary), row 2 shortened, caption no
  longer points to the internal file `TODO_VERIFY.md`.
- Section 4: Remarks 4.1 (discrete nonhomogeneous data) and 4.2 (continuous) merged into one remark with two paragraphs
  (both labels kept). Remark on related literature shortened.
- Section 5: the "Gap" bullet duplicated Section 3.3 and is now a cross-reference; "no gap" results and the
  consequences for conforming elements are one paragraph each. **Remark 5.2 and Proposition 5.3 kept**: Remark 5.2 is
  what makes (S) and the absence of a gap rigorous for every numerical test (including e^t-1, e^t-1-t and the
  indicator two-phase coefficient), and Proposition 5.3 is the only u_h -> u result; neither is redundant. Only the
  surrounding commentary was compressed; the caveats (exact discrete minimisers, no rate, two-phase not covered) kept.

## Section 6 (numerics)
- Kept: **Table 2** (new compact table of observed orders): fitted orders for all radial cases, plus p = 1.2, 1.1 and the
  manufactured solution, with iteration counts on the coarsest/finest mesh. Finest-mesh error values moved to one sentence.
- Dropped tables (data remain in `results/`):
  - nodal vs true error (`sc.txt`) -> one sentence in Protocol (underestimation by factors 7-25, spurious orders 2.1-2.5);
  - iteration counts per mesh -> column of Table 2 plus observation (5);
  - p close to 1 (`pto1.txt`, `pto1_hist.txt`) -> paragraph with residual ranges, iteration-error caveat (3% for p=1.1)
    and orders in Table 2;
  - manufactured solution (no separate file in `results/`; values copied from the old Table 6) -> pairwise orders and errors in one sentence, row in Table 2;
  - iterations vs r (`rman.txt`, `rdead.txt`) -> counts inline; the stopping-criterion caveat (r = 100) kept;
  - Lyapunov check (`lyapunov.json`) -> one paragraph with the extreme values.
- Dropped figures: dead-core profile (`exp_dead_core.png`; the flatness is stated quantitatively) and Lyapunov decay
  (`lyapunov.pdf`; replaced by "appears geometric"). Kept: convergence orders (Fig. 1) and two-phase sector (Fig. 2,
  scaled to 0.85\textwidth). **Authors decide** whether the dead-core figure should come back as a third figure.
- The remark about the original `lyapunov_check.py` script was removed (implementation detail).
- "Comparison with the literature on rates" and "Limitations" merged into one subsection (both labels kept); the item on
  iteration counts depending on h and r was removed as it repeats observation (5) and Remark 4.3.
- References in the text to `NOTES.md` and `results/exactcheck.txt` replaced by "scripts and raw data are distributed with
  the code (Software availability)".

## Conclusions
- The long recap (which repeated abstract and introduction) reduced to one sentence; the five open problems kept as an
  inline list.

## Unchanged
- TODO placeholders (addresses, e-mail of A. Aragon, funding).
- All citation keys that remain used; British spelling (minimise, discretisation) throughout the text.

## For the authors to decide
1. New title; removal of the table of contents.
2. Whether 23 pages is acceptable. Further cuts at the same level of care would have to touch proofs (excluded here) or
   drop a literature block (e.g. the Cianchi-Maz'ya 2011 row cited "from the abstract only", or the Ball-Mizel /
   Giaquinta-Hong-Marcellini examples, which concern u-dependent or anisotropic integrands).
3. Whether to restore the dead-core figure.
