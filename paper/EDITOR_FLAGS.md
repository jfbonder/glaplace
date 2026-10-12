# Editor flags (not changed in the text; for the authors / the mathematical referee)

No theorem, lemma, proposition, corollary, definition, assumption, displayed equation or proof was modified
(checked mechanically: all `theorem/lemma/proposition/corollary/maintheorem/definition/assumption/proof/equation/align`
environments are byte-identical to `new/main.tex`). Only notation in a *literature summary* was changed: the second
functional in the reproduction of Glowinski's theorem (Section 2.6) is now `\Gamma` instead of `G_2`, to avoid a clash
with the Young function `G`.

Points the authors should check before submission:

1. **Unverified sources.** Several literature statements still rest on sources the authors say they have not read in full:
   Glowinski (1984) itself (Section 2.6 relies on the reproduction in Glowinski-Pan-Tai), Cianchi-Maz'ya 2011 ("abstract
   only"), the boundary part of Lieberman's Thm 1.7, Fan 2007, the notion of solution in Cianchi-Maz'ya 2018, and whether
   De Filippis-Mingione Thms 1.4-1.5 cover e^t-1-t. A journal referee will likely object to citing a paper "from the abstract
   only" (Table 1, row 2). Either verify or drop that row/citation. The caveats were kept, only shortened.
2. **Table 2, "theory" column for p = 1.2, 1.1.** These rows were previously only in text/Table 5 without a theory label.
   I labelled them "B(b), L" by combining Remark 5.2 ((S) holds for the radial data) with Theorem B(b) and (3.1) for t^p/p.
   Please confirm.
3. **Abstract wording.** The abstract says strong convergence holds "under a monotonicity condition on g(x,t)/t"; the
   precise hypothesis (also C^1, g(x,0)=0 and a uniform positive lower bound) is in the introduction and Theorem B.
   Acceptable for an abstract in my view, but check you are happy with it.
4. **Figure 1 caption** referred to `convergence_order.csv`, which is not in `results/`. The reference was removed; make
   sure the raw data are in the released repository.
5. **Remark 4.3 (Del Pezzo-Martinez).** The claim that their injectivity hypothesis fails for their lifted gradient is a
   statement about another paper; it is kept (shortened) but should be double-checked.
6. **Fits that use unconverged runs** (e^t-1 two finest meshes; all p = 1.1, 1.2 runs) are marked in Table 2 and in the
   limitations paragraph; the p = 1.1 orders contain an iteration component of a few percent. Kept as is.
