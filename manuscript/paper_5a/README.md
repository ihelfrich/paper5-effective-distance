# Paper 5a — The Untaken Limit

**A Pole at θ = -2 in the Head-Mayer Closed Form for Intra-National Distance**

Sole-author SSRN preprint by Ian T. Helfrich. Sibling of the Gonchar-Helfrich
*Effective Distance* panel paper in the parent directory.

## The contribution in one paragraph

The Head-Mayer (2002) closed-form expression for intra-national distance,
`d_ii = (2/(θ+2))^(1/θ) R`, has a pole at θ = -2 and is complex-valued for
θ < -2. The underlying integral `∫_0^R r^(θ+1) dr` diverges in 2D for θ ≤ -2.
The structurally consistent CES integration parameter θ = δ(1-σ) for empirically
defensible σ ∈ [4, 8] falls entirely in [-7, -3], below the pole. CEPII evaluates
the formula at θ = -1, one unit above the pole, and justifies the choice as
"the usual coefficient estimated from gravity models" — i.e., the reduced-form
slope, not the structural integration parameter. The applied gravity literature
has built on this conflation since 2005. The paper documents the structural
inconsistency, proves the 2D divergence, computes the regularization-dependent
implied d_ii, and estimates the consequence for the home-bias coefficient in
gravity regressions. A material share of the conventional border-puzzle
multiplier is a measurement artifact rather than a structural barrier.

## How to compile

```bash
export PATH="/Library/TeX/texbin:$PATH"
xelatex -interaction=nonstopmode paper_5a.tex
bibtex paper_5a
xelatex -interaction=nonstopmode paper_5a.tex
xelatex -interaction=nonstopmode paper_5a.tex
```

Output: `paper_5a.pdf` (target: ~25 pages when empirics drop in).

## Status

| Section | State |
|---|---|
| 1. Introduction | Drafted (with HM 2014 line 1287 collision quote; HM 2002 line 502 conflation quote) |
| 2. Closed-form and domain of validity | Drafted (with numerical table at the pole) |
| 3. 2D integral divergence theorem | Drafted (with proof and dimension remark) |
| 4. Regularization sweep (nb19c) | Drafted, awaiting figure insert |
| 5. Gravity-regression empirical evidence | Awaiting nb18d output |
| 6. Implications | Drafted |
| 7. Conclusion | Drafted |
| Appendix A. Notation | Drafted |
| Appendix B. Numerical verification | Drafted |

## Files

- `paper_5a.tex` — main manuscript
- `paper_5a_refs.bib` — minimal bibliography (full bibliography in parent `refs.bib`)
- `paper_5a.pdf` — compiled output (8 pages currently; ~25 when figures and tables drop)

## Novelty check

The contribution has been verified novel against:

- Head & Mayer (2002), CEPII WP 2002-01 — *they performed the conflation themselves*
- Mayer & Zignago (2011), CEPII WP 2011-25 — *they institutionalized the conflation*
- Head & Mayer (2014), Handbook of International Economics Vol. 4 — *they retained the parameter-name collision deliberately (line 1287: "we retain the symbol to emphasize the similarity in resulting terms") and flagged the functional-form question as future research (line 4533) without pursuing it*
- Rauch (2016), RIE 24(5):1167-1177 — *independently derives θ = -1 from a gravity-in-physics analogy; does not push beyond*
- Coughlin & Novy (2021), IER — *finds spatial-aggregation effects on border estimates; different mechanism*
- Eaton-Kortum, Caliendo-Parro, Costinot-Donaldson literature — *uses a separate, oppositely-signed θ parameter; the wrong-θ collision in HM 2014 is its own contribution*

Full novelty assessment: `/Volumes/HELFRICH-GD/Paper5_EffectiveDistance_Outputs/literature_intel/novelty_assessment_2026-05-21.md`.

## Next moves

1. nb18d completes → drop empirical numbers into Section 5
2. Generate figures: regularization sweep curve (Section 4) and home-bias-vs-d_min curve (Section 5)
3. SSRN posting target: within 2 weeks
4. Scoop-protection email to Thierry Mayer or Yoto Yotov once preprint PDF is in defensible shape

---

*Last updated: 2026-05-21*
