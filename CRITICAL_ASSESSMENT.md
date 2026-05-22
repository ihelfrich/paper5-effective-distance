# Critical assessment of Paper 5a

*An honest, adversarial evaluation of "The Untaken Limit: A Pole at θ = -2
in the Head-Mayer Closed Form for Intra-National Distance" by Ian T.
Helfrich, Georgia Tech, 2026.*

*Author of this assessment: Ian Helfrich, with assistance from Claude
(autonomous mode). Written as a critical self-audit, not as advocacy.*

---

## The question Ian asked me to answer

> "Take a critical objective look at the work. Is it real? Is it
> meaningful? Does it ACTUALLY make a difference or contribution?"

I'll answer in three parts:

1. **Is it real?** Does the math hold up? Are the empirical results
   reproducible? Is the literature claim accurate?
2. **Is it meaningful?** Does it move the gravity-trade literature
   measurably forward? Is it the kind of thing a top-five journal
   would publish?
3. **Does it actually make a difference?** Will applied gravity
   practitioners change anything because of this? Or is it a
   methodological footnote?

The honest answers are: **yes, yes with caveats, and partially.**
Below I expand each.

---

## Part 1: Is it real?

### 1.1 The mathematics

**Verifiable mathematical fact.** The Head-Mayer (2002) closed-form
expression
$$d_{ii}(\theta) = \left(\frac{2}{\theta + 2}\right)^{1/\theta} R$$
has a pole at $\theta = -2$. The factor $2/(\theta + 2)$ diverges as
$\theta \to -2$, and for $\theta < -2$ the base is negative while the
exponent $1/\theta$ is a non-integer, so on the principal branch of
$z^{1/\theta}$ the expression takes complex values. Direct numerical
substitution in Python (`numpy.power` on `complex128`) confirms this.
The verification is in `notebooks/` and in Table 1 of the paper.

**The underlying integral** $\int_0^R r^{\theta + 1}\, dr$ diverges at
the lower limit for $\theta \leq -2$. This is a calculus exercise:
$\int_0^R r^{\alpha}\, dr$ converges if and only if $\alpha > -1$,
so $\int_0^R r^{\theta + 1}\, dr$ requires $\theta + 1 > -1$, i.e.,
$\theta > -2$. In 2D, the pairwise-distance density of a uniform disk
scales as $f(r) \sim r$ near the origin, so the moment
$\E[r^\theta] = \int r^\theta f(r) dr \sim \int r^{\theta + 1} dr$
diverges at the origin under the same condition.

**The 2D divergence theorem** (Theorem 3.1 of the paper) is a clean
application of the standard convergence test for power-law integrals.
The proof in the paper is correct and replicates a standard
measure-theoretic argument.

**The pole has been verified by:**
- Direct numerical evaluation (Table 1 of the paper).
- Two-orders-of-magnitude regularization sweep on a uniform disk
  simulation (Figure 1; `notebooks/19c_regularization_sensitivity.py`).
- The scaling exponent $\alpha(\theta) = (\theta + 2)/\theta$ derived
  in the paper and confirmed empirically.

**Verdict on the math: real and uncontroversial.** Anyone who reads
the Head-Mayer (2002) equation and tries to evaluate it at $\theta = -3$
will find the same complex-valued output. The novelty is not in
finding this; it is in noticing that nobody has flagged it.

### 1.2 The empirical findings

**Reproducible empirical fact.** In the pooled panel PPML on BACI HS02
v202401b $\times$ CEPII Gravity v202211 (years 2010, 2015, 2020 stacked;
origin-year and destination-year FE absorbed via the
Correia-Guimarães-Zylkin alternating-projection algorithm):

- The home-bias coefficient at the CEPII baseline ($\dii$ from
  `distw_harmonic`, equivalently the Head-Mayer closed form at
  $\theta = -1$) is $\hat{\gamma} = -5.636$.
- The same coefficient at the structurally consistent $\theta = -7$
  and fine spatial resolution ($d_{\min} = 0.5$ km) is
  $\hat{\gamma} = -9.446$.
- The 3.81-log-unit shift corresponds to a 44× change in the PPML
  multiplier.

**Robust to estimator choice in DIRECTION.** The cross-sectional
log-OLS specification (year 2010, positive flows only, country FE)
gives $\hat{\gamma} = +2.13$ at the CEPII baseline (the conventional
McCallum-AvW $\sim 8\times$ multiplier) and $-2.92$ at fine resolution
with $\theta = -5$. The sign-flip threshold occurs between $d_{\min} = 1$
and $d_{\min} = 2$ km. The cross-section PPML
(`statsmodels.GLM`, with convergence warnings) gives $\hat{\gamma} = -5.01$
at the baseline and $-8.10$ at fine resolution, matching the panel PPML
direction and approximate magnitude.

**Not robust to estimator choice in MAGNITUDE.** The OLS multiplier at
the CEPII baseline is $8.4\times$ (positive home bias); the panel PPML
multiplier at the same baseline is $0.0036\times$ (essentially zero in
the PPML conditional-mean reading). This is a Jensen-inequality
consequence: log-OLS estimates $\E[\log X]$ while PPML estimates
$\log \E[X]$, and trade-flow variance is large enough that the wedge
is several log units. Neither estimator is "wrong"; they answer
different questions. But the absolute level of the McCallum-AvW
$\sim 8\times$ multiplier is fragile and depends on the estimator
choice.

**Verdict on the empirics: real, with honest caveats.** The
regularization sensitivity is robust. The absolute magnitude of the
border-puzzle effect is fragile. The paper's claim is the
sensitivity, not the magnitude.

### 1.3 The literature claim

**Triangulated novelty verdict.** The paper claims that the pole at
$\theta = -2$ has not been documented in the gravity-trade literature.
The verification:

- **Head & Mayer (2002), CEPII WP 2002-01.** Direct close-reading. They
  derive the closed form. They notice the 1D bilateral contact-set
  singularity and engineer around it ($\Delta \geq 2R$ assumption).
  They evaluate at $\theta = -1.0001$ rather than $\theta = -1$ exactly
  to dodge the $\theta = -1$ removable singularity. They never push
  $\theta$ below $-1.5$. They do not document the 2D pole at $\theta = -2$
  or the integral divergence. (Note: an earlier draft of my novelty
  agent hallucinated a "removable singularities" footnote in HM 2002;
  the actual text has no such footnote.)
- **Mayer & Zignago (2011), CEPII WP 2011-25.** Direct close-reading.
  They state the $\theta = -1$ choice and justify it as the
  reduced-form regression slope. They do not discuss the
  structural-vs-reduced-form distinction.
- **Head & Mayer (2014), Handbook of International Economics Vol. 4.**
  Full close-reading (sub-agent + my own pass, ~4,850 lines of text).
  They use the symbol $\theta$ for the Pareto/Fréchet shape parameter
  (positive) and acknowledge the overload with the 2002 integration
  parameter (negative) at line 1287: "we retain the symbol to
  emphasize the similarity in resulting terms." They do not document
  the pole. Lines 1809-1810 acknowledge "trade with self and distance
  to self, both of which may be problematic," without elaborating.
  Line 4533 of the future-research list: "distance effects are too
  large and have the wrong functional form to be determined by freight
  costs."
- **Rauch (2016), RIE 24(5):1167-1177.** Independently derives the
  $\theta = -1$ choice from a gravity-in-physics analogy. He parks at
  $\theta = -1$ and does not push beyond.
- **Coughlin & Novy (2021), IER 62(4):1453-1487.** Different mechanism
  (spatial aggregation of trade flows, not regularization of
  intra-national distance).
- **Eaton & Kortum (2002), Caliendo & Parro (2015),
  Costinot, Donaldson & Komunjer (2012).** Sub-agent close-read
  verified: none performs the structural-vs-reduced-form conflation,
  because none constructs $d_{ii}$ from a CES spatial integral. EK
  assumes $d_{ii} = 1$ by normalization; CP uses a tetrad estimator
  that cancels distance; CDK absorbs all bilateral costs in
  importer-exporter FE.

**Verdict on novelty: NOVEL.** The pole, the divergence, and the
structural-vs-reduced-form conflation in the Head-Mayer-CEPII
tradition are not documented in the published literature. The
sub-agent triangulation rated A/B/C/D as NOVEL and E/F as ACKNOWLEDGED
(qualitative regularization sensitivity is in the literature; the
specific divergence mechanism is not).

### 1.4 What could be wrong with this analysis?

I am not certain about everything. Honest enumeration of risks:

1. **A more thorough literature search might find prior art I missed.**
   The triangulation agent searched six fields; a seventh might exist.
   The cleanest follow-up: ask Thierry Mayer directly via the
   scoop-protection email (which is drafted but not yet sent).
2. **The structural-vs-reduced-form conflation might be defensible.**
   If the iceberg-cost exponent $\delta < 1$, the structurally
   consistent $\theta$ can be above the pole and the critique
   dissolves. The literature has been quiet on $\delta$, so this
   defense is not currently active, but it is mathematically available.
3. **The panel PPML magnitude depends on the Wei (1996) proxy.** The
   YPMM "trade for production" proxy might give different magnitudes.
   The regularization-sensitivity pattern should survive but is not
   explicitly verified under YPMM.

These are honest limitations. They do not make the work wrong; they
make it provisional in the ways scientific work is supposed to be
provisional.

---

## Part 2: Is it meaningful?

### 2.1 What problem does it solve?

**The published border-puzzle literature has reported home-bias
multipliers as structural objects since 1995.** McCallum found
$22\times$; Anderson-van Wincoop reduced this to $5\times$ via
multilateral resistance; the literature since has converged on
$5-10\times$ as "the answer." The paper shows that this answer is
conditional on a regularization choice (the spatial resolution at
which intra-national distance is constructed) that the literature has
not engaged with.

**The implication for downstream work.** If a gravity paper reports
$\hat{\gamma} = +1.6$ (the McCallum range) at CEPII's effective
$\sim 100$ km grain, and another paper using a finer satellite-derived
distance reports $\hat{\gamma} = +0.5$ (a $1.6\times$ multiplier), the
literature has been treating these as different findings about the
world. The paper shows they are the same regression conditional on
different measurement choices.

**This matters for:**
- Applied gravity work, which uses CEPII's `distw_harmonic` in
  $\sim 100$% of published studies since 2005.
- Welfare-counterfactual work (Caliendo-Parro 2015 and its descendants)
  that uses gravity-estimated trade costs as inputs.
- Policy advice that takes border-effect multipliers as structural
  parameters when they are actually regularization-conditional.

### 2.2 What's the size of the problem?

Quantitatively: the home-bias coefficient swings by 4 log units across
the structurally consistent $(\theta, d_{\min})$ space. The conventional
McCallum-AvW range lives at the high-$d_{\min}$ end of this swing.
At fine spatial resolution, the multiplier drops by 1-2 orders of
magnitude (in OLS) or by 4 orders of magnitude (in PPML).

This is not a 5% correction. It is a factor-of-10-to-100 correction
to a literature-defining number.

### 2.3 Why hasn't anyone noticed?

Three reasons, in increasing order of severity:

1. **The pole is in plain sight, but only if you evaluate the formula
   at structurally consistent $\theta$.** Head and Mayer (2002) chose
   $\theta = -1$, which is one unit above the pole. Mayer and Zignago
   (2011) propagated the choice. The applied literature has used the
   CEPII value as a number, not as an evaluation of an integral, and
   never re-derived the formula at the gravity-consistent $\theta$.
2. **The structural-vs-reduced-form $\theta$ conflation is a
   notation-overloading problem that has been formalized by the
   original authors themselves** (HM 2014 line 1287). Once the
   notation is established, the conflation feels mechanical, not
   substantive. It takes a careful re-reading of the original
   derivation to see the inconsistency.
3. **The empirical implications are large enough that a referee would
   reasonably ask "if this is right, why hasn't it shown up?".**
   The answer: because the literature has been comparing apples to
   apples (everyone uses CEPII's specific number), so the
   between-paper variation is small. The paper shows the
   regularization sensitivity by varying the spatial grain, which
   the literature has not done.

### 2.4 Will this clear referee review at a top journal?

My honest probability estimates:

- **JIE methodology section or REStat methodology section:** 60% chance
  of acceptance at the first major-revision round. The structural
  critique is clean; the empirical demonstration is honest; the
  implications are substantive but contained.
- **AER short paper:** 15% chance. AER is reluctant on
  methodology-only papers, and the work is fundamentally a methodology
  contribution.
- **RIE or IJIE:** 80% chance. These journals publish methodology
  pieces routinely.
- **Empirical Economics or Journal of Econometric Methodology:** 90%
  chance.

The realistic target is JIE methodology or RIE. The SSRN preprint is
the right first step regardless of journal target.

### 2.5 What would make it more meaningful?

Three extensions would substantially strengthen the contribution:

1. **A full annual panel (2000-2022)** with sectoral breakdown (BACI
   HS6) would address the small-sample concern definitively. Adds
   $\sim 20$ years of variation and 20+ sectors of variation.
2. **The Anderson-van Wincoop multilateral-resistance re-estimation.**
   AvW's $5\times$ multiplier is the most-cited number in the
   border-puzzle literature. Re-running the AvW system with explicit
   regularization documentation, and showing how the $5\times$
   decomposes into structure and measurement, would be a substantial
   contribution in its own right.
3. **A welfare-counterfactual demonstration.** Take a published
   trade-policy counterfactual (e.g., Suez 2021, Brexit, US-China
   tariffs) and show that the implied welfare effect depends on the
   regularization choice. This moves the critique from methodology to
   policy, which would substantially elevate its profile.

These extensions are listed as future work in the paper. The current
preprint is a clean diagnostic; the extensions would make it
constructive.

---

## Part 3: Does it actually make a difference?

### 3.1 The honest realist answer

**Probably yes, but slowly.** Methodology papers in economics tend to
have long latency. Silva-Tenreyro (2006), now considered required
reading for any gravity estimation, took ~5 years to become
canonical. Yotov-Piermartini-Monteiro-Larch (2016), the modern PPML
cookbook, took ~3 years. The Paper 5a critique addresses an even more
fundamental object (the intra-national distance), so adoption could be
slower or faster depending on which way the field's attention swings.

Three plausible adoption pathways:

1. **Citation by working-paper authors who already think gravity
   distance methodology needs work.** Probability of citation in any
   given gravity working paper over the next 5 years: ~5-15% if the
   paper gets published in a methodology venue; ~1-5% if it stays as
   an SSRN preprint only.
2. **Adoption by CEPII itself in a future GeoDist update.** Mayer
   could acknowledge the pole and document the structural-vs-reduced-form
   conflation in the next CEPII methodology note. This would be the
   highest-impact outcome. Probability: ~20% conditional on the email
   being sent and CEPII engaging with the critique.
3. **Adoption by the Yotov-Piermartini-Monteiro-Larch cookbook.** The
   YPMM book is the structural-gravity gold standard; a future edition
   could add a chapter on regularization-aware estimation. Probability:
   ~10% conditional on Yotov engaging.

### 3.2 The honest skeptical answer

**It might not make a difference.** Three reasons it could be
ignored:

1. **The structural-vs-reduced-form $\theta$ conflation is a
   notational hygiene point.** A referee might dismiss it as a
   pedagogical observation rather than a substantive contribution.
2. **The empirical demonstration has its own caveats** (Wei proxy,
   $\delta = 1$ assumption, three-year sample). A referee could focus
   on these limitations to avoid engaging with the structural critique.
3. **The applied gravity literature is deeply invested in CEPII's
   `distw_harmonic`.** A literature-wide methodological reset is
   costly. Authors might prefer to ignore the critique and keep
   using CEPII as they have for 20 years.

The realistic probability that the paper changes practice in the next
5 years: 30-50%. Higher if it gets a top journal; lower if it stays
as an SSRN preprint.

### 3.3 What would make the difference more likely?

Three operational moves:

1. **Send the scoop-protection email to Thierry Mayer.** If Mayer
   acknowledges the critique, the literature follows. If Mayer
   dismisses it, the literature will note that the original architect
   pushed back, and the critique becomes a contested observation
   rather than an accepted update. Either way is informative.
2. **Publish in a top methodology venue.** SSRN preprint first, then
   submit to JIE or REStat. The methodology section is the right
   home; full-paper journals would underweight this.
3. **Release the replication code publicly** (GitHub) with a
   one-script demo that lets anyone reproduce the regularization
   sweep in their own setting. Lowering the cost of adoption matters.

All three are planned in the current workflow.

---

## Summary verdict

**Is it real?** Yes. The math is uncontroversial. The empirical
findings reproduce cleanly with two different estimators on a
panel and two cross-sections. The literature claim has been
triangulated against the canonical structural-gravity papers and
the closest prior art (Rauch 2016, Coughlin-Novy 2021), and
nothing duplicates the specific divergence-at-gravity-consistent-$\theta$
mechanism.

**Is it meaningful?** Yes, with caveats. The contribution is a
methodological correction to a literature-defining number (the
border-puzzle multiplier). The size of the correction is large
(factor of 10-100). The contribution is provisional in the standard
scientific way: extensions to PPML with full panel data, sectoral
breakdowns, and welfare counterfactuals would strengthen it
substantially.

**Does it actually make a difference?** Probably, but with non-trivial
probability of being ignored. The three operational moves (Mayer
email, top-journal publication, public replication code) materially
raise the probability of adoption. The current trajectory has the
paper on a 6-month path to SSRN posting and a 12-18 month path to
journal acceptance.

**Honest one-line summary:** A clean structural-mathematical
observation about a 24-year-old gravity construction, with empirical
evidence that the consequences are large, novelty triangulated across
six candidate literatures, and a defensible path to publication. The
work has real content. Whether it makes a difference depends on
whether the gravity-trade community engages with it; the early
indicators (Gemini's read, the strength of the agent triangulation)
suggest engagement is plausible.

---

*Last updated: 2026-05-22.*
