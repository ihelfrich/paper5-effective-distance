# The border puzzle was partly a measurement artifact

*A plain-English companion to "The Untaken Limit: A Pole at θ = -2 in
the Head-Mayer Closed Form for Intra-National Distance" by Ian T.
Helfrich, Georgia Tech, 2026.*

---

## What economists meant by the "border puzzle"

In 1995, the economist John McCallum looked at trade flows between
Canadian provinces and US states. He found that two Canadian provinces
trade with each other about 22 times more than they trade with a US
state of similar size and similar distance away. Twenty-two times.
This was the "border puzzle." Borders should matter, sure, but a
factor of 22 was hard to explain by tariffs alone (which were small
under NAFTA) or by language and culture (which Canadian provinces
also share with American states).

In 2003, James Anderson and Eric van Wincoop wrote a famous paper
showing that McCallum's 22x was inflated by a missing term in the
gravity model: the "multilateral resistance" each country has to
trade with the rest of the world, not just its immediate partner.
Correcting for this, the multiplier dropped to about 5x. Still a lot,
but more reasonable. The Anderson-van Wincoop paper is one of the
most cited papers in international economics.

For the next 20 years, the "border puzzle" was a settled question.
Estimating gravity regressions with the AvW correction gave a
home-bias multiplier somewhere between 5x and 10x, depending on the
sample and specification. The literature moved on to other things.

## The thing nobody noticed

There's one ingredient in every gravity regression that controls how
much of the border puzzle survives: a number called the
"intra-national distance." This is the distance from a country's
producers to its consumers, measured inside the country's borders.
For the US, this is the typical distance from, say, a factory in
Ohio to a consumer in California: maybe 2,000 km. For Belgium, it's
maybe 100 km.

Where does this number come from? Almost every applied paper uses
the same source: a database called CEPII GeoDist, built by Thierry
Mayer and Soledad Zignago in 2005 and updated in 2011. CEPII gives
you a single number per country, computed from a closed-form
mathematical formula that Keith Head and Thierry Mayer derived in
2002. The formula:

> d_ii = (2/(θ+2))^(1/θ) × R

where R is the equivalent-disk radius of the country (the radius
of a circle with the same area as the country), and θ is a
parameter that the analyst has to choose.

CEPII chooses θ = -1. They are explicit about this: Mayer and
Zignago (2011, page 7) write that "θ equal to -1 ... corresponds
to the usual coefficient estimated from gravity models of bilateral
trade flows." In other words, they choose θ to match the slope
you see in a typical gravity regression of log-trade on log-distance.

Here's the wrinkle. The θ in the closed-form formula and the θ
that gravity regressions estimate are *not the same thing*. The
formula's θ is the **structural integration parameter**, derived
from the underlying CES (constant elasticity of substitution) model
of consumer preferences. The structural θ is mathematically equal
to δ × (1 - σ), where δ is the "iceberg cost exponent" (typically
1) and σ is the elasticity of substitution between varieties
(typically 4 to 8 in the trade literature). So the structural θ
is typically between -7 and -3.

The "θ" that gravity regressions estimate is a reduced-form
regression slope. That number is around -1, give or take.

The two θ's happen to have the same name. They are not the same
quantity.

CEPII chose the reduced-form value (-1) to evaluate the structural
formula. That is the wrong-θ conflation. Head and Mayer themselves
did it first, in their 2002 paper. They wrote that "hundreds of
gravity equation estimates of θ" support the choice θ ≈ -1, citing
the regression slope as evidence for the structural parameter.

## Why this matters

The closed-form formula has a pole at θ = -2. That means: the
formula blows up at θ = -2, and for θ < -2, the formula is
complex-valued (it produces a number with an imaginary part).
At the structurally consistent θ values (-7 to -3), the formula
doesn't give you a real number at all.

Now, CEPII doesn't actually use the closed form. They compute the
intra-national distance numerically, by averaging pairwise distances
between major cities in each country, weighted by city population.
This gives them a real number even at θ values where the closed
form would be complex. But here's the catch: the numerical
calculation also depends on the spatial resolution of the city
sample. With more cities at finer resolution, the implied distance
gets *smaller*. With fewer cities at coarser resolution, the
implied distance gets *larger*. There's no fixed point. The
"intra-national distance" of the United States is not a number;
it's a function of how finely you sample.

I show in the paper that on a uniform-disk simulation, the implied
intra-national distance at the structurally consistent θ = -5 ranges
from about 6 km (at fine spatial resolution) to about 640 km (at
coarse resolution). The CEPII value at θ = -1 happens to be about
500 km, which is near the high end of this range. CEPII's choice
is consistent with a coarse-resolution evaluation of an object that,
if you take the math seriously, doesn't actually exist as a single
number.

## What happens to the border puzzle

The home-bias multiplier in a gravity regression depends directly
on the intra-national distance. If you swap CEPII's number for a
finer-resolution one, the regression's prediction of intra-national
trade goes up (because the predicted trade volume is inversely
related to distance), and the home dummy adjusts downward to keep
the prediction close to the actual trade volume.

In a standard gravity regression on 2010 trade data (the cleanest
sample I have), the home-bias multiplier moves from:

- **0.05x** at fine spatial resolution (1 km grid), meaning intra-national
  trade is *less* than what gravity predicts.
- **1.0x** at moderate resolution (around 10-20 km grid), meaning
  intra-national trade matches the gravity prediction exactly.
- **8.4x** at the CEPII baseline (~500 km equivalent grid), the
  canonical McCallum-Anderson-van Wincoop "border puzzle"
  multiplier.

The "border puzzle" of 8.4x in the literature corresponds to a
specific coarse-resolution choice of how to construct intra-national
distance. If you take the structural math seriously and use a finer
resolution, the puzzle disappears (or even reverses sign).

This doesn't mean borders don't matter. It means that the **specific
quantification** of how much they matter is conditional on a
measurement choice that the literature has not engaged with.

## What to do about it

Three things.

**First**, papers that report border-puzzle multipliers should
report the spatial resolution at which they computed intra-national
distance. A multiplier of 10x at 100 km grain and 5x at 10 km grain
are not contradictory results; they are the same regression under
different measurement choices.

**Second**, the convention of treating intra-national distance as
a fixed property of a country's geometry, rather than as a
discretization choice, has misrepresented the literature's
empirical content. New work should acknowledge the
structural-versus-reduced-form distinction in the integration
parameter θ.

**Third**, the future of structural gravity estimation may be in
techniques that do not require pinning intra-national distance to
a single number. The Caliendo-Parro (2015) tetrad-estimator
methodology cancels distance entirely. The Eaton-Kortum (2002)
framework assumes intra-national distance equals 1 without
integrating over space. The Costinot-Donaldson-Komunjer (2012)
approach absorbs all bilateral trade costs in importer-exporter
fixed effects. These are not workarounds; they are honest
acknowledgments that the structural integral does not exist as a
fixed object when the structural θ is in its empirically
defensible range.

## A note on what's NEW here

The pole at θ = -2 in the Head-Mayer (2002) closed form has been
sitting in the published literature for 24 years. The
structural-versus-reduced-form θ conflation is on lines 502-505 of
Head and Mayer's 2002 paper, in plain English. CEPII's choice of
θ = -1 has been documented since 2005. Nobody appears to have
connected these dots and asked what happens at the structurally
consistent θ.

This paper makes that connection. The math is not new. The data is
not new. The novelty is reading what is already there.

---

*Full technical paper at SSRN [link to come]. Replication code at
github.com/ihelfrich/paper5-effective-distance [to come]. Companion
teaching guide volumes 1 and 2 cover the math, the empirical pipeline,
and the GIS construction for readers who want to work through the
derivation themselves.*

*Comments, critiques, and corrections welcome at ianthelfrich@gmail.com.*
