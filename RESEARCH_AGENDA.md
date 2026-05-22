# Research agenda — questions the raster-distance pipeline lets us actually answer

*Drafted 2026-05-20, autonomous-run mode. This is a brainstorm document: I am writing down every direction I would seriously consider, not just the safe ones. The "do this next" filtering happens at the bottom.*

The infrastructure being built for Paper 5 (raster-resolution effective distance, directional asymmetry, Wasserstein velocity) is more general-purpose than Paper 5 alone. The same machinery answers many other open questions in international trade and economic geography, several of which are more important than the home-bias / border-effect question that motivated the build.

## Eleven questions I think are worth actually answering

### 1. The Disdier-Head puzzle: why hasn't the distance elasticity fallen?

Disdier and Head (2008) collected ~1,000 gravity estimates and showed that the distance elasticity has been remarkably stable around -1 since 1870, even as transport costs have collapsed (steam shipping, containerization, air freight, the internet). This is "missing globalization." Numerous explanations exist (heterogeneous-firm composition, intermediate goods, services). None of them are fully convincing.

A possibility the existing literature does not consider: *effective distance has not actually fallen as much as geographic distance suggests*. As internal density redistributes (urbanization, suburbanization, coastal migration in some countries, interior development in China), the producer-consumer effective distance might have stayed roughly constant in the CES sense even as nominal transport costs collapsed. This is testable with the Li et al. nightlights series 1992-2024 plus GHS-POP 2000-2020.

**Specifically.** For each country-pair, compute $d_{ij}^{\eff}(t, \theta = -5)$ year-by-year for 1992-2022 using harmonized nightlights as the activity proxy. Compute the time series of bilateral effective distance. Compare against the time series of (a) geographic centroid distance (constant) and (b) implied trade cost from gravity residuals. The prediction: $d_{ij}^{\eff}$ has *not* fallen, or has fallen by far less than transport-cost reductions would imply.

If this works, it is a serious paper. The Disdier-Head puzzle has hundreds of citations. A clean resolution would be QJE-tier.

### 2. The China shock at sub-national resolution

The Autor-Dorn-Hanson (2013) "China shock" literature measures local exposure to Chinese imports by weighting national import shares with local industry composition. This is the standard shift-share approach. What it misses: *the geography of who actually trades with whom*. A US PUMA's exposure to China depends not just on its industry mix but on where in China its industry's suppliers are located, and where in the US that PUMA sits relative to the ports through which Chinese goods enter.

The pipeline lets you compute, for each US PUMA $p$ and each Chinese province $c$, the effective distance $d_{pc}^{\eff}$ at industry-relevant $\theta$. The "real" China shock to PUMA $p$ is then weighted not just by national-level Chinese imports but by the bilateral *effective* trade cost. This changes the implied exposure measure for inland PUMAs vs coastal ones; for North-East rustbelt PUMAs vs South-East ones.

**The strong version:** the resolution to the "China shock placebo" debate (Borusyak-Hull-Jaravel 2025) may live in the geography. If shocks identified by the standard ADH instrument are correlated with effective-distance-driven exposure, that's an explanation for the residual placebo correlation that doesn't require throwing out the result.

### 3. Sanctions-induced trade rerouting

Sanctions on Russia (2014, 2022) created a natural experiment in effective-distance change. Russia's exports rerouted through Turkey, China, the UAE. The effective distance from Russia to Europe rose sharply (not because Russia moved, but because the intermediary network changed). The effective distance from Russia to China shrank (more routes, more capacity).

The pipeline lets you compute $d_{ij}^{\eff}$ before and after the sanctions episode, and ask: does the change in $d_{ij}^{\eff}$ predict the change in $X_{ij}$? This is a textbook event study where the treatment (sanctions) operates through a measured channel (effective distance) that nobody else has measured at this resolution.

**Why this matters:** the sanctions literature is policy-relevant and well-funded. A paper that quantifies the "rerouting cost" using effective distance has direct policy implications for the design of future sanctions regimes.

### 4. Sea-level rise and trade

The IPCC AR6 puts coastal flooding under various SLR scenarios. The 2°C warming scenario displaces ~50M people from low-elevation coastal zones. As populations migrate inland, internal density shifts away from ports.

The pipeline lets you compute, under each SLR scenario, the post-displacement effective distance $d_{ij}^{\eff}(2050, \text{SSP2-4.5})$ and so on. The economic question: how does projected trade volume shift if internal density adjusts to SLR but ports remain in their current locations? This is a "structural projection" exercise of the kind Conte-Desmet-Nagy-Rossi-Hansberg (2021) do at 1° resolution, but at 1km and with effective-distance machinery.

**Side benefit:** this generates a publication-quality climate-economics paper that uses none of Paper 5's home-bias machinery but all of its raster infrastructure. Diversifies the publication portfolio.

### 5. Civil war and the fragmentation of "country"

When Syria fragments into government-held, opposition-held, ISIS-held, and Kurdish-held zones (2014-2018), the standard country-level gravity equation treats Syria as a coherent unit. It is not. Internal trade between Aleppo and Damascus is essentially zero during certain windows.

The pipeline lets you compute *fragmented effective distances*: for a country with internal conflict zones, decompose the population into accessible and inaccessible cells, recompute $d_{ii}^{\eff}$ with the conflict zones treated as separate units. The prediction: the gravity-implied "trade collapse" of conflict countries is partly an artifact of treating $d_{ii}$ as continuous when it has been topologically broken.

This connects to the Krugman-Venables agglomeration literature and to recent work on "fragile states" trade.

### 6. The agglomeration-trade duality

Krugman-Venables theory predicts that trade integration and internal agglomeration are jointly determined. As trade costs fall, mass concentrates near borders/ports. As mass concentrates, the effective distance to partners shifts. There is a fixed-point structure here.

The pipeline lets you measure the agglomeration side ($\vec{c}_i(t)$ trajectory from GHS-POP / nightlights time series) and the trade side ($X_{ij}(t)$ from BACI). The question: is the *speed* of internal agglomeration predicted by exposure to specific trading partners' effective distance? This is the bilateral version of the Krugman-Venables prediction, with directionally-specific identification.

This becomes Paper 5b's Prediction P2 (Wasserstein velocity) but in reverse: instead of asking whether velocity predicts trade growth, ask whether trade growth predicts velocity. The two should be jointly identified.

### 7. Trade in services: what is "distance" for digital exports?

The services-gravity literature (Anderson-Milot-Yotov 2018, Borchert et al. 2024) consistently finds smaller distance elasticities for services than goods. This is unsurprising — banking, consulting, software cross borders with minimal physical friction. But it is also under-theorized: what *is* the distance object for a service?

A possible answer: for digitally-deliverable services, distance is *time zone overlap* and *language barrier* and *regulatory-regime distance*, not km. For tourism services, distance is *travel cost*. For business consulting, distance is *meeting-day cost* which loads on time-zone differential.

The pipeline can be re-purposed: replace the population raster with a *time-zone raster* (or a regulatory-regime raster, or a language raster), and compute the CES effective distance over that. The question: do services trade flows respond to these alternative distance objects in the way the goods literature responds to geographic distance?

### 8. The environmental footprint of trade

Shipping emissions are roughly proportional to ton-km. The standard "carbon footprint of trade" measures multiply trade volume by an emissions intensity per ton-km. But the relevant ton-km is not the great-circle distance; it is the actual shipping distance, weighted by the producer and consumer mass densities. The pipeline computes exactly this object.

**Prediction:** the carbon footprint of trade between two countries depends not just on trade volume but on the *effective* distance, which can differ from geographic distance by 30-50% for some pairs (US-China, EU-South-America, intra-EU). Existing footprint calculations are biased in a measurable direction.

This is a clean, narrow paper. Targets *Nature Climate Change* or *Environmental Science & Policy*.

### 9. The "great trade collapse" of 2008-09 was geographically uneven

Trade fell ~20% in late 2008. The decline was concentrated in durables, in long-supply-chain manufacturing, in countries dependent on trade credit. But it was also geographically uneven in a way that hasn't been fully explained: some country pairs collapsed more than others, conditional on industry mix.

Hypothesis: pairs with *longer effective distance* (more intermediate handling, longer transport spells) suffered larger collapses because longer supply chains amplify shocks. This is the Bems-Johnson-Yi (2010, 2011) intermediate-goods explanation re-cast as a measurement question.

The pipeline computes the effective-distance variation across pairs at the right resolution. A panel regression on the trade-collapse data ($\Delta X_{ij, 2008-09}$ on $d_{ij}^{\eff}$ controlling for industry mix and standard gravity) would test this directly.

### 10. The "missing demand" puzzle in agglomeration economics

The agglomeration literature (Glaeser, Henderson) finds that big cities have wage premiums, productivity premiums, output premiums. The standard explanation is local knowledge spillovers, thick labor markets, intermediate-input variety. But there is a "missing demand" piece: big cities are also where the consumers are. A producer's effective distance to consumer mass is shorter if they sit in a large city than if they sit at the country's geographic median.

The pipeline computes producer-to-consumer effective distance at sub-national resolution. Question: do firms in cities with shorter effective distance to *the rest of the country* enjoy a productivity premium that survives the standard agglomeration controls? This is a within-country horse-race between knowledge-spillover and demand-proximity explanations of the city wage premium.

### 11. Power-law agglomerations and the harmonic-mean property

A theoretical-empirical hybrid. The bilateral d_eff(θ=-5) collapse for adjacent countries (CHN-RUS, MEX-USA in the nb13 v3 result) is a property of the generalized mean at $\theta < 0$. As $\theta \to -\infty$, $M_\theta(d) \to \min(d)$. At $\theta = -5$ with spatially-concentrated atom clouds, the min-distance pair (border-crossing) dominates.

This connects to the structure of city-size distributions (Zipf's law, power-law tails). For a country with power-law city sizes, the harmonic-mean distance is essentially the distance between the two largest cities. For a country with uniform spatial distribution, the harmonic-mean distance is much larger.

**The deeper question:** does the gravity-relevant effective distance depend systematically on the *Zipf coefficient* of the country's city-size distribution? This is a theoretical claim about how the shape of urban hierarchies shows up in trade elasticities. Testable on the Bluhm-Krause (2022) corrected nightlights data.

If true, it generates a sharp prediction: countries with steeper Zipf coefficients should show systematically smaller $d^{\eff}_{ii}$ at $\theta < 0$, and therefore look "smaller" to gravity than their geographic area suggests. Pareto-tail-corrected effective distance is a different object than uncorrected effective distance, and that difference is testable.

## Strategic ranking

| Q | Question | Difficulty | Time to publishable draft | Tier |
|---|---|---|---|---|
| 1 | Disdier-Head puzzle | 4 | 4-6 months | QJE/AER aspiration |
| 2 | Sub-national China shock | 4 | 6-9 months | AER/QJE |
| 3 | Sanctions rerouting | 3 | 3-4 months | JIE/EER |
| 4 | Sea-level rise + trade | 4 | 6-8 months | NCC + JEEM |
| 5 | Civil-war fragmentation | 3 | 4-6 months | JIE/REStat |
| 6 | Agglomeration-trade duality | 5 | 12+ months | QJE/REStud |
| 7 | Services gravity | 3 | 3-5 months | JIE/RIE |
| 8 | Environmental footprint | 2 | 2-3 months | NCC/Env Sci |
| 9 | 2008-09 trade collapse | 3 | 3-5 months | JIE/EER |
| 10 | Missing demand in cities | 4 | 8-12 months | JEEA/JUE |
| 11 | Zipf and effective distance | 3 (theory+empirics) | 4-6 months | JIE/JEcG |

## My honest recommendation

If I had to commit to one paper *outside Paper 5* to start now, in priority order:

1. **Question 1 (Disdier-Head).** It's the question with the largest existing audience, the cleanest test design, and the most rigorous answer the data infrastructure can produce. The result either confirms missing globalization is a measurement artifact (huge) or rejects it (clean null with a quantified upper bound on how much measurement can explain). Both are publishable.

2. **Question 3 (sanctions rerouting).** Most policy-relevant, fastest to publication, lowest risk of being scooped. The 2022 Russia sanctions are recent enough that the trade rerouting data is just becoming clean. The window is open for ~12-18 months.

3. **Question 8 (environmental footprint).** Quickest paper. Two months of work for a *Nature Climate Change*-tier publication. Diversifies the portfolio out of pure economics. The methodology hub is shared with Paper 5.

I would not pursue Questions 6 (agglomeration-trade duality) or 10 (missing demand) right now even though they are the deepest. Both require longer commitment than the current tenure-track positioning allows. Bank the infrastructure, come back to them after Paper 5 lands.

## What I would NOT do

- Don't try to replicate Allen-Arkolakis 2014 or 2022 immediately. The road-network data isn't on disk, and the model setup takes 6+ months to internalize. Use the AA framework as a citation target, not as a replication.
- Don't try to extend Eaton-Kortum 2002 immediately. The Ricardian comparative-advantage object requires productivity data that's harder to get than the trade-cost object.
- Don't pursue purely theoretical work (no empirics) while the data infrastructure is still being built. The leverage right now is empirical, not theoretical.
- Don't write a literature review paper. Liz and I have ideas; the audience wants results.

## Connection to Liz's interests

Liz's V2 sent over earlier indicates focus on the US panel and on the bilateral specification. Questions 2 (sub-national China shock) and 9 (2008-09 trade collapse) are the strongest fits with that. Question 1 (Disdier-Head) is a stretch fit but the global panel scope might appeal.

Worth a direct conversation with her about whether she wants to be lead on a joint Paper 5b or on one of the new directions above. Don't assume.

## Open methodological pieces this list requires

- **Time-series of $d_{ij}^{\eff}(t)$.** Currently the pipeline is cross-sectional. Question 1 needs the full time series from 1992. Need to extend the country-raster builder to iterate over the Li et al. annual nightlights.
- **Sub-national to sub-national bilateral.** Currently $d_{ij}^{\eff}$ is country-pair. Question 2 needs PUMA-to-province granularity. Need to extend the region masker to handle GADM admin-1 / admin-2 polygons.
- **Conflict zones as polygon overlays.** Question 5 needs conflict-control-zone polygons (ACLED, ETH GROW Up) to overlay on the country masks.
- **Sea-level rise scenarios as raster modifications.** Question 4 needs a way to "remove" cells from the country mask based on elevation × SLR.

All four are 2-4 week engineering items. None block starting the conversations.
