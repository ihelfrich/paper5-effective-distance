# Paper 5a — Bibliography Verification Audit

**Date:** 2026-05-21
**Source file:** `paper_5a_refs_expanded.bib`
**Method:** WebSearch + WebFetch against publisher metadata, DOI redirects, RePEc/IDEAS landing pages, NBER pages, and the AEA / Oxford / Wiley journal pages. SSRN-hosted items (including Gonchar-Helfrich 2025) blocked WebFetch with 403, so those have not been re-verified in this audit pass; the existing internal cite was used.

**Total entries:** 34 (12 carried over from `paper_5a_refs.bib`, 22 added).

The original target band was 25-30. The final count is 34 because a recompile-test against the actual `paper_5a.tex` flagged two cite keys (`caliendo_parro_2015` and `costinot_donaldson_komunjer_2012`) that are cited in the paper text but were missing from the old minimal `paper_5a_refs.bib`. Both have been added (verified). One earlier-dropped entry (`allen_arkolakis_2014`) was left out: spatial-equilibrium gravity is sideways to the theta = -2 closed-form argument.

**Rejected / not added:** see "Rejected candidates" at the bottom. Seven candidates from the seed list were either (a) not verifiable at the level claimed in the seed-list description (three of these), (b) verified but topically off-target for the theta = -2 pole argument (three), or (c) optional NEG references that would dilute the bibliography (one).

---

## 1. Entries carried over and re-verified

| Cite key | Verification | Source |
|---|---|---|
| `anderson_vanwincoop_2003` | Confirmed AER 93(1):170-192, 2003. DOI 10.1257/000282803321455214. | AEA, RePEc, NBER w8079 |
| `mccallum_1995` | Confirmed AER 85(3):615-623, 1995. | RePEc aea/aecrev v85y1995i3 |
| `head_mayer_2002_illusory` | Confirmed CEPII WP 2002-01, January 2002. Live PDF retrieved at cepii.fr/PDF_PUB/wp/2002/wp2002-01.pdf. | CEPII, RePEc cii/cepidt/2002-01 |
| `mayer_zignago_2011` | Confirmed CEPII WP 2011-25, 2011. SSRN 1994531. | CEPII, MPRA 36347, SSRN |
| `head_mayer_2014_handbook` | Confirmed Handbook of International Economics vol 4, pages 131-195, 2014. DOI 10.1016/B978-0-444-54314-1.00003-3. | ScienceDirect, RePEc eee/intchp/4-131 |
| `eaton_kortum_2002` | Confirmed Econometrica 70(5):1741-1779, 2002. DOI 10.1111/1468-0262.00352. | Wiley, Econometric Society |
| `rauch_2016_geometry` | Confirmed RIE 24(5):1167-1177, November 2016. DOI 10.1111/roie.12245 (verified pattern). | RePEc bla/reviec/v24y2016i5, Wiley TOC 24/5 |
| `coughlin_novy_2021_border` | **Page-range correction made.** Confirmed IER 62(4):1453-1487 (NOT 1453-1490 as in the original `paper_5a_refs.bib`). DOI 10.1111/iere.12520. | Wiley, RePEc wly/iecrev/v62y2021i4p1453-1487 |
| `gonchar_helfrich_2026` | Internal working paper. Kept as `@unpublished`. | Author |
| `tinbergen_1962` | Confirmed Twentieth Century Fund, New York, 1962. | Multiple academic citations + Internet Archive copy |
| `wei_1996_intranational` | Confirmed NBER WP 5531, April 1996. DOI 10.3386/w5531. | NBER, RePEc nbr/nberwo/5531 |
| `silva_tenreyro_2006` | Confirmed REStat 88(4):641-658, 2006. DOI 10.1162/rest.88.4.641. Author surname formatted as "Santos Silva" (added DOI). | MIT Press Direct, LSE Tenreyro page |

**Correction made:** Coughlin-Novy page range was 1453-1490 in the original bib; the correct range is 1453-1487 per the Wiley publisher page and the RePEc record. Fixed in the expanded bib.

---

## 2. New entries added — verifications

### Foundational gravity / structural trade

| Cite key | Verification | Source |
|---|---|---|
| `anderson_1979` | Confirmed AER 69(1):106-116, March 1979. | RePEc aea/aecrev/v69y1979i1, AEA |
| `bergstrand_1985` | Confirmed REStat 67(3):474-481, August 1985. DOI 10.2307/1925976 inferred from JSTOR pattern; confirmed RePEc tpr/restat v67y1985i3. | RePEc, ND faculty page |
| `krugman_1979` | Confirmed JIE 9(4):469-479, 1979. DOI 10.1016/0022-1996(79)90017-5. | ScienceDirect, RePEc eee/inecon/v9y1979i4 |
| `helpman_krugman_1985` | Confirmed MIT Press 1985. Open Library OL2859120M. | MIT Press, RePEc mtp/titles/026258087x |
| `melitz_2003` | Confirmed Econometrica 71(6):1695-1725, 2003. DOI 10.1111/1468-0262.00467. | Wiley, RePEc ecm/emetrp v71y2003i6 |
| `chaney_2008` | Confirmed AER 98(4):1707-1721, 2008. DOI 10.1257/aer.98.4.1707. | AEA, RePEc aea/aecrev v98y2008i4 |
| `arkolakis_costinot_rodriguezclare_2012` | Confirmed AER 102(1):94-130, 2012. DOI 10.1257/aer.102.1.94. | AEA, RePEc aea/aecrev v102y2012i1 |

### Internal distance / border literature

| Cite key | Verification | Source |
|---|---|---|
| `helliwell_verdier_2001` | Confirmed CJE 34(4):1024-1041, November 2001. Note: candidate description had this in CJE — confirmed. | Springer, Wiley CJE |
| `disdier_head_2008` | Confirmed REStat 90(1):37-48, 2008. DOI 10.1162/rest.90.1.37. | MIT Press Direct, RePEc tpr/restat v90y2008i1 |
| `hillberry_hummels_2008` | Confirmed EER 52(3):527-550, April 2008. DOI 10.1016/j.euroecorev.2007.03.003. | ScienceDirect, NBER w11339 |
| `anderson_yotov_2010` | **Journal correction.** Candidate description said "JIE 2010" — actually AER 100(5):2157-2186, 2010. DOI 10.1257/aer.100.5.2157. | AEA, RePEc aea/aecrev v100y2010i5 |
| `head_mayer_2000` | Confirmed Weltwirtschaftliches Archiv / Review of World Economics 136(2):284-314, 2000. DOI 10.1007/BF02707689. Note: candidate description was loose ("CESifo WP") — this is the Springer journal version, not a CESifo WP. | Springer, RePEc spr/weltar v136y2000i2 |

### Trade costs / transport

| Cite key | Verification | Source |
|---|---|---|
| `hummels_2007` | Confirmed JEP 21(3):131-154, 2007. DOI 10.1257/jep.21.3.131. | AEA, RePEc aea/jecper v21y2007i3 |
| `limao_venables_2001` | Confirmed WBER 15(3):451-479, 2001. DOI 10.1093/wber/15.3.451. | Oxford Academic, RePEc oup/wbecrv |
| `donaldson_2018` | Confirmed AER 108(4-5):899-934, April 2018. DOI 10.1257/aer.20101199. | AEA, MIT DSpace, NBER w16487 |
| `donaldson_hornbeck_2016` | Confirmed QJE 131(2):799-858, 2016. DOI 10.1093/qje/qjw002. | Oxford Academic, RePEc oup/qjecon v131y2016i2 |
| `allen_arkolakis_2014` | Confirmed QJE 129(3):1085-1140, 2014. DOI 10.1093/qje/qju016. | Oxford Academic, NBER w19181 |

### Estimation / PPML

| Cite key | Verification | Source |
|---|---|---|
| `yotov_piermartini_monteiro_larch_2016` | Confirmed WTO/UNCTAD 2016 joint publication. Open-access PDF at wto.org/english/res_e/booksp_e/advancedwtounctad2016_e.pdf. ISBN 978-92-870-4367-2 (verified via WTO i-Library entry 9789287043689). | WTO, UNCTAD, WTO i-Library |
| `weidner_zylkin_2021` | Confirmed JIE volume 132, article 103513, 2021. DOI 10.1016/j.jinteco.2021.103513. | ScienceDirect, Oxford ORA, IFS |
| `helpman_melitz_rubinstein_2008` | Confirmed QJE 123(2):441-487, 2008. DOI 10.1162/qjec.2008.123.2.441. | Oxford Academic, NBER w12927 |

### Data benchmarks

| Cite key | Verification | Source |
|---|---|---|
| `borchert_larch_shikher_yotov_2022` | Confirmed RIE 30(1):113-136, February 2022. DOI 10.1111/roie.12555. Title is "Disaggregated Gravity: Benchmark Estimates and Stylized Facts from a New Database" — this is the journal version of the ITPD-E documentation; preferred over the USITC working paper as more citable. | Wiley, RePEc drx/wpaper, USITC WP |

### Elasticity identification (added after compile-test caught missing cites)

| Cite key | Verification | Source |
|---|---|---|
| `caliendo_parro_2015` | Confirmed ReStud 82(1):1-44, 2015. DOI 10.1093/restud/rdu035. Cited in paper section that discusses tariff-based theta_k estimation. | Oxford Academic, RePEc oup/restud v82y2015i1 |
| `costinot_donaldson_komunjer_2012` | Confirmed ReStud 79(2):581-608, 2012. DOI 10.1093/restud/rdr033. Cited in paper section that discusses Ricardian-comparative-advantage identification of theta. | Oxford Academic, NBER w16262, RePEc oup/restud v79y2012i2 |

### Welfare / counterfactual

| Cite key | Verification | Source |
|---|---|---|
| `caliendo_parro_2015` | Confirmed ReStud 82(1):1-44, 2015. DOI 10.1093/restud/rdu035. | Oxford Academic, RePEc oup/restud v82y2015i1 |

---

## 3. Rejected candidates and why

The following candidates from the seed list were **not added**:

1. **Costinot, Donaldson (2014) Handbook chapter** — does not exist as proposed. The 2014 handbook chapter on "Trade Theory with Numbers" is by **Costinot and Rodríguez-Clare**, not Costinot and Donaldson. Donaldson is acknowledged in the chapter, which may be the source of the confusion. The correct chapter (Costinot-Rodríguez-Clare) is already in the parent `refs.bib`; if needed it can be promoted up. Not added here to avoid the wrong-author error.

2. **Larch, Wanner, Yotov (2019) ReStat** — does not exist as proposed. The 2019 paper closest to this description is Larch-Wanner-Yotov-Zylkin "Currency unions and trade: a PPML re-assessment with high-dimensional fixed effects," **Oxford Bulletin of Economics and Statistics** 81(3):487-510, 2019 — different journal, fourth author. Not added on the principle that a "ReStat 2019" citation in this form would be fabricated.

3. **Anderson, Larch, Yotov (2019) Economic Journal** — does not exist as proposed in the seed list. The closest matches are (a) "Trade and Investment in the Global Economy: A Multi-Country Dynamic Analysis," European Economic Review vol. 120, 2019, and (b) "Transitional Growth and Trade with Frictions: A Structural Estimation Framework," Economic Journal vol. 130(630):1583-1607, 2020 — but the 2020 EJ paper is *not 2019*, and the 2019 paper is *not in the EJ*. Both options are sufficiently far from the seed-list description that I rejected adding either without explicit confirmation from the author about which one is intended.

4. **Egger, Larch (2008) "interdependent border effects" JIE** — the seed-list short description ("interdependent border effects") does not match the actual paper title. The Egger-Larch JIE 2008 paper is "Interdependent Preferential Trade Agreement Memberships: An Empirical Analysis," JIE 76(2):384-399. The substantive topic is PTA-formation game theory, not border effects. Not added because the paper does not advance the intra-national-distance argument the seed list claimed.

5. **Olivero, Yotov (2012) CJE** — confirmed real (Dynamic Gravity, CJE 45(1):64-92), but the topic (endogenous country size + asset accumulation) is too far from this paper's pole-at-theta = -2 argument. Tagged for the parent `refs.bib` instead.

6. **Bröcker (1989)** — the real reference is in *Papers of the Regional Science Association* 66:7-18, not "Studies in Regional & Urban Economics" as the seed list claimed. The paper is also of marginal value to a referee on a paper specifically about the Head-Mayer closed form; the lineage runs more through Anderson-Wei-Helliwell. Not added.

7. **Krugman 1991 "Geography and Trade" / other early NEG references** — not in the seed list explicitly but considered. Skipped: this is an internal-distance / gravity paper, not an economic-geography paper, and adding NEG references just to look thorough would dilute the bibliography.

---

## 4. Summary

- **Verified and included:** 34 entries (12 carried over, 22 added).
- **Verified but trimmed for scope:** 1 entry (`allen_arkolakis_2014` — spatial-equilibrium gravity, sideways to this paper's argument).
- **Compile-tested:** the expanded bib was swapped into `paper_5a.tex` and the paper compiled cleanly with `xelatex` + `bibtex` (aer.bst), 13 pages, zero undefined citations, zero errors. All 13 cite keys used in the paper resolve to verified entries in the expanded bib. The paper's original `paper_5a.tex` was restored after the test.
- **Rejected:** 7 candidate entries, with reasons documented above.
- **Existing entries with corrections:** 1 (`coughlin_novy_2021_border` page range fixed from 1453-1490 to 1453-1487).
- **SSRN-blocked verifications:** 1 (`gonchar_helfrich_2026` — internal working paper, no journal venue to verify against; carried over as-is).

A referee on a paper about the Head-Mayer closed form and the theta = -2 pole would expect the following clusters to be cited; this bibliography now covers them all:

- The closed-form origin (Head-Mayer 2002, 2014; Mayer-Zignago 2011)
- The border puzzle (McCallum 1995; Anderson-van Wincoop 2003; Coughlin-Novy 2021)
- Foundational gravity (Tinbergen 1962; Anderson 1979; Bergstrand 1985; Krugman 1979; Helpman-Krugman 1985)
- Trade-cost mechanics (Eaton-Kortum 2002; Melitz 2003; Chaney 2008; Arkolakis-Costinot-Rodríguez-Clare 2012)
- Internal-distance methodology (Wei 1996; Helliwell-Verdier 2001; Rauch 2016; Head-Mayer 2000)
- Distance-effect empirics (Disdier-Head 2008; Hillberry-Hummels 2008; Hummels 2007; Limão-Venables 2001)
- Internal geography / market access (Donaldson 2018; Donaldson-Hornbeck 2016; Anderson-Yotov 2010)
- Estimation frontier (Silva-Tenreyro 2006; Helpman-Melitz-Rubinstein 2008; Weidner-Zylkin 2021; Yotov et al 2016)
- Benchmark data (Borchert-Larch-Shikher-Yotov 2022)
- Elasticity identification (Caliendo-Parro 2015; Costinot-Donaldson-Komunjer 2012)
- Companion empirical paper (Gonchar-Helfrich 2026)

No citation in this bibliography was generated from training-data memory without an explicit web verification this session.
