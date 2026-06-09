# Ship plan, written 2026-06-09

Single source of truth for getting Paper 5a posted and the effective-distance
panel released. Supersedes the scattered next-steps lists in
GH_REPO_INIT_PLAN.md and SESSION_NOTES_2026-05-20.md.

## Where things stand (verified today)

- Paper 5a manuscript is complete: 17 pages, 6 sections plus appendix,
  35 verified bibliography entries, figures final (May 22). Recompiled
  clean from source today with a 3-pass xelatex build, zero undefined
  references.
- CRITICAL_ASSESSMENT.md (May 22) verdict stands: math verified, empirics
  reproduce under two estimators, novelty triangulated against six
  literatures.
- Repo is public on GitHub. The `_private/` folder was accidentally
  committed in the initial commit and pushed; on 2026-06-09 it was purged
  from history and force-pushed. The folder now lives only on local disk.
- LICENSE (MIT), DATA_LICENSE.md (CC-BY-4.0 for derived data), and
  CITATION.cff added 2026-06-09.

## Phase 1: post Paper 5a (target: this week)

Steps only Ian can do, roughly one hour total:

1. Read the compiled `manuscript/paper_5a/paper_5a.pdf` one final time.
2. Send the scoop-protection email. Drafts (local only):
   `_private/scoop_protection_email_draft.md` (Mayer) and
   `_private/scoop_protection_email_draft_v0.2_2026-05-21.md` (Yotov).
   Pick one recipient; sending both versions to both people would be odd.
3. Post to SSRN (free, screening typically under 3 business days).
   Metadata is pre-filled in `_private/SSRN_METADATA.md`, abstract updated
   2026-06-09 to match the current panel-PPML version of the paper.
   Networks: ERN International Trade eJournal and Econometric Modeling:
   International Economics eJournal.

Steps Claude can do once the SSRN ID exists:

4. Add the SSRN URL to the paper footnote, recompile, commit.
5. Tag v1.0.0, create the GitHub release with paper PDF plus visual
   abstract, enable the GitHub-to-Zenodo webhook so the release gets a DOI.
6. Submit the RePEc listing so citations track where economists look.

## Phase 2: the effective-distance panel (joint with Liz)

This is the Gonchar-Helfrich dataset play: the full panel extension named
in CRITICAL_ASSESSMENT.md as the highest-value follow-on, and the scale-up
of the 2025 SSRN seed paper "Trade in the Spotlight" (SSRN 5202676).

Known blockers, in order:

1. WorldPop raster pull was at 42 percent on May 20. Resume and finish.
2. Notebook 13 (real-country d_ii at all theta values) stalled on a
   geopandas environment problem. Fix the venv, re-run for the full
   country set.
3. The top-N atom subsampling bias documented in METHODOLOGY_NOTES.md
   needs the planned 5 km block-coarsening robustness run before any
   public data release.

Release stack (decided 2026-06-09 after live venue research):

- Zenodo as primary deposit: 50 GB limit fits, concept DOI plus version
  DOIs fits an annually updated panel, GitHub integration auto-archives
  releases. Data CC-BY-4.0, code MIT.
- Descriptor paper: Nature Scientific Data as primary target (Data
  Descriptor format; APC roughly $2,690; median 162 days submission to
  acceptance; the Li et al. 2020 harmonized nightlights descriptor is the
  structural model and has 500+ citations since 2020). Fallback venue:
  International Economics (CEPII-affiliated Elsevier journal).
- Note from the venue research: CEPII's own gravity database paper was
  never journal-published. It is a working paper with on the order of
  1,800 citations. The citations come from the data being good, free, and
  documented. Venue prestige is secondary to that.

## Decision log

- 2026-06-09: "Start over" rejected. The May 22 critical self-audit found
  the paper real, novel, and publishable; restarting would discard
  verified work. The correct fix was shipping infrastructure, not new
  research.
- 2026-06-09: Ship order is 5a first, panel second. The panel has hard
  data blockers; 5a has none.
