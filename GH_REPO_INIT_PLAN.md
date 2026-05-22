# GitHub repo initialization plan

When Ian gives the green light ("push it"), run these commands in order
from `/Users/ian/1A_Helfrich_ThesisResearch_May2024/Paper5_EffectiveDistance/`.

## 1. Initialize local repo

```bash
cd /Users/ian/1A_Helfrich_ThesisResearch_May2024/Paper5_EffectiveDistance/
git init
git config user.email "ianthelfrich@gmail.com"
git config user.name "Ian Helfrich"
```

## 2. Stage files

```bash
# Stage everything EXCEPT what's in .gitignore
git add .gitignore README.md manuscript/ src/ tests/ notebooks/ \
  learning_guide/ docs/ pyproject.toml refs.bib Cargo.toml \
  crates/ workflow/ AUTONOMOUS_BLOCK_DELIVERABLES.md \
  ARCHITECTURE.md DATA_INVENTORY.md METHODOLOGY_NOTES.md \
  NUMERICAL_ISSUES.md PRE_ANALYSIS_PLAN.md PROJECT_BRIEF.md \
  RESEARCH_AGENDA.md SCOPING_DOC.md SPRINT_STATUS.md

# Check what's about to be committed
git status

# Inspect the diff before commit
git diff --cached --stat
```

## 3. Initial commit

```bash
git commit -m "Initial commit: Paper 5a (the wrong-theta critique) + companion materials

Includes:
- Paper 5a v0.5 (13-page SSRN preprint, manuscript/paper_5a/)
- Visual abstract (3-page executive summary)
- Volume 1 (35-page) and Volume 2 (50-page) of teaching guide
- Bug-fixed gravity regression notebooks (18d, 18e)
- Spec B-prime regularization sweep data (96 regressions)
- Pre-emptive referee response document
- SSRN posting metadata
- Bibliography: 34 verified entries
- Literature triangulation reports
- GitHub Pages stub (docs/)

Status: SSRN posting target within 2 weeks. PPML robustness check
running but with convergence warnings; ppmlhdfe-style implementation
is the natural follow-up.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

## 4. Create GitHub remote

```bash
gh repo create ihelfrich/paper5-effective-distance \
  --public \
  --description "Paper 5a: The Untaken Limit — A Pole at theta = -2 in the Head-Mayer Closed Form for Intra-National Distance. SSRN preprint + replication code + teaching guide." \
  --homepage "https://ihelfrich.github.io/paper5-effective-distance/" \
  --source . \
  --remote origin \
  --push
```

The `--push` flag will push the initial commit immediately. If you'd
rather verify first, drop `--push` and run `git push -u origin main`
separately.

## 5. Enable GitHub Pages

```bash
# Set the docs/ directory as the Pages source
gh api -X POST /repos/ihelfrich/paper5-effective-distance/pages \
  -f source.branch=main \
  -f source.path=/docs
```

Wait ~60 seconds for Pages to propagate. Verify at:
`https://ihelfrich.github.io/paper5-effective-distance/`

## 6. Verify the live site

```bash
# Wait for Pages to propagate
sleep 60

# Open in browser (macOS)
open "https://ihelfrich.github.io/paper5-effective-distance/"

# Verify the PDF link works
curl -I "https://ihelfrich.github.io/paper5-effective-distance/paper_5a.pdf"
```

## 7. Tag the initial release

```bash
git tag -a v0.5.0 -m "v0.5: pre-SSRN posting preprint with PPML preliminary results"
git push origin v0.5.0

# Optional: create a GitHub release with the PDF as a binary asset
gh release create v0.5.0 \
  manuscript/paper_5a/paper_5a.pdf \
  manuscript/paper_5a/visual_abstract.pdf \
  learning_guide/master_guide_vol2.pdf \
  --title "v0.5: Pre-SSRN preprint" \
  --notes "13-page SSRN preprint with PPML preliminary results. Visual abstract and Volume 2 teaching guide attached. Awaiting SSRN posting."
```

## 8. After SSRN posting

Once Ian assigns the SSRN ID:

```bash
# Edit paper_5a.tex to add the SSRN URL in the title footnote
# Recompile
xelatex paper_5a.tex && bibtex paper_5a && xelatex paper_5a.tex && xelatex paper_5a.tex

# Edit docs/index.html to link to the SSRN URL
# Commit and push
git add manuscript/paper_5a/paper_5a.tex manuscript/paper_5a/paper_5a.pdf docs/
git commit -m "Add SSRN ID to paper and Pages"
git push

# Tag the SSRN-posted version
git tag -a v1.0.0 -m "v1.0: SSRN preprint posted"
git push origin v1.0.0
```

## 9. Zenodo DOI

```bash
# Connect the repo to Zenodo (one-time, via Zenodo web UI):
# https://zenodo.org/account/settings/github/
# Toggle on for the paper5-effective-distance repo.
# Then any GitHub release automatically gets a Zenodo DOI.

# After v1.0.0 release, fetch the Zenodo DOI and add to:
#   - manuscript/paper_5a/paper_5a.tex (replace [DOI to be assigned])
#   - docs/index.html
#   - README.md
```

## Verification checklist before pushing

- [ ] `.gitignore` is excluding `.venv/`, `__pycache__/`, raw data files
- [ ] `git status` shows expected files staged
- [ ] No API keys, Zotero local paths, or PII in committed files
- [ ] LaTeX `paper_5a.tex` does not contain "[acknowledgments to be added]"
      placeholder text (or it's acceptable to leave for v0.5)
- [ ] LaTeX `paper_5a.tex` does not contain "[DOI to be assigned]" in a
      way that would break (it's OK as a placeholder; will be fixed at v1.0)
- [ ] README.md has the right contact email and affiliation
- [ ] docs/index.html links work (or will work once the repo exists)

## Rollback plan if something goes wrong

```bash
# Delete the remote (irreversible)
gh repo delete ihelfrich/paper5-effective-distance --yes

# Reset the local repo
rm -rf .git
```

---

*Document prepared 2026-05-21 during the autonomous block. Awaiting
Ian's "push it" go-ahead. Per CLAUDE.md §3.6: "Commit or push only when
the user asks."*
