# Publish to GitHub

## Repository identity

| | |
|---|---|
| **Name** | `two-strategies-to-target-kras-mutants` |
| **Title** | Two Therapeutic Strategies to Target KRAS Mutants |
| **Tagline** | An open, end-to-end computational platform that screens 13.6 billion molecules and designs ADAR RNA editors against the same pancreatic cancer driver |
| **GitHub "About" text** | From 13.6 billion molecules to a 50-compound shortlist, and from a G>A hotspot to ADAR-recruiting guide RNAs: two reproducible computational routes against KRAS G12D in pancreatic cancer. Nextflow + Python. Code and protocols only, no data. |
| **Topics** | kras, kras-g12d, pancreatic-cancer, drug-discovery, virtual-screening, nextflow, adar, rna-editing, cheminformatics, bioinformatics, machine-learning, molecular-dynamics, boltz |

## Option 1. One command (recommended)

From the repository root:

```bash
gh auth status                      # must show your account
scripts/push_to_github.sh           # private first; add --public when you are ready
```

The script audits the files, initializes git, commits, creates `github.com/<you>/two-strategies-to-target-kras-mutants` with the description and topics above, and pushes.

## Option 2. Manual commands

```bash
cd <this folder>

# 1. audit: refuses data, results, secrets, personal paths, em/en dashes, large files
scripts/check_release.sh

# 2. what exactly will be published? (largest files last)
git init -b main
git add -A --dry-run | wc -l
git ls-files -co --exclude-standard | xargs du -k | sort -n | tail -5

# 3. commit
# git identity must already be set (git config --global user.name / user.email)
git add -A
git commit -m "Two Therapeutic Strategies to Target KRAS Mutants: end-to-end platform"

# 4. create the repo and push
gh repo create truong128/two-strategies-to-target-kras-mutants --private \
  --description "From 13.6 billion molecules to a 50-compound shortlist, and from a G>A hotspot to ADAR-recruiting guide RNAs: two reproducible computational routes against KRAS G12D in pancreatic cancer. Nextflow + Python. Code and protocols only, no data." \
  --source . --remote origin --push

# 5. topics, then make it public when ready
for t in kras kras-g12d pancreatic-cancer drug-discovery virtual-screening nextflow adar rna-editing cheminformatics bioinformatics machine-learning molecular-dynamics boltz; do
  gh repo edit truong128/two-strategies-to-target-kras-mutants --add-topic "$t"; done
gh repo edit truong128/two-strategies-to-target-kras-mutants --visibility public --accept-visibility-change-consequences
```

Without the GitHub CLI: create an empty repository in the web UI, then

```bash
git remote add origin https://github.com/truong128/two-strategies-to-target-kras-mutants.git
git push -u origin main
```

## Release polish (optional)

```bash
git tag -a v1.0.0 -m "First public release"
git push origin v1.0.0
gh release create v1.0.0 --title "v1.0.0" --notes "First public release: Tracks A, A2, B, C. Code and protocols only."
```

## Before you flip it public

1. Read `README.md` once as a stranger would; the repository has no data, so every claim must stand on method alone.
2. Confirm you are allowed to publish the third-party-derived code (RESTORE 2.0 handling, Enamine-related scripts); no vendor data or publisher files are included.
3. Decide on the license (MIT by default) and your author line in `LICENSE` and `CITATION.cff`.
4. Delete nothing from `.gitignore`.
