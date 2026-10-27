#!/usr/bin/env bash
# Publish this repository to GitHub with the GitHub CLI.
#
#   scripts/push_to_github.sh                      # private repo truong128/two-strategies-to-target-kras-mutants (safe default)
#   scripts/push_to_github.sh --public             # public
#   scripts/push_to_github.sh --name my-repo-name  # other name
#
# Steps: audit -> git init -> identity check -> commit -> gh repo create + push -> topics.
# Re-running after a first push only commits new changes and pushes.
set -euo pipefail
cd "$(dirname "$0")/.."

NAME="two-strategies-to-target-kras-mutants"
VIS="--private"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --public) VIS="--public"; shift ;;
    --private) VIS="--private"; shift ;;
    --name) NAME="$2"; shift 2 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done

DESC="From 13.6 billion molecules to a 50-compound shortlist, and from a G>A hotspot to ADAR-recruiting guide RNAs: two reproducible computational routes against KRAS G12D in pancreatic cancer. Nextflow + Python. Code and protocols only, no data."
TOPICS=(kras kras-g12d pancreatic-cancer drug-discovery virtual-screening nextflow adar rna-editing cheminformatics bioinformatics machine-learning molecular-dynamics boltz)

command -v gh >/dev/null || { echo "Install the GitHub CLI: https://cli.github.com" >&2; exit 1; }
gh auth status >/dev/null 2>&1 || { echo "Run: gh auth login" >&2; exit 1; }
OWNER="$(gh api user --jq .login)"

echo "== 1/5 audit"
scripts/check_release.sh

echo "== 2/5 git init"
[[ -d .git ]] || git init -q -b main
git checkout -q -B main 2>/dev/null || true

echo "== 3/5 identity"
if [[ -z "$(git config user.name || true)" || -z "$(git config user.email || true)" ]]; then
  echo "git user.name / user.email are not set. Set them once, then re-run:" >&2
  echo "  git config --global user.name  \"Your Name\"" >&2
  echo "  git config --global user.email \"<your-email>\"" >&2
  exit 1
fi
echo "committing as $(git config user.name) <$(git config user.email)>"

echo "== 4/5 commit"
git add -A
git status --short | head -20
echo "... $(git ls-files | wc -l) files staged"
if ! git diff --cached --quiet; then
  git commit -q -m "Two Therapeutic Strategies to Target KRAS Mutants: end-to-end platform for ultra-large docking, ADAR guide design and peptide binders

Code, configuration, documentation, method figures and tables only. No data or results."
fi

echo "== 5/5 publish to github.com/${OWNER}/${NAME} (${VIS#--})"
if git remote get-url origin >/dev/null 2>&1; then
  git push -u origin main
else
  gh repo create "${OWNER}/${NAME}" "$VIS" --description "$DESC" --source . --remote origin --push
fi
for t in "${TOPICS[@]}"; do gh repo edit "${OWNER}/${NAME}" --add-topic "$t" >/dev/null 2>&1 || true; done
echo
echo "Done: https://github.com/${OWNER}/${NAME}"
