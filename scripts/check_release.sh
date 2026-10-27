#!/usr/bin/env bash
# Pre-push audit. Fails (exit 1) if the files git would publish contain data, results, secrets,
# personal absolute paths, em or en dashes, or anything unexpectedly large.
#
#   scripts/check_release.sh
set -euo pipefail
cd "$(dirname "$0")/.."
[[ -d .git ]] || git init -q
MAX_FILE_KB="${MAX_FILE_KB:-1024}"
MAX_TOTAL_KB="${MAX_TOTAL_KB:-8192}"
fail=0
bad() { printf '\033[1;31mFAIL\033[0m %s\n' "$*"; fail=1; }
ok()  { printf '\033[1;32mOK\033[0m   %s\n' "$*"; }

mapfile -t FILES < <(git ls-files -co --exclude-standard | sort -u)
echo "files that would be published: ${#FILES[@]}"
[[ ${#FILES[@]} -gt 0 ]] || { bad "nothing to publish"; exit 1; }

# 1. size
total=0
for f in "${FILES[@]}"; do
  [[ -f "$f" ]] || continue
  kb=$(( ($(stat -c %s "$f") + 1023) / 1024 )); total=$((total + kb))
  (( kb > MAX_FILE_KB )) && bad "large file (${kb} KB): $f"
done
(( total > MAX_TOTAL_KB )) && bad "total size ${total} KB exceeds ${MAX_TOTAL_KB} KB" || ok "total size ${total} KB"

# 2. forbidden data-like names
forbidden='\.(pdb|pdbqt|cif|mmcif|sdf|mol2|prmtop|inpcrd|rst7|dcd|nc|xtc|trr|bam|bai|cram|fastq|fq|vcf|bcf|xlsx|xls|zip|tar|tgz|7z|ckpt|pt|pth|safetensors|joblib|pkl|pickle|h5|onnx|gtf|bz2|cxsmiles)(\.gz)?$|(^|/)(results|runs|work|data|structures|reference|source_material|targets|external)/'
allow='^tracks/b_adar_rna_editing/examples/synthetic_only/'
hits=0
for f in "${FILES[@]}"; do
  if [[ "$f" =~ $forbidden ]] && ! [[ "$f" =~ $allow ]]; then bad "data-like file: $f"; hits=1; fi
done
(( hits == 0 )) && ok "no data-like files"

# 3. content scan (text files only, skip this script)
text=()
for f in "${FILES[@]}"; do
  [[ "$f" == "scripts/check_release.sh" || "$f" == *.png || "$f" == *.svg ]] && continue
  [[ -f "$f" ]] && grep -Iq . "$f" 2>/dev/null && text+=("$f")
done
scan() {  # scan <label> <regex> [extra grep flags]
  local label="$1" re="$2" out
  out="$(grep -nE "$re" "${text[@]}" 2>/dev/null | grep -v 'users\.noreply\.github\.com' | head -5 || true)"
  if [[ -n "$out" ]]; then bad "$label"; printf '%s\n' "$out" | sed 's/^/       /'; else ok "$label: none"; fi
}
scan "personal absolute paths" '(/home/[A-Za-z0-9_.-]+/|/Users/[A-Za-z0-9_.-]+/|/media/[A-Za-z0-9_.-]+/|[A-Z]:\\Users\\)'
scan "tokens and keys" '(gh[pousr]_[A-Za-z0-9]{20,}|hf_[A-Za-z0-9]{20,}|AKIA[0-9A-Z]{16}|BEGIN [A-Z ]*PRIVATE KEY|api[_-]?key *[:=] *[A-Za-z0-9]{16,})'
scan "em or en dashes" $'(\xe2\x80\x94|\xe2\x80\x93)'
scan "PDB coordinate records" '^(ATOM  |HETATM)'
scan "email addresses" '[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}'

if (( fail )); then echo; echo "Release check FAILED. Fix the lines above before pushing."; exit 1; fi
echo; echo "Release check PASSED: ${#FILES[@]} files, ${total} KB."
