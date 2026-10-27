#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p external
if [[ -d external/AdarEdit/.git ]]; then
  echo 'Existing AdarEdit checkout retained (pinned by local git revision).'
else
  git clone https://github.com/Scientific-Computing-Lab/AdarEdit.git external/AdarEdit
fi
(cd external/AdarEdit && git rev-parse HEAD) > external/AdarEdit.commit
printf 'Repo checkout: %s\n' "$(cat external/AdarEdit.commit)"
echo 'Install its Linux x86_64 PyTorch/PyG requirements in a separate environment; see external/AdarEdit/README.md.'
echo 'Do NOT label its predicted native-site editability as modified-guide editing efficiency.'
