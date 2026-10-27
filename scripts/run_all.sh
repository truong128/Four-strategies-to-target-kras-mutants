#!/usr/bin/env bash
# One entry point for the whole platform.
#
#   scripts/run_all.sh smoke              # tests for every track that needs no data or GPU
#   scripts/run_all.sh setup              # bootstrap envs + download public inputs
#   scripts/run_all.sh track-b            # ADAR guide design (plan only)
#   scripts/run_all.sh track-a LIB=PATH   # ultra-large docking campaign on the library at PATH
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
case "${1:-}" in
  smoke)   "$HERE/run_track_a.sh" smoke; "$HERE/run_track_b.sh" smoke ;;
  setup)   "$HERE/bootstrap.sh" all; "$HERE/fetch_public_inputs.sh" all ;;
  track-a) shift; LIB="${1#LIB=}"; "$HERE/run_track_a.sh" screen --library "$LIB" ;;
  track-a2) shift; "$HERE/run_track_a2.sh" "$@" ;;
  track-b) "$HERE/run_track_b.sh" preflight; "$HERE/run_track_b.sh" plan ;;
  track-c) shift; "$HERE/run_track_c.sh" "${1:-prepare}" ;;
  *) echo "usage: run_all.sh smoke | setup | track-a LIB=PATH | track-a2 | track-b | track-c [prepare|design]" >&2; exit 2 ;;
esac
