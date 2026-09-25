#!/usr/bin/env bash
# Validate every workflow entry point (phase*/terra/wdl/) with womtool (Cromwell's validator,
# the engine Terra uses).
#
# Usage: scripts/validate_wdls.sh [womtool.jar]
# Without an argument, uses a `womtool` on PATH.
set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
cd "$root"

if [ $# -ge 1 ]; then
  womtool=(java -jar "$1")
else
  womtool=(womtool)
fi

failed=0
for wdl in phase*/terra/wdl/*/*/*.wdl; do
  if "${womtool[@]}" validate "$wdl" > /dev/null 2>&1; then
    echo "ok    $wdl"
  else
    echo "FAIL  $wdl"
    "${womtool[@]}" validate "$wdl" 2>&1 | head -5 | sed 's/^/      /'
    failed=$((failed + 1))
  fi
done

echo "$failed failed"
[ "$failed" -eq 0 ]
