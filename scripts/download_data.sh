#!/usr/bin/env bash
# Download the 39 participants of the final sample from OpenNeuro (ds006979)
# into data/ds006979 (about 25 GB). Files already downloaded are skipped, so the
# script can be rerun after an interruption.
#
# Use dataset version 1.1.0 or later: in versions 1.0.x, the files of sub-10 and
# sub-49 were not the recordings used in the analysis (README, "sub-10 and sub-49").
#
# Usage (from the project folder):
#   bash scripts/download_data.sh

cd "$(dirname "$0")/.." || exit 1
BASE=https://s3.amazonaws.com/openneuro.org/ds006979
for s in 01 02 07 08 09 10 11 12 13 14 16 17 18 19 21 22 23 24 28 30 32 34 35 40 \
         47 48 49 50 51 54 58 59 61 62 64 69 70 73 76; do
  d=data/ds006979/sub-$s/eeg; mkdir -p "$d"
  for f in eeg.vhdr eeg.vmrk eeg.eeg events.tsv; do
    out=$d/sub-${s}_task-visualworkingmemorytask_$f
    [ -s "$out" ] && [ "$(wc -c < "$out")" -gt 1000 ] && continue
    curl -sSfL --retry 3 -o "$out" "$BASE/sub-$s/eeg/sub-${s}_task-visualworkingmemorytask_$f" \
      || echo "  failed: $out (rerun the script)"
  done
  echo "sub-$s done"
done

echo "Participants: $(ls -d data/ds006979/sub-* | wc -l) (expected 39)"
