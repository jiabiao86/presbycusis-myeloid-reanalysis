#!/usr/bin/env bash
# Merge the chunked full CellChatDB.mouse run and finalise one object per age.
set -u
cd "$(dirname "$0")"

NBOOT=${NBOOT:-100}
mkdir -p cellchat_full/logs cellchat_full/final

for age in 3 12 24; do
  Rscript cellchat_full_merge.R "$age" "$NBOOT" \
    > "cellchat_full/logs/merge_${age}M.log" 2>&1
  echo "merged ${age}M"
done
echo "ALL MERGES DONE"
