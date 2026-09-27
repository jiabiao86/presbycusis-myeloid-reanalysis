#!/usr/bin/env bash
# Drive the full CellChatDB.mouse run (all 2,019 interactions) for ages 3/12/24.
# Chunks are merged afterwards with cellchat_full_merge.R.
set -u
cd "$(dirname "$0")"

CHUNK_SIZE=${CHUNK_SIZE:-505}
NBOOT=${NBOOT:-100}
CONCURRENCY=${CONCURRENCY:-4}
N_INTERACTIONS=2019

mkdir -p cellchat_full/logs cellchat_full/chunks

for age in 3 12 24; do
  start=1
  while [ "$start" -le "$N_INTERACTIONS" ]; do
    end=$((start + CHUNK_SIZE - 1))
    if [ "$end" -gt "$N_INTERACTIONS" ]; then
      end=$N_INTERACTIONS
    fi
    while [ "$(jobs -rp | wc -l | tr -d ' ')" -ge "$CONCURRENCY" ]; do
      sleep 5
    done
    Rscript cellchat_full_chunk.R "$age" "$start" "$end" "$NBOOT" \
      > "cellchat_full/logs/chunk_${age}M_${start}-${end}.log" 2>&1 &
    start=$((end + 1))
  done
done

wait
echo "ALL CHUNKS DONE"
