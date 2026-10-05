#!/bin/bash
# NOTE (sanitasi Sep 2026): path absolut asal (home user asal,
# /mgpfs/scratch/$USER) diganti $HOME dan /mgpfs/scratch/$USER.
# Sesuaikan PROJECT_DIR / tool paths di bawah dengan akun HPC masing-masing.
set -e
export PATH=$HOME/.conda/envs/training_qc/bin:$PATH
export LD_LIBRARY_PATH=$HOME/.conda/envs/training_qc/lib:$LD_LIBRARY_PATH
NP=$HOME/acacia_project/tools/NextPolish/NextPolish/nextPolish

BASE=$HOME/polishrun
ROUNDS=${1:-5}
GENOME=$BASE/assembly/nd.asm.fasta
SGS=$BASE/sgs.fofn

mkdir -p "$BASE/rounds" "$BASE/cfgs"

for ((r=1; r<=ROUNDS; r++)); do
  OUT="$BASE/rounds/round_$r"
  mkdir -p "$OUT"
  CFG="$BASE/cfgs/round_$r.cfg"
  printf '[General]\njob_type = local\njob_prefix = nextPolish_r%s\ntask = best\nrewrite = yes\ndeltmp = yes\nrerun = 5\nparallel_jobs = 6\nmultithread_jobs = 5\ngenome = %s\ngenome_size = auto\nworkdir = %s/01_rundir\npolish_options = -p {multithread_jobs}\n\n[sgs_option]\nsgs_fofn = %s\nsgs_options = -max_depth 100 -minimap2\n' "$r" "$GENOME" "$OUT" "$SGS" > "$CFG"

  echo "=== ROUND $r/$ROUNDS : input=$GENOME ==="
  "$NP" "$CFG"
  OUTFA=$(ls "$OUT/01_rundir/"*.nextpolish.fasta 2>/dev/null | head -1)
  if [ -z "$OUTFA" ]; then
    echo "ROUND $r FAILED : no output fasta in $OUT/01_rundir/"
    ls "$OUT/01_rundir/" 2>/dev/null | head -20
    exit 1
  fi
  GENOME="$OUTFA"
done
echo "=== ALL $ROUNDS ROUNDS DONE ==="
echo "FINAL: $GENOME"