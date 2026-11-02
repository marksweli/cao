#!/usr/bin/env bash
# Generate the whole multi-fidelity dataset with GNU Octave (Linux),
# robust against crashes and hangs.
#
#   ./run_all.sh [N_WORKERS]          (default 2 workers)
#
#   1. builds the LF/HF meshes with Gmsh (if missing)
#   2. starts N_WORKERS Octave processes running
#        run_batch({'config/HF.yaml','config/LF.yaml'})
#      they share the 10 HF + 19 LF cases through per-case lock files;
#   3. a watchdog relaunches a worker whenever it crashes or its log has not
#      changed for STALL_MIN minutes (hang); the killed case resumes from its
#      last checkpoint, finished cases are never recomputed;
#   4. when all 29 cases are done: build_datasets + plot_results.
#
# Environment: REDBKIT_PATH (default ../../redbKIT), OCTAVE (default octave),
#              STALL_MIN (default 20)
# In MATLAB the same is obtained with
#   >> generate_meshes; run_batch({'config/HF.yaml','config/LF.yaml'}); build_datasets; plot_results
# (re-run run_batch after an interruption).
set -u
cd "$(dirname "$0")"
NW=${1:-2}
OCTAVE=${OCTAVE:-octave}
STALL_MIN=${STALL_MIN:-20}
export REDBKIT_PATH=${REDBKIT_PATH:-$(cd ../.. && pwd)/redbKIT}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1   # one core per worker
mkdir -p logs results
N_TOTAL=29

log() { echo "[$(date '+%F %T')] $*" | tee -a logs/watchdog.log; }
n_done() { ls results/LF/Re*/result.mat results/HF/Re*/result.mat 2>/dev/null | wc -l; }

if [ ! -f mesh/LF_coarse.msh ] || [ ! -f mesh/HF_fine.msh ]; then
  log "generating meshes"
  xvfb-run -a $OCTAVE --no-gui --eval "graphics_toolkit qt; generate_meshes" >> logs/mesh.log 2>&1 \
    || $OCTAVE --no-gui --eval "generate_meshes" >> logs/mesh.log 2>&1
fi

start_worker() {   # $1 = worker id
  local w=$1
  echo "===== start $(date '+%F %T') =====" >> logs/worker$w.log
  $OCTAVE --no-gui --eval "run_batch({'config/HF.yaml','config/LF.yaml'})" >> logs/worker$w.log 2>&1 &
  echo $! > logs/worker$w.pid
  log "worker $w started (pid $(cat logs/worker$w.pid)), $(n_done)/$N_TOTAL cases done"
}

for w in $(seq 1 $NW); do start_worker $w; sleep 20; done

while [ "$(n_done)" -lt $N_TOTAL ]; do
  sleep 60
  for w in $(seq 1 $NW); do
    pid=$(cat logs/worker$w.pid)
    if ! kill -0 "$pid" 2>/dev/null; then
      [ "$(n_done)" -ge $N_TOTAL ] && break
      log "worker $w (pid $pid) exited"
      sleep 30; start_worker $w
    elif [ -n "$(find logs/worker$w.log -mmin +$STALL_MIN)" ]; then
      log "worker $w (pid $pid) silent for $STALL_MIN min -> killed and restarted from checkpoint"
      kill -9 "$pid" 2>/dev/null; sleep 5; start_worker $w
    fi
  done
done
for w in $(seq 1 $NW); do wait "$(cat logs/worker$w.pid)" 2>/dev/null; done
log "all $N_TOTAL cases done -> building datasets and figures"
$OCTAVE --no-gui --eval "build_datasets" >> logs/post.log 2>&1
xvfb-run -a $OCTAVE --no-gui --eval "graphics_toolkit qt; plot_results" >> logs/post.log 2>&1 \
  || $OCTAVE --no-gui --eval "plot_results" >> logs/post.log 2>&1
log "finished"
