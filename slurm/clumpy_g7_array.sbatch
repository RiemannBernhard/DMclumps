#!/bin/bash -l
# =============================================================================
# Generic SLURM ARRAY job: sweep a CLUMPY -g7 draw over a (seed x DM-mass) grid,
# one array task per combination, running the clumpy.sif Apptainer image.
#
# Default grid: 9 seeds x {100, 500} GeV = 18 tasks (--array=0-17).
# Each task writes to its own scratch dir and stages results back to a per-task
# subfolder, so tasks never collide.
#
# Prerequisites on the cluster (one-time):
#   1. Build/copy the image:  apptainer build clumpy.sif clumpy.def
#                             scp clumpy.sif  LOGIN_NODE:~/
#   2. Have the repo (with clumpy_config/) on the cluster, e.g. ~/DMclumps.
#   3. Find your SLURM account:  sacctmgr show assoc -n user=$USER format=Account
#
# Submit:   sbatch slurm/clumpy_g7_array.sbatch
# Override: SEEDS="4448 4450" MASSES="100" NSIDE=512 sbatch slurm/clumpy_g7_array.sbatch
#           (adjust --array=0-<N-1> to match len(SEEDS)*len(MASSES))
# =============================================================================

#SBATCH --job-name=clumpy-g7-array
#SBATCH --array=0-17%6                # 18 tasks; %6 caps concurrent tasks at 6
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8             # OpenMP threads per task
#SBATCH --mem=32gb
#SBATCH --time=12:00:00              # per task; full-sky nside=2048 is slow
#SBATCH --account=CHANGE_ME          # <-- sacctmgr show assoc -n user=$USER format=Account
#SBATCH --output=clumpy-%A_%a.out    # %A = array job id, %a = array task id
#SBATCH --error=clumpy-%A_%a.err
##SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=USERNAME@example.org

set -euo pipefail

# ---- Grid definition (override via environment) ----------------------------
read -r -a SEEDS  <<< "${SEEDS:-4448 4450 4452 4454 4456 4457 4459 4462 4465}"
read -r -a MASSES <<< "${MASSES:-100 500}"
NSIDE=${NSIDE:-2048}

NM=${#MASSES[@]}
NTOTAL=$(( ${#SEEDS[@]} * NM ))

# Map this array task id -> (seed, mass)
i=${SLURM_ARRAY_TASK_ID:-0}
if (( i >= NTOTAL )); then
    echo "Array index $i >= grid size $NTOTAL — nothing to do (fix --array=0-$((NTOTAL-1)))." >&2
    exit 0
fi
SEED=${SEEDS[$(( i / NM ))]}
MDM=${MASSES[$(( i % NM ))]}

# ---- Apptainer -------------------------------------------------------------
module load apptainer 2>/dev/null \
  || module load tools/apptainer 2>/dev/null \
  || module load system/apptainer 2>/dev/null \
  || module load singularity 2>/dev/null \
  || true
RUNNER=$(command -v apptainer || command -v singularity) || {
    echo "ERROR: apptainer/singularity not found. Run 'module avail apptainer' and fix above." >&2
    exit 1
}

# ---- Inputs / outputs ------------------------------------------------------
SIF=${SIF:-$HOME/clumpy.sif}
PROJECT=${PROJECT:-$HOME/DMclumps}
WORKDIR=/scratch/${USER}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$WORKDIR"

[[ -f "$SIF" ]] || { echo "ERROR: SIF not found: $SIF" >&2; exit 1; }
[[ -f "$PROJECT/clumpy_config/clumpy_params_g7.txt" ]] \
    || { echo "ERROR: param file not found under PROJECT=$PROJECT" >&2; exit 1; }

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

echo "Host    : $(hostname)"
echo "Array   : ${SLURM_ARRAY_JOB_ID:-?} task ${SLURM_ARRAY_TASK_ID:-?} of grid size $NTOTAL"
echo "Combo   : seed=$SEED  mDM=${MDM}GeV  nside=$NSIDE  threads=$OMP_NUM_THREADS"
echo "Scratch : $WORKDIR"

# ---- Run CLUMPY in the container -------------------------------------------
# Full-sky 360 deg => RING mode: --gSIM_THETA_ORTH_SIZE_DEG=d
"$RUNNER" exec \
    --bind "$PROJECT":/work,"$WORKDIR":/out \
    --pwd /work \
    "$SIF" \
    clumpy -g7 -i clumpy_config/clumpy_params_g7.txt \
        --gSIM_SEED="$SEED" --gPP_DM_MASS_GEV="$MDM" --gSIM_HEALPIX_NSIDE="$NSIDE" \
        --gSIM_THETA_ORTH_SIZE_DEG=d --gSIM_THETA_SIZE_DEG=360 \
        --gSIM_OUTPUT_DIR=/out

# ---- Stage results back, one folder per (seed, mass) -----------------------
DEST="$PROJECT/output/seed${SEED}_mDM${MDM}GeV"
mkdir -p "$DEST"
cp -av "$WORKDIR"/*.drawn "$WORKDIR"/*.fits "$WORKDIR"/*.root "$DEST"/ 2>/dev/null || true
echo "Done. Catalog(s) in: $DEST"
