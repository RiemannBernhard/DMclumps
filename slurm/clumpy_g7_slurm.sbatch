#!/bin/bash -l
# =============================================================================
# Generic SLURM HPC batch job: run the CLUMPY Apptainer container.
#
# CLUMPY is a single-node, CPU job (OpenMP, not MPI), so this is an SMP job:
# 1 task, several threads. With no --partition it auto-routes to `smp`.
#
# Prerequisites on the cluster (one-time):
#   1. Build the image on a host with Docker + Apptainer, then copy it over:
#        apptainer build clumpy.sif clumpy.def
#        scp clumpy.sif  LOGIN_NODE:~/            # or build from a registry/tarball
#   2. Have the repo (with clumpy_config/) on the cluster, e.g. ~/DMclumps.
#   3. Find your SLURM account:  sacctmgr show assoc -n user=$USER format=Account
#      and put it in --account below (jobs are billed to it).
#
# Submit:   sbatch slurm/clumpy_g7_slurm.sbatch
# Override: SIF=~/clumpy.sif PROJECT=~/DMclumps NSIDE=512 SEED=4448 \
#             sbatch slurm/clumpy_g7_slurm.sbatch
# =============================================================================

#SBATCH --job-name=clumpy-g7
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8           # OpenMP threads for CLUMPY/ROOT
#SBATCH --mem=32gb                  # nside=2048 full-sky maps are memory-hungry
#SBATCH --time=12:00:00             # full-sky J-map is slow; raise or lower NSIDE
#SBATCH --account=CHANGE_ME         # <-- sacctmgr show assoc -n user=$USER format=Account
#SBATCH --output=clumpy-%j.out
#SBATCH --error=clumpy-%j.err
# Email is optional; some sites restrict it to institutional addresses:
##SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=USERNAME@example.org
# For a GPU node you would add (CLUMPY does not need it):  --partition=gpu -G h100:1

set -euo pipefail

# ---- Apptainer -------------------------------------------------------------
# The cluster may expose Apptainer as a module or have it on PATH; try a few names.
module load apptainer 2>/dev/null \
  || module load tools/apptainer 2>/dev/null \
  || module load system/apptainer 2>/dev/null \
  || module load singularity 2>/dev/null \
  || true
command -v apptainer >/dev/null 2>&1 || command -v singularity >/dev/null 2>&1 || {
    echo "ERROR: apptainer/singularity not found. Run 'module avail apptainer' on the cluster" >&2
    echo "       and set the correct 'module load ...' line above." >&2
    exit 1
}
RUNNER=$(command -v apptainer || command -v singularity)

# ---- Inputs / outputs (override via environment) ---------------------------
SIF=${SIF:-$HOME/clumpy.sif}                 # the image copied to the cluster
PROJECT=${PROJECT:-$HOME/DMclumps}           # repo containing clumpy_config/
WORKDIR=/scratch/${USER}/${SLURM_JOB_ID}     # fast scratch for the run
mkdir -p "$WORKDIR"

# CLUMPY draw parameters (override via env)
SEED=${SEED:-4448}
MDM=${MDM:-100}
NSIDE=${NSIDE:-2048}

[[ -f "$SIF" ]]                       || { echo "ERROR: SIF not found: $SIF" >&2; exit 1; }
[[ -f "$PROJECT/clumpy_config/clumpy_params_g7.txt" ]] \
    || { echo "ERROR: param file not found under PROJECT=$PROJECT" >&2; exit 1; }

# Threads = cores SLURM allocated
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

echo "Host        : $(hostname)"
echo "Account/Job : ${SLURM_JOB_ACCOUNT:-?} / ${SLURM_JOB_ID}"
echo "Image       : $SIF"
echo "Project     : $PROJECT  (mounted at /work)"
echo "Scratch out : $WORKDIR  (mounted at /out)"
echo "Params      : seed=$SEED  mDM=${MDM}GeV  nside=$NSIDE  threads=$OMP_NUM_THREADS"

# ---- Run CLUMPY in the container -------------------------------------------
# Bind the repo at /work (read inputs) and scratch at /out (write results).
# --pwd /work sets the container CWD and PWD (CLUMPY reads getenv("PWD")).
# Full-sky 360 deg map => RING mode: --gSIM_THETA_ORTH_SIZE_DEG=d
"$RUNNER" exec \
    --bind "$PROJECT":/work,"$WORKDIR":/out \
    --pwd /work \
    "$SIF" \
    clumpy -g7 -i clumpy_config/clumpy_params_g7.txt \
        --gSIM_SEED="$SEED" --gPP_DM_MASS_GEV="$MDM" --gSIM_HEALPIX_NSIDE="$NSIDE" \
        --gSIM_THETA_ORTH_SIZE_DEG=d --gSIM_THETA_SIZE_DEG=360 \
        --gSIM_OUTPUT_DIR=/out

# ---- Stage results back to the project --------------------------------------
DEST="$PROJECT/output"
mkdir -p "$DEST"
cp -av "$WORKDIR"/*.drawn "$WORKDIR"/*.fits "$WORKDIR"/*.root "$DEST"/ 2>/dev/null || true
echo "Done. Catalog(s) in: $DEST"
