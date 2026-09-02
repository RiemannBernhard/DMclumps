#!/bin/bash -l
# =============================================================================
# SLURM ARRAY job: draw a statistically significant ENSEMBLE of Milky-Way-type
# galaxies with CLUMPY in -g8 mode (mode g7 "drawn subhaloes 2D" + list haloes).
#
# One array task == one row of clumpy_config/galaxy_ensemble.txt, i.e. one full
# realization of the clumpy Galactic DM halo (different seed / DM mass / clump
# census / clump mass window / line-of-sight orientation). This is the mode and
# workflow used in the paper (DMclumps_and_NS_MNRAS_erratum, skymaps drawn in
# the -g8 mode, N_side = 2048, chi chi -> b bbar) and reproduces the released
# *.drawn catalogues under /data/data/CLUMPY/data/.
#
# Prerequisites on the cluster (one-time):
#   1. Build/copy the image:  apptainer build clumpy.sif clumpy.def
#                             scp clumpy.sif  LOGIN_NODE:~/
#   2. Have the repo (with clumpy_config/) on the cluster, e.g. ~/DMclumps.
#   3. Find your SLURM account:  sacctmgr show assoc -n user=$USER format=Account
#
# Submit (size --array to the number of DATA lines in the manifest):
#   N=$(grep -vcE '^\s*(#|$)' clumpy_config/galaxy_ensemble.txt)
#   sbatch --array=0-$((N-1))%6 slurm/clumpy_g8_ensemble.sbatch
# Override the manifest / image / project:
#   MANIFEST=clumpy_config/my_list.txt SIF=~/clumpy.sif PROJECT=~/DMclumps \
#       sbatch --array=0-30%6 slurm/clumpy_g8_ensemble.sbatch
# =============================================================================

#SBATCH --job-name=clumpy-g8-ens
#SBATCH --array=0-30%6                # 31 rows in galaxy_ensemble.txt; %6 caps concurrency
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8             # OpenMP threads per task
#SBATCH --mem=32gb                    # full-sky nside=2048 maps are memory-hungry
#SBATCH --time=240:00:00              # per task; raise or lower NSIDE as needed
#SBATCH --output=clumpy-g8-%A_%a.out  # %A = array job id, %a = array task id
#SBATCH --error=clumpy-g8-%A_%a.err
##SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=mdeliyer@uni-koeln.de

set -euo pipefail

# ---- Inputs / outputs (override via environment) ---------------------------
SIF=${SIF:-$HOME/TestArea/DMclumps/clumpy.sif}
PROJECT=${PROJECT:-$HOME/TestArea/DMclumps}
MANIFEST=${MANIFEST:-clumpy_config/galaxy_ensemble.txt}
PARAMS=${PARAMS:-clumpy_config/clumpy_params_g8.txt}

[[ -f "$SIF" ]]                 || { echo "ERROR: SIF not found: $SIF" >&2; exit 1; }
[[ -f "$PROJECT/$MANIFEST" ]]   || { echo "ERROR: manifest not found: $PROJECT/$MANIFEST" >&2; exit 1; }
[[ -f "$PROJECT/$PARAMS" ]]     || { echo "ERROR: param file not found: $PROJECT/$PARAMS" >&2; exit 1; }

# ---- Pick the manifest row for this array task -----------------------------
# Strip comments/blank lines, then take the (task_id+1)-th remaining line.
i=${SLURM_ARRAY_TASK_ID:-0}
ROW=$(grep -vE '^\s*(#|$)' "$PROJECT/$MANIFEST" | sed -n "$((i + 1))p")
if [[ -z "$ROW" ]]; then
    echo "Array index $i has no matching data line in $MANIFEST — nothing to do." >&2
    echo "Set --array=0-\$((N-1)) with N = $(grep -vcE '^\s*(#|$)' "$PROJECT/$MANIFEST")." >&2
    exit 0
fi
read -r ID SEED MDM NCLUMPS MMIN MMAXFR PSI THETA NSIDE GALLIST LABEL <<< "$ROW"

# ---- Apptainer -------------------------------------------------------------
RUNNER=$(command -v apptainer || command -v singularity) || {
    echo "ERROR: apptainer/singularity not found. Run 'module avail apptainer'." >&2
    exit 1
}

WORKDIR=/scratch/${USER}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$WORKDIR"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

# ---- Build the optional per-row overrides (only when != -1) ----------------
OPTS=()
[[ "$NCLUMPS" != "-1" ]] && OPTS+=( "--gMW_SUBS_N_INM1M2=$NCLUMPS" )
[[ "$MMIN"    != "-1" ]] && OPTS+=( "--gDM_SUBS_MMIN=$MMIN" )
[[ "$MMAXFR"  != "-1" ]] && OPTS+=( "--gDM_SUBS_MMAXFRAC=$MMAXFR" )

# GALLIST: which external MW-type galaxy (inclination) to inject via gLIST_HALOES.
#   default -> keep the param-file list;  none/-1 -> add nothing;  <path> -> that list file.
case "$GALLIST" in
    default|DEFAULT|"" ) : ;;                                       # keep clumpy_params_g8.txt default
    none|NONE|-1 )       OPTS+=( "--gLIST_HALOES=-1" ) ;;
    * )                  OPTS+=( "--gLIST_HALOES=$GALLIST" )
                         [[ -f "$PROJECT/$GALLIST" ]] \
                             || echo "WARN: list-halo file not found: $PROJECT/$GALLIST" >&2 ;;
esac

echo "Host    : $(hostname)"
echo "Array   : ${SLURM_ARRAY_JOB_ID:-?} task ${SLURM_ARRAY_TASK_ID:-?}"
echo "Galaxy  : $ID ($LABEL)"
echo "Params  : seed=$SEED mDM=${MDM}GeV nClumps=$NCLUMPS Mmin=$MMIN MmaxFrac=$MMAXFR"
echo "Orient. : LOS psi=$PSI theta=$THETA  nside=$NSIDE  threads=$OMP_NUM_THREADS"
echo "ListHalo: $GALLIST"
echo "Scratch : $WORKDIR"

# ---- Run CLUMPY -g8 in the container ---------------------------------------
# Full-sky 360 deg => RING mode: --gSIM_THETA_ORTH_SIZE_DEG=d
"$RUNNER" exec \
    --bind "$PROJECT":/work,"$WORKDIR":/out \
    --pwd /work \
    "$SIF" \
    clumpy -g8 -i "$PARAMS" \
        --gSIM_SEED="$SEED" --gPP_DM_MASS_GEV="$MDM" --gSIM_HEALPIX_NSIDE="$NSIDE" \
        --gSIM_PSI_OBS_DEG="$PSI" --gSIM_THETA_OBS_DEG="$THETA" \
        --gSIM_THETA_ORTH_SIZE_DEG=d --gSIM_THETA_SIZE_DEG=360 \
        ${OPTS[@]+"${OPTS[@]}"} \
        --gSIM_OUTPUT_DIR=/out

# ---- Stage results back, one folder per galaxy realization -----------------
DEST="$PROJECT/output/${ID}_seed${SEED}_mDM${MDM}GeV_LOS${PSI}_${THETA}_nside${NSIDE}"
mkdir -p "$DEST"
cp -av "$WORKDIR"/*.drawn "$WORKDIR"/*.fits "$WORKDIR"/*.root "$DEST"/ 2>/dev/null || true
echo "Done. Catalogue(s) for $ID in: $DEST"
