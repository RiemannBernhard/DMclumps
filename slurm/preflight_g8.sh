#!/usr/bin/env bash
# =============================================================================
# Pre-flight validator for the CLUMPY -g8 ensemble array job.
#
# Catches, before you burn cluster time, the failure modes that have actually
# occurred in this project:
#   * duplicate realizations  -- rows that differ only in a parameter which does
#     NOT change the clump draw (DM mass, line of sight, injected galaxy list),
#     so they produce byte-identical .drawn catalogues. The previous g7 grid
#     shipped 15 files that held only 9 distinct realizations.
#   * heterogeneous ensembles -- rows with different NCLUMPS / MMIN / MMAXFRAC
#     pooled as if they were one statistical population.
#   * missing param or list-halo files -- the array dies instantly on the first,
#     and silently changes physics on the second (the driver only warns).
#   * --array range not covering the manifest.
#   * a --time request the partition will refuse at submit.
#
# Usage
#   slurm/preflight_g8.sh [-m MANIFEST] [-p PARAMS] [-s SBATCH] [-P PROJECT] [-i SIF]
#   slurm/preflight_g8.sh --post OUTPUT_DIR     # after the run: checksum-dedup .drawn
#
# Exit status: 0 = ready to submit, 1 = blocking error(s) found.
# =============================================================================
set -uo pipefail

PROJECT=${PROJECT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
MANIFEST=${MANIFEST:-clumpy_config/galaxy_ensemble_100.txt}
PARAMS=${PARAMS:-clumpy_config/clumpy_params_g8.txt}
SBATCH=${SBATCH:-slurm/clumpy_g8_ensemble.sbatch}
SIF=${SIF:-$HOME/TestArea/DMclumps/clumpy.sif}
GB_PER_REALIZATION=${GB_PER_REALIZATION:-0.5}
POST_DIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -m|--manifest) MANIFEST=$2; shift 2 ;;
        -p|--params)   PARAMS=$2;   shift 2 ;;
        -s|--sbatch)   SBATCH=$2;   shift 2 ;;
        -P|--project)  PROJECT=$2;  shift 2 ;;
        -i|--sif)      SIF=$2;      shift 2 ;;
        --post)        POST_DIR=$2; shift 2 ;;
        -h|--help)     sed -n '2,26p' "${BASH_SOURCE[0]}"; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

ERRORS=0
WARNINGS=0
err()  { echo "  [FAIL] $*"; ERRORS=$((ERRORS+1)); }
warn() { echo "  [WARN] $*"; WARNINGS=$((WARNINGS+1)); }
ok()   { echo "  [ ok ] $*"; }
sec()  { echo; echo "== $* =="; }

# ---------------------------------------------------------------------------
# --post: checksum the produced catalogues and report duplicate realizations.
# ---------------------------------------------------------------------------
if [[ -n "$POST_DIR" ]]; then
    echo "== post-run check: $POST_DIR =="
    [[ -d "$POST_DIR" ]] || { echo "  [FAIL] not a directory: $POST_DIR"; exit 1; }
    mapfile -t DRAWN < <(find "$POST_DIR" -name '*.drawn' -type f | sort)
    if [[ ${#DRAWN[@]} -eq 0 ]]; then
        echo "  [FAIL] no .drawn files found under $POST_DIR"; exit 1
    fi
    TMP=$(mktemp); trap 'rm -f "$TMP"' EXIT
    for f in "${DRAWN[@]}"; do printf '%s  %s\n' "$(md5sum < "$f" | cut -d" " -f1)" "$f"; done > "$TMP"
    NTOT=${#DRAWN[@]}
    NUNI=$(cut -d' ' -f1 "$TMP" | sort -u | wc -l)
    echo "  catalogues found        : $NTOT"
    echo "  distinct realizations   : $NUNI"
    if [[ "$NUNI" -lt "$NTOT" ]]; then
        echo "  [WARN] $((NTOT-NUNI)) duplicate catalogue(s) -- these add NO statistics:"
        cut -d' ' -f1 "$TMP" | sort | uniq -d | while read -r h; do
            echo "    hash ${h:0:12}:"
            grep "^$h" "$TMP" | awk '{print "      " $2}'
        done
        exit 1
    fi
    echo "  [ ok ] every catalogue is a distinct realization"
    exit 0
fi

echo "CLUMPY -g8 ensemble pre-flight"
echo "  project : $PROJECT"

# ---------------------------------------------------------------------------
sec "files and container"
RUNNER=$(command -v apptainer || command -v singularity || true)
[[ -n "$RUNNER" ]] && ok "container runtime: $RUNNER" \
                   || warn "apptainer/singularity not on PATH here (fine on a login node without the module loaded)"
[[ -r "$SIF" ]] && ok "image: $SIF" || warn "image not readable here: $SIF"

MF="$PROJECT/$MANIFEST"; PF="$PROJECT/$PARAMS"; SB="$PROJECT/$SBATCH"
[[ -r "$MF" ]] && ok "manifest: $MANIFEST" || { err "manifest missing: $MF"; echo; echo "cannot continue"; exit 1; }
[[ -r "$PF" ]] && ok "params:   $PARAMS"   || err "param file missing: $PF  (EVERY array task will exit 1 immediately)"
[[ -r "$SB" ]] && ok "sbatch:   $SBATCH"   || warn "driver script not found: $SB"

# ---------------------------------------------------------------------------
sec "manifest structure"
DATA=$(grep -vE '^[[:space:]]*(#|$)' "$MF")
N=$(printf '%s\n' "$DATA" | grep -c . )
ok "data rows: $N"

BADCOLS=$(printf '%s\n' "$DATA" | awk 'NF!=11 {print NR": "NF" columns"}')
[[ -z "$BADCOLS" ]] && ok "all rows have 11 columns" \
                    || { err "rows with wrong column count:"; printf '%s\n' "$BADCOLS" | sed 's/^/         /'; }

DUPID=$(printf '%s\n' "$DATA" | awk '{print $1}' | sort | uniq -d)
[[ -z "$DUPID" ]] && ok "all IDs unique" || { err "duplicate IDs: $(echo $DUPID | tr '\n' ' ')"; }

BADSEED=$(printf '%s\n' "$DATA" | awk '$2 !~ /^[0-9]+$/ {print $1" ("$2")"}')
[[ -z "$BADSEED" ]] && ok "all seeds are integers" || err "non-integer seeds: $BADSEED"
ZEROSEED=$(printf '%s\n' "$DATA" | awk '$2==0 {print $1}')
[[ -z "$ZEROSEED" ]] && ok "no clock-seeded (SEED=0) rows" \
                     || warn "SEED=0 is clock-seeded and NOT reproducible: $(echo $ZEROSEED | tr '\n' ' ')"

# ---------------------------------------------------------------------------
sec "realization independence  (the check that matters)"
# The clump draw depends on SEED + NCLUMPS + MMIN + MMAXFRAC only. Rows sharing
# that key produce identical .drawn catalogues regardless of DM mass, line of
# sight or injected galaxy list.
KEYS=$(printf '%s\n' "$DATA" | awk '{print $2"|"$4"|"$5"|"$6}')
NKEY=$(printf '%s\n' "$KEYS" | sort -u | wc -l)
echo "  rows: $N   distinct clump realizations: $NKEY"
if [[ "$NKEY" -lt "$N" ]]; then
    warn "$((N-NKEY)) row(s) repeat a realization already drawn by another row."
    echo "         These cost full runtime and add ZERO clump statistics."
    printf '%s\n' "$KEYS" | sort | uniq -d | while IFS='|' read -r s nc mn fr; do
        ids=$(printf '%s\n' "$DATA" | awk -v s="$s" -v nc="$nc" -v mn="$mn" -v fr="$fr" \
              '$2==s && $4==nc && $5==mn && $6==fr {printf "%s ", $1}')
        echo "         seed=$s nclumps=$nc mmin=$mn mmaxfrac=$fr  ->  $ids"
    done
else
    ok "every row is an independent realization"
fi

# ---------------------------------------------------------------------------
sec "statistical homogeneity"
NCFG=$(printf '%s\n' "$DATA" | awk '{print $4"|"$5"|"$6}' | sort -u | wc -l)
if [[ "$NCFG" -eq 1 ]]; then
    ok "single clump-population configuration -- all $NKEY realizations are poolable"
else
    warn "$NCFG different clump-population configurations present."
    echo "         Rows with different NCLUMPS/MMIN/MMAXFRAC are different statistical"
    echo "         populations and must NOT be pooled into one ensemble:"
    printf '%s\n' "$DATA" | awk '{print $4"|"$5"|"$6}' | sort | uniq -c | sort -rn |
        while read -r c k; do
            IFS='|' read -r nc mn fr <<< "$k"
            ns=$(printf '%s\n' "$DATA" | awk -v nc="$nc" -v mn="$mn" -v fr="$fr" \
                 '$4==nc && $5==mn && $6==fr {print $2}' | sort -u | wc -l)
            echo "         NCLUMPS=$nc MMIN=$mn MMAXFRAC=$fr : $c rows, $ns distinct seeds"
        done
    LARGEST=$(printf '%s\n' "$DATA" | awk '{print $4"|"$5"|"$6"|"$2}' | sort -u |
              awk -F'|' '{print $1"|"$2"|"$3}' | sort | uniq -c | sort -rn | head -1 | awk '{print $1}')
    echo "         Largest poolable ensemble: $LARGEST realizations."
fi

# ---------------------------------------------------------------------------
sec "external galaxy lists (GALLIST)"
NLIST=0
while read -r _ _ _ _ _ _ _ _ _ gl _; do
    case "$gl" in
        none|NONE|-1|default|DEFAULT|"") ;;
        *) NLIST=$((NLIST+1))
           [[ -r "$PROJECT/$gl" ]] || err "list-halo file missing: $gl  (driver only WARNS; the run would use different physics)" ;;
    esac
done <<< "$DATA"
NPOLLUTE=$(printf '%s\n' "$DATA" | awk '$10!="none" && $10!="NONE" && $10!="-1"' | grep -c . )
[[ "$NLIST" -eq 0 ]] && ok "no per-row list-halo files referenced"
if [[ "$NPOLLUTE" -gt 0 ]]; then
    warn "$NPOLLUTE row(s) inject external objects into the .drawn catalogue."
    echo "         Those are not Galactic subhaloes; filter on the ntuple 'type' branch"
    echo "         before any pulsar-clump association."
else
    ok "all rows are pure Galactic draws (GALLIST=none)"
fi

# ---------------------------------------------------------------------------
sec "array sizing"
if [[ -r "$SB" ]]; then
    ARR=$(grep -oE '^#SBATCH[[:space:]]+--array=[0-9]+-[0-9]+' "$SB" | head -1 | sed 's/.*--array=//')
    if [[ -n "$ARR" ]]; then
        LO=${ARR%-*}; HI=${ARR#*-}; NSLOT=$((HI-LO+1))
        if [[ "$NSLOT" -eq "$N" ]]; then ok "--array=$ARR covers all $N rows"
        elif [[ "$NSLOT" -lt "$N" ]]; then err "--array=$ARR runs $NSLOT of $N rows -- $((N-NSLOT)) never run. Use --array=0-$((N-1))"
        else warn "--array=$ARR has $NSLOT slots for $N rows; the extra tasks exit 0 harmlessly"
        fi
    else warn "could not parse an --array range from $SBATCH"
    fi
    echo "  submit with: sbatch --array=0-$((N-1))%6 $SBATCH"
fi

# ---------------------------------------------------------------------------
sec "partition limits and disk"
REQ=$(grep -oE '^#SBATCH[[:space:]]+--time=[0-9:-]+' "$SB" 2>/dev/null | head -1 | sed 's/.*--time=//')
[[ -n "$REQ" ]] && echo "  requested --time: $REQ"
if command -v sinfo >/dev/null 2>&1; then
    echo "  partition time limits:"; sinfo -h -o "    %P %l" 2>/dev/null | sort -u | head -10
    echo "         compare by hand: a 240:00:00 (10 day) request is refused at submit"
    echo "         by partitions with a shorter MaxTime."
else
    warn "sinfo not available here -- check the partition MaxTime on the cluster:"
    echo "         sinfo -o \"%P %l\"    (a 240:00:00 request needs MaxTime >= 10 days)"
fi
NEED=$(awk -v n="$N" -v g="$GB_PER_REALIZATION" 'BEGIN{printf "%.0f", n*g}')
echo "  estimated output: ~${NEED} GB for $N realizations (~${GB_PER_REALIZATION} GB each)"
OUT="$PROJECT/output"; mkdir -p "$OUT" 2>/dev/null
AVAIL=$(df -Pk "$OUT" 2>/dev/null | awk 'NR==2{printf "%.0f", $4/1048576}')
if [[ -n "$AVAIL" ]]; then
    echo "  free space on $(df -P "$OUT" 2>/dev/null | awk 'NR==2{print $6}'): ${AVAIL} GB"
    [[ "$AVAIL" -lt "$NEED" ]] && err "not enough free space: need ~${NEED} GB, have ${AVAIL} GB"
fi

# ---------------------------------------------------------------------------
sec "driver robustness"
if [[ -r "$SB" ]] && grep -q 'cp -av.*|| true' "$SB"; then
    warn "stage-out uses 'cp ... || true', which defeats 'set -euo pipefail'."
    echo "         A failed copy after a multi-day run would report success and lose"
    echo "         the catalogue. Verify the .drawn landed, e.g.:"
    echo "           cp -av \"\$WORKDIR\"/*.drawn \"\$DEST\"/ || { echo 'stage-out failed' >&2; exit 1; }"
    echo "           ls \"\$DEST\"/*.drawn >/dev/null || { echo 'no catalogue produced' >&2; exit 1; }"
fi

echo
echo "==================================================================="
if [[ "$ERRORS" -gt 0 ]]; then
    echo "NOT READY: $ERRORS blocking error(s), $WARNINGS warning(s)."
    exit 1
fi
echo "READY: 0 blocking errors, $WARNINGS warning(s)."
echo "After the run, verify the catalogues really are distinct:"
echo "  $0 --post \$PROJECT/output"
exit 0
