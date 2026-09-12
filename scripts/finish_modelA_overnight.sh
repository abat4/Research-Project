#!/usr/bin/env bash

set -u

ROOT="$HOME/Research-Project"
LOOSE="$ROOT/calculations/step2_3_bare_clusters/model_A/03_optimization/BS1_loose_preopt"
FINAL="$ROOT/calculations/step2_3_bare_clusters/model_A/03_optimization/BS1_final_refinement"
ANALYSIS="$ROOT/analysis"
LOOSE_OUT="$LOOSE/A_BS1_loose.out"

stamp() {
    date -Is
}

echo "========================================"
echo "MODEL A AUTOMATIC OVERNIGHT PIPELINE"
echo "Started: $(stamp)"
echo "========================================"

echo
echo "[1] Waiting for current LOOSE optimization..."

while pgrep -f '[n]wchem A_BS1_loose.nw' >/dev/null; do
    echo "$(stamp) - loose optimization still running"
    sleep 60
done

echo "$(stamp) - loose NWChem process ended"

cd "$LOOSE" || exit 1
date -Is > optimization_end_time.txt

echo
echo "[2] Checking whether LOOSE optimization converged..."

grep '^@' A_BS1_loose.out | tail -20 > loose_final_opt_table.txt

if ! grep -qiE 'optimization converged|geometry optimization converged' A_BS1_loose.out; then
    echo "LOOSE_NOT_CONVERGED" > overnight_status.txt
    echo "$(stamp) - LOOSE optimization did NOT report convergence."
    echo "Automatic pipeline stopping safely."
    tail -200 A_BS1_loose.out > loose_failure_tail.txt
    exit 2
fi

echo "LOOSE_CONVERGED" > overnight_status.txt
echo "$(stamp) - LOOSE optimization converged."

LASTXYZ=$(ls A_BS1_loose-*.xyz 2>/dev/null | sort -V | tail -1)

if [ -z "${LASTXYZ:-}" ] || [ ! -s "$LASTXYZ" ]; then
    echo "ERROR_NO_LOOSE_XYZ" >> overnight_status.txt
    exit 3
fi

if [ ! -s A_BS1_loose.movecs ]; then
    echo "ERROR_NO_LOOSE_MOVECS" >> overnight_status.txt
    exit 4
fi

echo "Latest loose geometry: $LASTXYZ"

mkdir -p "$FINAL"

cp "$LASTXYZ" "$FINAL/start.xyz"
cp A_BS1_loose.movecs "$FINAL/start.movecs"

sha256sum \
    "$FINAL/start.xyz" \
    "$FINAL/start.movecs" \
    > "$FINAL/starting_files.sha256"

echo
echo "[3] Creating DEFAULT final-refinement input..."

cat > "$FINAL/A_BS1_final.nw" <<'EONW'
echo

start A_BS1_final

title "Model A BS1 final default geometry refinement"

memory total 1 gb

charge -2

geometry units angstrom noautoz nocenter noautosym
    load format xyz start.xyz
end

basis spherical
    * library def2-svp
end

dft
    odft
    mult 1
    xc becke88 perdew86

    iterations 150

    convergence fast energy 1e-5 density 1e-4 gradient 1e-3

    vectors input start.movecs output A_BS1_final.movecs

    mulliken
end

driver
    default
    maxiter 30
    xyz A_BS1_final
end

task dft optimize
EONW

cd "$FINAL" || exit 1

cat > job_notes.md <<'EONOTES'
Model A final DEFAULT geometry refinement.

Starting structure:
Converged geometry from the preceding LOOSE BP86/def2-SVP
preoptimization.

Electronic state:
Representative BS1 broken-symmetry determinant.

Charge: -2
Multiplicity: 1
Method: BP86/def2-SVP

The purpose of this run is to refine the preoptimized geometry using
NWChem DEFAULT geometry convergence criteria.

The final structure is accepted only if the DEFAULT optimization itself
reports convergence and the broken-symmetry Fe spin pattern is retained.
EONOTES

echo
echo "[4] Starting DEFAULT refinement: $(stamp)"

export OMP_NUM_THREADS=1
date -Is > optimization_start_time.txt

/usr/bin/mpirun.openmpi -np 4 nwchem A_BS1_final.nw \
    > A_BS1_final.out 2>&1

RUN_STATUS=$?

date -Is > optimization_end_time.txt

echo "NWChem exit status: $RUN_STATUS"

echo
echo "[5] Checking final refinement..."

if ! grep -qiE 'optimization converged|geometry optimization converged' A_BS1_final.out; then

    echo "FINAL_REFINEMENT_NOT_CONVERGED" > overnight_status.txt

    LASTFINAL=$(ls A_BS1_final-*.xyz 2>/dev/null | sort -V | tail -1)

    if [ -n "${LASTFINAL:-}" ]; then
        cp "$LASTFINAL" latest_unconverged_geometry.xyz
    fi

    {
        echo "MODEL A FINAL REFINEMENT — NOT CONVERGED"
        echo
        echo "NWChem exit status: $RUN_STATUS"
        echo
        echo "OPTIMIZATION TABLE"
        grep '^@' A_BS1_final.out | tail -30
        echo
        echo "LAST ENERGY"
        grep 'Total DFT energy' A_BS1_final.out | tail -1
        echo
        echo "LAST S2"
        grep '<S2>' A_BS1_final.out | tail -1
    } > model_A_final_summary.txt

    cp model_A_final_summary.txt \
       "$ANALYSIS/model_A_final_summary.txt"

    echo "$(stamp) - DEFAULT refinement did not converge."
    echo "Pipeline stopped safely."
    exit 5
fi

echo "FINAL_REFINEMENT_CONVERGED" > overnight_status.txt

LASTFINAL=$(ls A_BS1_final-*.xyz 2>/dev/null | sort -V | tail -1)

if [ -z "${LASTFINAL:-}" ] || [ ! -s "$LASTFINAL" ]; then
    echo "ERROR_NO_FINAL_XYZ" >> overnight_status.txt
    exit 6
fi

cp "$LASTFINAL" optimized_model_A.xyz

LINE=$(grep -n "Spin Density - Mulliken Population Analysis" \
       A_BS1_final.out | tail -1 | cut -d: -f1)

if [ -n "${LINE:-}" ]; then
    sed -n "${LINE},$((LINE+40))p" \
        A_BS1_final.out > optimized_spin_density.txt
fi

{
    echo "MODEL A FINAL RESULT"
    echo
    echo "Status: DEFAULT GEOMETRY OPTIMIZATION CONVERGED"
    echo
    echo "Final geometry file: optimized_model_A.xyz"
    echo
    echo "FINAL OPTIMIZATION TABLE"
    grep '^@' A_BS1_final.out | tail -30
    echo
    echo "FINAL DFT ENERGY"
    grep 'Total DFT energy' A_BS1_final.out | tail -1
    echo
    echo "FINAL S2"
    grep '<S2>' A_BS1_final.out | tail -1
    echo
    echo "FINAL SPIN DENSITY"
    cat optimized_spin_density.txt 2>/dev/null || true
} > model_A_final_summary.txt

sha256sum \
    optimized_model_A.xyz \
    A_BS1_final.nw \
    A_BS1_final.out \
    A_BS1_final.movecs \
    model_A_final_summary.txt \
    > final_checksums.sha256

cp model_A_final_summary.txt \
   "$ANALYSIS/model_A_final_summary.txt"

cp optimized_model_A.xyz \
   "$ANALYSIS/model_A_optimized.xyz"

cat >> "$ROOT/notes/research_log.md" <<'EOLOG'

## Model A bare-cluster geometry optimization completed

The representative BS1 state was first preoptimized using the LOOSE
NWChem geometry criteria. The converged preoptimized structure was then
refined using DEFAULT geometry convergence criteria at BP86/def2-SVP.

The final geometry, electronic energy, <S^2>, local spin-density
analysis, and checksums are preserved in the Model A final-refinement
directory and the analysis directory.

EOLOG

echo
echo "========================================"
echo "MODEL A AUTOMATIC PIPELINE COMPLETE"
echo "Finished: $(stamp)"
echo "========================================"
echo
cat model_A_final_summary.txt
