#!/usr/bin/env bash

set -u

MODEL="$1"

ROOT="$HOME/Research-Project"
CALC="$ROOT/calculations/step2_3_bare_clusters"
ANALYSIS="$ROOT/analysis"

case "$MODEL" in
    B)
        SRC="$ROOT/structures/model_B_CODH_NiFeS/model_B_CODH_NiFeS.xyz"
        NAME="model_B"
        ;;
    C)
        SRC="$ROOT/structures/model_C_CODH_Fe_control/model_C_CODH_Fe_control.xyz"
        NAME="model_C"
        ;;
    *)
        echo "Usage: run_BC_deadline.sh B"
        echo "   or: run_BC_deadline.sh C"
        exit 2
        ;;
esac

DIR="$CALC/$NAME/deadline_pipeline"

mkdir -p "$DIR"
mkdir -p "$ANALYSIS"

cd "$DIR" || exit 1

cp "$SRC" start.xyz

CHARGE=-2
MPI=2

echo "MODEL=$MODEL" > pipeline_status.txt
echo "CHARGE=$CHARGE" >> pipeline_status.txt
echo "START=$(date -Is)" >> pipeline_status.txt
echo "STAGE=SPIN_SCREEN" >> pipeline_status.txt

cat > modeling_assumptions.txt <<EOF2
Model ${MODEL} deadline-mode electronic-structure assumptions

Overall charge: -2
Charge status: computational modeling assumption, not assigned by Step 1.

Tested multiplicities:
2
4
6

Purpose:
Screen several plausible spin-polarized states rapidly and select the
lowest-energy state that successfully completes under the same protocol.

Geometry convergence protocol:
Gmax < 0.005 a.u.
Grms < 0.003 a.u.
Displacement thresholds nonrestrictive.

Geometry optimization wall-time limit:
25 minutes.

The same charge and screening protocol are applied to Models B and C
to preserve a matched Ni-versus-Fe comparison.
EOF2

sha256sum start.xyz > starting_structure.sha256

echo "multiplicity,energy_hartree,S2,status" > spin_screen.csv

export OMP_NUM_THREADS=1

for MULT in 2 4 6
do

cat > "${MODEL}_m${MULT}.nw" <<EOF2
echo

start ${MODEL}_m${MULT}

title "Model ${MODEL} multiplicity ${MULT} screening"

memory total 700 mb

charge ${CHARGE}

geometry units angstrom noautoz nocenter noautosym
    load format xyz start.xyz
end

basis spherical
    * library def2-svp
end

dft
    odft
    mult ${MULT}
    xc becke88 perdew86

    iterations 100

    convergence fast energy 1e-4 density 5e-4 gradient 5e-3

    vectors output ${MODEL}_m${MULT}.movecs

    mulliken
end

task dft energy
EOF2

echo "$(date -Is) starting multiplicity ${MULT}" >> pipeline.log

timeout --signal=TERM --kill-after=30s 10m \
/usr/bin/mpirun.openmpi -np ${MPI} \
nwchem "${MODEL}_m${MULT}.nw" \
> "${MODEL}_m${MULT}.out" 2>&1

EXIT=$?

ENERGY=$(grep "Total DFT energy" "${MODEL}_m${MULT}.out" \
         | tail -1 | awk '{print $5}')

S2=$(grep "<S2>" "${MODEL}_m${MULT}.out" \
     | tail -1 | awk '{print $3}')

if [ -n "${ENERGY:-}" ]; then
    echo "${MULT},${ENERGY},${S2:-NA},SUCCESS" >> spin_screen.csv
else
    if [ "$EXIT" -eq 124 ]; then
        STATUS="TIMEOUT"
    else
        STATUS="FAILED"
    fi

    echo "${MULT},,NA,${STATUS}" >> spin_screen.csv
fi

done

BEST=$(python3 - <<'PY'
import csv

successful = []

with open("spin_screen.csv") as f:
    for row in csv.DictReader(f):
        if row["status"] != "SUCCESS":
            continue
        try:
            successful.append(
                (int(row["multiplicity"]),
                 float(row["energy_hartree"]))
            )
        except ValueError:
            pass

if successful:
    print(min(successful, key=lambda x: x[1])[0])
PY
)

if [ -z "${BEST:-}" ]; then
    echo "FAILED_NO_SPIN_STATE" >> pipeline_status.txt
    echo "No spin-screen calculation produced a final energy."
    exit 10
fi

echo "WINNING_MULTIPLICITY=$BEST" >> pipeline_status.txt
echo "STAGE=FAST_OPTIMIZATION" >> pipeline_status.txt

cp "${MODEL}_m${BEST}.movecs" winner.movecs

cat > "${MODEL}_fast_opt.nw" <<EOF2
echo

start ${MODEL}_fast_opt

title "Model ${MODEL} deadline geometry optimization"

memory total 700 mb

charge ${CHARGE}

geometry units angstrom noautoz nocenter noautosym
    load format xyz start.xyz
end

basis spherical
    * library def2-svp
end

dft
    odft
    mult ${BEST}
    xc becke88 perdew86

    iterations 100

    convergence fast energy 1e-4 density 5e-4 gradient 5e-3

    vectors input winner.movecs output ${MODEL}_fast_opt.movecs

    mulliken
end

driver
    gmax 0.005
    grms 0.003

    xrms 1.0
    xmax 1.0

    trust 0.05
    maxiter 3

    xyz ${MODEL}_fast_opt
end

task dft optimize
EOF2

date -Is > optimization_start.txt

timeout --signal=TERM --kill-after=30s 25m \
/usr/bin/mpirun.openmpi -np ${MPI} \
nwchem "${MODEL}_fast_opt.nw" \
> "${MODEL}_fast_opt.out" 2>&1

OPT_EXIT=$?

date -Is > optimization_end.txt

LASTXYZ=$(ls "${MODEL}_fast_opt"-*.xyz 2>/dev/null \
          | sort -V | tail -1)

if [ -n "${LASTXYZ:-}" ] && [ -s "$LASTXYZ" ]; then
    cp "$LASTXYZ" final_geometry.xyz
    GEOM_SOURCE="$LASTXYZ"
else
    cp start.xyz final_geometry.xyz
    GEOM_SOURCE="start.xyz"
fi

LASTLINE=$(grep '^@' "${MODEL}_fast_opt.out" \
           | awk '$2 ~ /^[0-9]+$/ {line=$0} END {print line}')

if [ -n "${LASTLINE:-}" ]; then
    GMAX=$(echo "$LASTLINE" | awk '{print $5}')
    GRMS=$(echo "$LASTLINE" | awk '{print $6}')

    OPT_STATUS=$(python3 - "$GMAX" "$GRMS" <<'PY'
import sys

gmax=float(sys.argv[1])
grms=float(sys.argv[2])

if gmax <= 0.005 and grms <= 0.003:
    print("CONVERGED_UNDER_PROJECT_CRITERIA")
else:
    print("TIME_LIMITED_BEST_GEOMETRY")
PY
)

else
    GMAX="NA"
    GRMS="NA"
    OPT_STATUS="NO_ACCEPTED_OPTIMIZATION_STEP"
fi

echo "OPT_STATUS=$OPT_STATUS" >> pipeline_status.txt
echo "OPT_GMAX=$GMAX" >> pipeline_status.txt
echo "OPT_GRMS=$GRMS" >> pipeline_status.txt
echo "GEOMETRY_SOURCE=$GEOM_SOURCE" >> pipeline_status.txt

if [ -s "${MODEL}_fast_opt.movecs" ]; then
    cp "${MODEL}_fast_opt.movecs" final_start.movecs
else
    cp winner.movecs final_start.movecs
fi

echo "STAGE=FINAL_SINGLE_POINT" >> pipeline_status.txt

cat > "${MODEL}_final_sp.nw" <<EOF2
echo

start ${MODEL}_final_sp

title "Model ${MODEL} final fixed-geometry single point"

memory total 700 mb

charge ${CHARGE}

geometry units angstrom noautoz nocenter noautosym
    load format xyz final_geometry.xyz
end

basis spherical
    * library def2-svp
end

dft
    odft
    mult ${BEST}
    xc becke88 perdew86

    iterations 150

    convergence fast energy 1e-5 density 1e-4 gradient 1e-3

    vectors input final_start.movecs output ${MODEL}_final.movecs

    mulliken
end

task dft energy
EOF2

date -Is > final_sp_start.txt

timeout --signal=TERM --kill-after=30s 15m \
/usr/bin/mpirun.openmpi -np ${MPI} \
nwchem "${MODEL}_final_sp.nw" \
> "${MODEL}_final_sp.out" 2>&1

SP_EXIT=$?

date -Is > final_sp_end.txt

FINAL_E=$(grep "Total DFT energy" "${MODEL}_final_sp.out" \
          | tail -1 | awk '{print $5}')

FINAL_S2=$(grep "<S2>" "${MODEL}_final_sp.out" \
           | tail -1 | awk '{print $3}')

LINE=$(grep -n "Spin Density - Mulliken Population Analysis" \
       "${MODEL}_final_sp.out" | tail -1 | cut -d: -f1)

if [ -n "${LINE:-}" ]; then
    sed -n "${LINE},$((LINE+90))p" \
    "${MODEL}_final_sp.out" > final_spin_density.txt
fi

{
echo "MODEL ${MODEL} DEADLINE PIPELINE RESULT"
echo
echo "Charge assumption: ${CHARGE}"
echo
echo "SPIN SCREEN"
cat spin_screen.csv
echo
echo "Selected multiplicity: ${BEST}"
echo
echo "GEOMETRY"
echo "Status: ${OPT_STATUS}"
echo "Gmax: ${GMAX}"
echo "Grms: ${GRMS}"
echo "Geometry file: final_geometry.xyz"
echo "Geometry source: ${GEOM_SOURCE}"
echo
echo "FINAL SINGLE-POINT ENERGY"
echo "${FINAL_E:-NO_FINAL_ENERGY}"
echo
echo "FINAL S2"
echo "${FINAL_S2:-NA}"
echo
echo "FINAL SPIN DENSITY"
cat final_spin_density.txt 2>/dev/null || true
} > "${MODEL}_deadline_summary.txt"

sha256sum \
start.xyz \
final_geometry.xyz \
"${MODEL}_deadline_summary.txt" \
> final_checksums.sha256

cp final_geometry.xyz \
"$ANALYSIS/${NAME}_optimized.xyz"

cp "${MODEL}_deadline_summary.txt" \
"$ANALYSIS/${NAME}_deadline_summary.txt"

cp spin_screen.csv \
"$ANALYSIS/${NAME}_spin_screen.csv"

echo "SUCCESS" >> pipeline_status.txt
echo "END=$(date -Is)" >> pipeline_status.txt

echo
echo "============================="
echo "MODEL ${MODEL} PIPELINE FINISHED"
echo "============================="
cat "${MODEL}_deadline_summary.txt"
