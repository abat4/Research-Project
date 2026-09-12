#!/usr/bin/env bash

set -u

BS="$1"

ROOT="$HOME/Research-Project/calculations/step2_3_bare_clusters/model_A/02_broken_symmetry"
STRUCT="$HOME/Research-Project/structures/model_A_Fe4S4/model_A_Fe4S4.xyz"
REF="$HOME/Research-Project/calculations/step2_3_bare_clusters/model_A/01_high_spin_SVP/attempt3_clean/A_HS_SVP3_MASTER.movecs"

case "$BS" in
    BS2)
        S1="3.5"
        S2="-3.5"
        S3="3.5"
        S4="-3.5"
        PATTERN="+-+-"
        ;;
    BS3)
        S1="3.5"
        S2="-3.5"
        S3="-3.5"
        S4="3.5"
        PATTERN="+--+"
        ;;
    *)
        echo "Usage: run_modelA_bs.sh BS2 or BS3"
        exit 2
        ;;
esac

DIR="$ROOT/$BS"
mkdir -p "$DIR"
cd "$DIR" || exit 1

cp "$STRUCT" start.xyz
cp "$REF" A_HS_reference.movecs
sha256sum start.xyz > start_structure.sha256

cat > "${BS}_seed.nw" <<EONW
echo

start A_${BS}_seed

title "Model A ${BS} constrained spin seed"

memory total 800 mb

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

    convergence nolevelshifting
    convergence energy 1e-5 density 1e-4 gradient 1e-3

    vectors input A_HS_reference.movecs output ${BS}_seed.movecs

    cdft 1 1 spin ${S1} pop mulliken
    cdft 2 2 spin ${S2}
    cdft 3 3 spin ${S3}
    cdft 4 4 spin ${S4}

    mulliken
end

set dft:cdft_maxiter 200

task dft energy
EONW

cat > "${BS}_relax.nw" <<EONW
echo

start A_${BS}_relax

title "Model A ${BS} unconstrained broken-symmetry screening"

memory total 800 mb

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

    vectors input ${BS}_seed.movecs output ${BS}.movecs

    mulliken
end

task dft energy
EONW

cat > job_notes.md <<EONOTES
# Model A ${BS}

Target Fe spin pattern: ${PATTERN}

Fe1 seed target: ${S1}
Fe2 seed target: ${S2}
Fe3 seed target: ${S3}
Fe4 seed target: ${S4}

Method: BP86/def2-SVP
Charge: -2
Multiplicity: 1
MPI processes: 2
Memory: 800 MB per MPI process

Workflow:
1. Constrained local-spin seed
2. Unconstrained broken-symmetry single point

Only the unconstrained calculation is used for state-energy comparison.
EONOTES

export OMP_NUM_THREADS=1

{
    echo "STATE=$BS"
    echo "TARGET=$PATTERN"
    echo "PIPELINE_START=$(date -Is)"
    echo "STAGE=SEED"
} > pipeline_status.txt

date -Is > seed_start_time.txt

/usr/bin/mpirun.openmpi -np 2 nwchem "${BS}_seed.nw" > "${BS}_seed.out" 2>&1
SEED_STATUS=$?

date -Is > seed_end_time.txt
echo "SEED_EXIT_STATUS=$SEED_STATUS" >> pipeline_status.txt

if ! grep -q "Total DFT energy" "${BS}_seed.out"; then
    echo "FAILED_SEED_NO_FINAL_ENERGY" >> pipeline_status.txt
    exit 1
fi

if [ ! -s "${BS}_seed.movecs" ]; then
    echo "FAILED_SEED_NO_MOVECS" >> pipeline_status.txt
    exit 1
fi

LINE=$(grep -n "Spin Density - Mulliken Population Analysis" "${BS}_seed.out" | tail -1 | cut -d: -f1)

if [ -n "${LINE:-}" ]; then
    sed -n "${LINE},$((LINE+40))p" "${BS}_seed.out" > "${BS}_seed_spin_density.txt"
fi

echo "STAGE=RELAX" >> pipeline_status.txt

date -Is > relax_start_time.txt

/usr/bin/mpirun.openmpi -np 2 nwchem "${BS}_relax.nw" > "${BS}_relax.out" 2>&1
RELAX_STATUS=$?

date -Is > relax_end_time.txt
echo "RELAX_EXIT_STATUS=$RELAX_STATUS" >> pipeline_status.txt

if ! grep -q "Total DFT energy" "${BS}_relax.out"; then
    echo "FAILED_RELAX_NO_FINAL_ENERGY" >> pipeline_status.txt
    exit 1
fi

if [ ! -s "${BS}.movecs" ]; then
    echo "FAILED_RELAX_NO_MOVECS" >> pipeline_status.txt
    exit 1
fi

LINE=$(grep -n "Spin Density - Mulliken Population Analysis" "${BS}_relax.out" | tail -1 | cut -d: -f1)

if [ -n "${LINE:-}" ]; then
    sed -n "${LINE},$((LINE+40))p" "${BS}_relax.out" > "${BS}_spin_density.txt"
fi

{
    echo "MODEL A — ${BS}"
    echo
    echo "TARGET PATTERN: ${PATTERN}"
    echo
    echo "FINAL ENERGY"
    grep "Total DFT energy" "${BS}_relax.out" | tail -1
    echo
    echo "S2"
    grep "<S2>" "${BS}_relax.out" | tail -1
    echo
    echo "FINAL SCF"
    grep -E "d=.*diis" "${BS}_relax.out" | tail -5
    echo
    echo "FINAL SPIN DENSITY"
    cat "${BS}_spin_density.txt"
} > "${BS}_result_summary.txt"

sha256sum \
    start.xyz \
    "${BS}_seed.nw" \
    "${BS}_seed.out" \
    "${BS}_seed.movecs" \
    "${BS}_relax.nw" \
    "${BS}_relax.out" \
    "${BS}.movecs" \
    > "${BS}_checksums.sha256"

echo "SUCCESS" >> pipeline_status.txt
echo "PIPELINE_END=$(date -Is)" >> pipeline_status.txt
