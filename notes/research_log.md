# Computational Research Log

## Software setup

NWChem was installed under Ubuntu 24.04 using WSL2.
Open Babel was installed under Ubuntu.
NWChem functionality was verified using a B3LYP/6-31G single-point calculation on water.

Water test energy:
-76.3780475832444 Hartree

## Step 1

Three starting molecular models were prepared.

Model A:
[Fe4S4(SCH3)4]2- benchmark cluster.

Model B:
Truncated CODH Ni-Fe-S C-cluster model derived from PDB 3B52.

Model C:
Model B with Ni replaced by Fe while retaining the same initial coordinates.

Original structures are retained unchanged in structures/.

## Steps 2-3

Goal:
Determine reasonable electronic/spin states of the bare metal clusters and optimize the viable states before adding CO2.

Primary DFT method:
BP86 / def2-TZVP in NWChem.

Broken-symmetry calculations will be used for the Fe-S clusters because local high-spin Fe centers may couple antiferromagnetically.

# Computational Research Log

## Objective
Determine how replacing Ni with Fe changes CO2 binding and activation
in Fe-S clusters relevant to primitive carbon fixation.

## Step 1
Three starting structures were generated:
A: [Fe4S4(SCH3)4]2-
B: 3B52-derived CODH NiFe4S4 first-shell model
C: Model B with Ni replaced by Fe

Original structures are preserved unchanged in structures/.

## Computational method
Primary program: NWChem
Primary geometry method: BP86/def2-TZVP
Spin treatment: unrestricted/broken-symmetry DFT
CO2 is excluded during bare-cluster calculations.

Functional: BP86
NWChem: xc becke88 perdew86
Basis: def2-TZVP



## Model A high-spin screening reference — completed

A maximum-spin single-point DFT calculation was completed for
[Fe4S4(SCH3)4]2- using BP86/def2-SVP.

Charge: -2
Multiplicity: 19
Total spin: S = 9
SCF iterations: 25
Final DFT energy: -8399.820443085609 Hartree
<S^2>: 90.0364
Exact S(S+1): 90.0000

The calculation satisfied the screening convergence criteria.
This state is retained as the maximum-spin reference wavefunction
for construction of broken-symmetry states and is not interpreted
as the physical ground state of the cluster.

## Model A BS1 broken-symmetry calculation

The first Model A broken-symmetry pattern was initialized with Fe1 and
Fe2 spin-up and Fe3 and Fe4 spin-down. Local-spin constraints were used
only to generate the starting electronic state. The constraints were
then removed and an unrestricted BP86/def2-SVP single-point calculation
was performed.

Final unconstrained energy:
TO BE FILLED

<S^2>:
TO BE FILLED

Final Fe spin pattern:
TO BE FILLED

Status:
TO BE FILLED

## Model A BS1 broken-symmetry calculation

The first broken-symmetry arrangement was initialized with Fe1 and Fe2
spin-up and Fe3 and Fe4 spin-down. Local-spin constraints were used
only to generate the starting determinant. The constraints were then
removed and an unrestricted BP86/def2-SVP single-point calculation was
performed.

Final energy:
TO BE FILLED

<S^2>:
TO BE FILLED

Final Fe spin pattern:
TO BE FILLED

Status:
TO BE FILLED

## Model A broken-symmetry screening completed

Three 2-up/2-down broken-symmetry arrangements were generated from
the converged maximum-spin reference. Local-spin constraints were used
only to construct the initial determinants, after which the constraints
were removed and each state was converged independently using
unrestricted BP86/def2-SVP DFT.

BS1 target pattern: ++--
BS2 target pattern: +-+-
BS3 target pattern: +--+

The unconstrained energies, <S^2> values, and final Fe spin densities
were saved in analysis/model_A_BS_comparison.csv and the corresponding
BS result files. The states were compared before choosing structures
for geometry optimization.

## Model A broken-symmetry screening comparison

Three unconstrained 2-up/2-down broken-symmetry states were compared
using BP86/def2-SVP.

BS1: -8399.879385632847 Hartree
BS2: -8399.879385718777 Hartree
BS3: -8399.879387298866 Hartree

The total energy spread was only 0.004374 kJ/mol, smaller than the
numerical resolution implied by the screening SCF threshold. The three
broken-symmetry solutions were therefore treated as energetically
degenerate rather than assigning significance to their numerical
ordering.

BS1 was selected as a representative state for geometry optimization
because it had already been verified to retain large local Fe spin
densities with the ++-- pattern.
