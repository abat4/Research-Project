# Model A BS1 geometry optimization

System:
[Fe4S4(SCH3)4]2-

Purpose:
Optimize the bare Model A structure before CO2 is introduced.

Starting electronic state:
Representative unconstrained broken-symmetry BS1 state (++--).

Reason for choosing BS1:
BS1, BS2, and BS3 differed by only 0.004374 kJ/mol in the def2-SVP
screening calculations and were therefore treated as energetically
degenerate. BS1 was selected as a representative state.

Method:
BP86/def2-SVP

Charge:
-2

Multiplicity:
1 broken-symmetry determinant

SCF criteria:
energy = 1e-5 Hartree
density = 1e-4
gradient = 1e-3

Geometry optimizer:
NWChem DRIVER
default convergence criteria
maximum 40 optimization steps

Parallelization:
4 MPI processes
OMP_NUM_THREADS=1

Important outputs:
A_BS1_opt.out
A_BS1_opt.movecs
A_BS1_opt-###.xyz optimization trajectory
optimized.xyz after completion

This is the def2-SVP geometry optimization. Selected final structures
can later receive higher-level def2-TZVP single-point calculations.
