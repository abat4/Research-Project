# A_HS_SVP

Model:
Model A — [Fe4S4(SCH3)4]2-

Purpose:
Generate a maximum-spin reference wavefunction for broken-symmetry screening.

Charge:
-2

Multiplicity:
19

Functional:
BP86

Basis:
def2-SVP

Calculation:
Single-point DFT

Geometry optimization:
No

SCF stabilization:
40% damping
0.5 Hartree level shift

Reason for method choice:
An initial BP86/def2-TZVP calculation was substantially more expensive on
the available laptop hardware and did not reach stable convergence within
the attempted runtime. BP86/def2-SVP was therefore adopted for spin-state
screening and geometry optimization. Selected final structures may later
receive def2-TZVP single-point calculations.
EOF4
date -Is > start_time.txt
nwchem A_HS_SVP.nw > A_HS_SVP.out 2>&1
