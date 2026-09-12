# Model A BS2

Target Fe spin pattern: +-+-

Fe1 seed target: 3.5
Fe2 seed target: -3.5
Fe3 seed target: 3.5
Fe4 seed target: -3.5

Method: BP86/def2-SVP
Charge: -2
Multiplicity: 1
MPI processes: 2
Memory: 800 MB per MPI process

Workflow:
1. Constrained local-spin seed
2. Unconstrained broken-symmetry single point

Only the unconstrained calculation is used for state-energy comparison.
