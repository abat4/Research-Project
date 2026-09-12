# Steps 2-3 Computational Method

Program:
NWChem

Starting coordinates:
Structures generated during Step 1.

Exchange-correlation functional:
BP86

NWChem specification:
xc becke88 perdew86

Basis:
def2-TZVP

Model A:
Charge = -2.
Initial calculation = maximum-spin reference.
Maximum-spin multiplicity = 19.
Initial calculation is single-point only.

Model B:
CODH-derived Ni-Fe-S model.
Electronic/redox state must be explicitly defined before production optimization.

Model C:
Fe-substituted control derived from Model B.
Model B and Model C will use identical computational settings whenever they are directly compared.

Spin treatment:
Generate high-spin reference orbitals first.
Construct broken-symmetry states from those orbitals.
Compare converged physically reasonable electronic solutions.
Optimize viable low-energy states.

CO2:
Not present during Steps 2-3.
