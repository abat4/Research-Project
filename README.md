# Fe/Ni-S Cluster CO2 Activation Project

Computational study of how nickel incorporation into iron-sulfur clusters
influences CO2 binding and activation in models inspired by the CODH C-cluster.

## Models

### Model A
[Fe4S4(SCH3)4]2- cubane reference model.

### Model B
First-shell Ni-Fe-S CODH C-cluster model derived from PDB 3B52.

### Model C
Matched Fe-substituted control generated from Model B by replacing the Ni site
with Fe while initially retaining identical coordinates.

## Current workflow

1. Build and validate starting structures
2. Optimize bare-cluster geometries
3. Examine electronic/spin states
4. Add CO2
5. Analyze CO2 activation
6. Analyze charge transfer
7. Analyze orbital behavior
8. Test redox/protonation effects
9. Screen CO2-to-COOH reaction coordinates
10. Compare Ni-containing and Fe-only models

## Computational methods

Early Model A calculations used BP86/def2-SVP in NWChem, including
broken-symmetry calculations.

For rapid comparative geometry screening, GFN1-xTB with ALPB water is being
used consistently across Models A, B, and C.

Large NWChem scratch files and restart binaries are excluded from this
repository. Reproducible input files, textual outputs, optimized structures,
analysis files, and scripts are retained.

## Repository structure

- `structures/` - starting molecular structures and construction notes
- `calculations/` - calculation inputs, textual outputs, and optimized geometries
- `analysis/` - processed results and final structures
- `scripts/` - calculation and analysis scripts
- `notes/` - research notes
- `logs/` - selected workflow logs
