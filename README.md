# Step 1 starting structures: Fe–S / CODH C-cluster models

No DFT, geometry optimization, or spin-state assignment is included here. The files are starting geometries only.

| Model | Composition | Purpose |
|---|---|---|
| A | `[Fe4S4(SCH3)4]2-` | Simple Fe–S benchmark |
| B | CODH `NiFe4S4` C-cluster + Cys295/Cys333/Cys446/Cys476/Cys526/His261 truncated first shell | Biological CO2-activating center |
| C | Model B with Ni replaced by Fe | Test the specific effect of Ni |

## Directory layout

```text
structures/
    README.md
    build_structures.py
    validation_report.txt
    requirements.txt
    source_data/
        3B52_active_site_records.pdb
        SOURCE.txt
    model_A_Fe4S4/
        model_A_Fe4S4.xyz
        model_A_Fe4S4.pdb
        model_A_info.txt
    model_B_CODH_NiFeS/
        model_B_CODH_NiFeS.xyz
        model_B_CODH_NiFeS.pdb
        model_B_info.txt
    model_C_CODH_Fe_control/
        model_C_CODH_Fe_control.xyz
        model_C_CODH_Fe_control.pdb
        model_C_info.txt
```

## Reproduction

The builder uses Python plus the open-source RDKit package only. It can parse either the bundled active-site record subset or a full 3B52 PDB file.

```bash
python build_structures.py
# or, using a separately downloaded full 3B52 PDB:
python build_structures.py --pdb /path/to/3B52.pdb --out /path/to/structures
```

The bundled `source_data/3B52_active_site_records.pdb` contains the verbatim coordinate/LINK records needed for this construction, including both Cys295 and FE2 alternate locations and the deposited CO2 records for transparency. CO2 itself is not included in Model B because the requested composition specified the metal cluster plus the six protein ligands.

## Important 3B52 ambiguity

The current remediated 3B52 entry contains FE2 alternate positions: A (occupancy 0.10) and B (0.70). Cys295 also has A/B alternatives, with A at occupancy 0.80. The builder uses a deterministic highest-occupancy choice for each required atom, so Model B contains FE2(B) and Cys295(A); it does not average coordinates. Because the deposited LINK records specifically associate Cys295(A) with FE2(A) and His261 with FE2(B), the generated PDB does **not** invent a Cys295(A)–FE2(B) LINK. Their 2.337 Å separation is reported as a geometric contact only.

Five cysteine ligands are represented as thiolates. His261 is represented as neutral ND1-H / NE2-donor histidine. The total electronic charge for Models B/C is intentionally left unassigned in Step 1 because the PDB does not uniquely define the C-cluster redox/protonation state; this should be fixed together with the electronic-state setup before DFT.

## Validation

`validation_report.txt` checks atom counts, missing metals/sulfur, near-duplicate atoms, representative metal–ligand distances, H–heavy-atom distances, and short H/nonparent-heavy contacts. Model C is verified to have exactly zero coordinate displacement from Model B before/after the Ni→Fe identity replacement.
