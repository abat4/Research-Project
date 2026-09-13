from pathlib import Path
import csv

R = Path.home()/"Research-Project"
OUT = R/"analysis"

rows = [
    {
        "model":"B_Ni",
        "state":"baseline_-3",
        "metal_C_A":1.8362,
        "metal_C_WBO":0.6792,
        "CO2_or_COOH_charge":-0.6803,
        "OCO_deg":132.45,
        "CO1_A":1.2305,
        "CO2_A":1.2387,
        "frontier_gap_eV":0.2954198
    },
    {
        "model":"C_Fe",
        "state":"baseline_-3",
        "metal_C_A":3.7643,
        "metal_C_WBO":0.0000,
        "CO2_or_COOH_charge":-0.7633,
        "OCO_deg":131.45,
        "CO1_A":1.2319,
        "CO2_A":1.2331,
        "frontier_gap_eV":0.1634122
    },
    {
        "model":"B_Ni",
        "state":"reduced_-4",
        "metal_C_A":3.5052,
        "metal_C_WBO":0.0000,
        "CO2_or_COOH_charge":-0.8492,
        "OCO_deg":129.55,
        "CO1_A":1.2403,
        "CO2_A":1.2394,
        "frontier_gap_eV":0.1641049
    },
    {
        "model":"C_Fe",
        "state":"reduced_-4",
        "metal_C_A":3.5408,
        "metal_C_WBO":0.0000,
        "CO2_or_COOH_charge":-0.8391,
        "OCO_deg":129.97,
        "CO1_A":1.2394,
        "CO2_A":1.2400,
        "frontier_gap_eV":0.1865364
    },
    {
        "model":"B_Ni",
        "state":"COOH_-3",
        "metal_C_A":3.6344,
        "metal_C_WBO":0.0000,
        "CO2_or_COOH_charge":-0.3231,
        "OCO_deg":118.84,
        "CO1_A":1.3764,
        "CO2_A":1.2283,
        "frontier_gap_eV":""
    },
    {
        "model":"C_Fe",
        "state":"COOH_-3",
        "metal_C_A":3.7793,
        "metal_C_WBO":0.0000,
        "CO2_or_COOH_charge":-0.3277,
        "OCO_deg":118.76,
        "CO1_A":1.2318,
        "CO2_A":1.3769,
        "frontier_gap_eV":""
    }
]

p = OUT/"final_results_summary.csv"

with open(p,"w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys())
    w.writeheader()
    w.writerows(rows)

print("Wrote:",p)

print("\n===== FINAL B/C SUMMARY =====")
for r in rows:
    print(
        f'{r["model"]:5s} {r["state"]:12s} '
        f'M-C={r["metal_C_A"]:6.3f} A   '
        f'WBO={r["metal_C_WBO"]:5.3f}   '
        f'q={r["CO2_or_COOH_charge"]:7.3f}   '
        f'OCO={r["OCO_deg"]:6.2f}'
    )

print("\n===== MATCHED-GEOMETRY INTERACTION =====")
print("Ni =  -93.82 kJ/mol")
print("Fe = -124.16 kJ/mol")
print("Fe more favorable by 30.34 kJ/mol")

print("\n===== BASELINE STRUCTURAL REORGANIZATION =====")
print("B cluster RMSD = 0.6202 A")
print("B max displacement = 1.0583 A")
print("C cluster RMSD = 0.8499 A")
print("C max displacement = 2.5105 A")
