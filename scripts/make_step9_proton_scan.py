import math
from pathlib import Path

R = Path.home()/"Research-Project"
OUT = R/"calculations"/"step9_protonation"
OUT.mkdir(parents=True, exist_ok=True)

cases = {
    "B": {
        "xyz": R/"analysis"/"model_B_CO2_xtb_qm4_reduced.xyz",
        "charges": R/"calculations"/"step4_5_co2"/"B"/"charges"
    },
    "C": {
        "xyz": R/"analysis"/"model_C_CO2_xtb_qm4_reduced.xyz",
        "charges": R/"calculations"/"step4_5_co2"/"C"/"charges"
    }
}

def read_xyz(p):
    L = p.read_text().splitlines()
    n = int(L[0])
    atoms = []
    for line in L[2:2+n]:
        s = line.split()
        atoms.append([s[0], *map(float, s[1:4])])
    return atoms

def read_charges(p):
    q = []
    for line in p.read_text().splitlines():
        try:
            q.append(float(line.split()[0]))
        except:
            pass
    return q

for model, data in cases.items():
    atoms = read_xyz(data["xyz"])
    q = read_charges(data["charges"])

    # Last 3 atoms = C, O, O
    C = atoms[-3]
    O1 = atoms[-2]
    O2 = atoms[-1]

    # Choose more negatively charged oxygen
    if q[-2] <= q[-1]:
        O = O1
        Oidx = len(atoms)-1
        oq = q[-2]
    else:
        O = O2
        Oidx = len(atoms)
        oq = q[-1]

    c = C[1:]
    o = O[1:]

    # Direction outward along C-O bond
    v = [o[i]-c[i] for i in range(3)]
    norm = math.sqrt(sum(x*x for x in v))
    v = [x/norm for x in v]

    D = OUT/model
    D.mkdir(parents=True, exist_ok=True)

    print(f"{model}: protonating O atom {Oidx}, q(O)={oq:.4f}")

    for r in [3.00, 2.50, 2.00, 1.60, 1.30, 1.00]:
        H = [o[i] + r*v[i] for i in range(3)]
        new = atoms + [["H", *H]]

        p = D/f"OH_{r:.2f}.xyz"

        with open(p, "w") as f:
            f.write(f"{len(new)}\n")
            f.write(f"{model} proton approach; O-H={r:.2f} A\n")
            for a,x,y,z in new:
                f.write(
                    f"{a:2s} {x:15.8f} {y:15.8f} {z:15.8f}\n"
                )

print("Created all Step 9 structures.")
