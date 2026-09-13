from pathlib import Path
import csv

R = Path.home()/"Research-Project"

cases = [
    {
        "model":"B",
        "state":"baseline_-3",
        "dir":R/"calculations"/"step4_5_co2"/"qm3"/"B",
        "xyz":R/"analysis"/"model_B_CO2_xtb_qm3_baseline.xyz"
    },
    {
        "model":"C",
        "state":"baseline_-3",
        "dir":R/"calculations"/"step4_5_co2"/"qm3"/"C",
        "xyz":R/"analysis"/"model_C_CO2_xtb_qm3_baseline.xyz"
    },
    {
        "model":"B",
        "state":"reduced_-4",
        "dir":R/"calculations"/"step4_5_co2"/"B",
        "xyz":R/"analysis"/"model_B_CO2_xtb_qm4_reduced.xyz"
    },
    {
        "model":"C",
        "state":"reduced_-4",
        "dir":R/"calculations"/"step4_5_co2"/"C",
        "xyz":R/"analysis"/"model_C_CO2_xtb_qm4_reduced.xyz"
    }
]

def read_xyz(p):
    L=p.read_text().splitlines()
    n=int(L[0])
    A=[]
    for line in L[2:2+n]:
        s=line.split()
        A.append((s[0],*[float(x) for x in s[1:4]]))
    return A

def read_charges(p):
    q=[]
    for line in p.read_text().splitlines():
        try:
            q.append(float(line.split()[0]))
        except:
            pass
    return q

def read_wbo(p):
    W={}
    for line in p.read_text().splitlines():
        s=line.split()
        if len(s) < 3:
            continue
        try:
            i=int(s[0])
            j=int(s[1])
            v=float(s[2])
        except:
            continue
        W[tuple(sorted((i,j)))] = v
    return W

def bond(W,i,j):
    return W.get(tuple(sorted((i,j))),0.0)

rows=[]

for case in cases:
    A=read_xyz(case["xyz"])
    n=len(A)

    q=read_charges(case["dir"]/"charges")
    W=read_wbo(case["dir"]/"wbo")

    # B/C cluster = first 64 atoms; CO2 = atoms 65-67.
    site=64
    C=65
    O1=66
    O2=67

    rows.append({
        "model":case["model"],
        "state":case["state"],
        "q_CO2":sum(q[-3:]),
        "q_site_metal":q[site-1],
        "WBO_site_C":bond(W,site,C),
        "WBO_C_O1":bond(W,C,O1),
        "WBO_C_O2":bond(W,C,O2)
    })

out=R/"analysis"/"step6_charge_bond_order.csv"

with open(out,"w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys())
    w.writeheader()
    w.writerows(rows)

print("===== STEP 6 CHARGE / BOND ORDER =====")
for r in rows:
    print()
    print(r["model"],r["state"])
    print("q(CO2)       =",round(r["q_CO2"],4))
    print("q(site metal)=",round(r["q_site_metal"],4))
    print("WBO M-C      =",round(r["WBO_site_C"],4))
    print("WBO C-O1     =",round(r["WBO_C_O1"],4))
    print("WBO C-O2     =",round(r["WBO_C_O2"],4))

print()
print("Saved:",out)
