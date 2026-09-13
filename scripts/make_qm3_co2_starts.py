import numpy as np
from pathlib import Path

R=Path.home()/"Research-Project"
OUT=R/"calculations"/"step4_5_co2"/"qm3"
OUT.mkdir(parents=True,exist_ok=True)

def readxyz(p):
    L=Path(p).read_text().splitlines()
    n=int(L[0])
    A=[]
    for line in L[2:2+n]:
        s=line.split()
        A.append([s[0],*map(float,s[1:4])])
    return A

def xyzarr(A):
    return np.array([[x,y,z] for _,x,y,z in A])

def write(p,A,comment):
    with open(p,"w") as f:
        f.write(f"{len(A)}\n{comment}\n")
        for a,x,y,z in A:
            f.write(f"{a:2s} {x:15.8f} {y:15.8f} {z:15.8f}\n")

def transform(P,Q,X):
    pc=P.mean(0)
    qc=Q.mean(0)
    H=(P-pc).T@(Q-qc)
    U,s,Vt=np.linalg.svd(H)
    M=Vt.T@U.T
    if np.linalg.det(M)<0:
        Vt[-1]*=-1
        M=Vt.T@U.T
    return (X-pc)@M.T+qc

start=readxyz(
    R/"structures"/"model_B_CODH_NiFeS"/"model_B_CODH_NiFeS.xyz"
)

B=readxyz(R/"analysis"/"model_B_xtb_qm3_baseline.xyz")
C=readxyz(R/"analysis"/"model_C_xtb_qm3_baseline.xyz")

co2=[]
pdb=R/"calculations"/"step4_5_co2"/"3B52.pdb"

for line in pdb.read_text().splitlines():
    if line.startswith(("ATOM  ","HETATM")) and line[17:20].strip()=="CO2":
        el=line[76:78].strip()
        if not el:
            el=line[12:16].strip()[0]
        co2.append([
            el,
            float(line[30:38]),
            float(line[38:46]),
            float(line[46:54])
        ])

co2=sorted(co2,key=lambda x:0 if x[0].upper()=="C" else 1)

if len(co2)!=3:
    raise RuntimeError(f"Found {len(co2)} CO2 atoms, expected 3")

idx=[i for i,a in enumerate(start) if a[0] in ("Fe","Ni","S")]
P=xyzarr(start)[idx]
X=np.array([[a[1],a[2],a[3]] for a in co2])

for name,target in [("B",B),("C",C)]:
    Q=xyzarr(target)[idx]
    T=transform(P,Q,X)

    complex_atoms=target+[
        [co2[i][0],*T[i]] for i in range(3)
    ]

    write(
        OUT/f"{name}_qm3_CO2_start.xyz",
        complex_atoms,
        f"{name} charge -3 CO2 pose transferred from 3B52"
    )

    print(name,"created:",OUT/f"{name}_qm3_CO2_start.xyz")
