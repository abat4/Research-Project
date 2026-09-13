import numpy as np
from pathlib import Path

R = Path.home()/"Research-Project"
S = R/"structures"
A = R/"analysis"
OUT = R/"calculations"/"step4_5_co2"

def readxyz(p):
    x = Path(p).read_text().splitlines()
    n = int(x[0])
    atoms = []
    for line in x[2:2+n]:
        q=line.split()
        atoms.append([q[0], *map(float,q[1:4])])
    return atoms

def writexyz(p, atoms, comment):
    with open(p,"w") as f:
        f.write(f"{len(atoms)}\n{comment}\n")
        for a,x,y,z in atoms:
            f.write(f"{a:2s} {x:15.8f} {y:15.8f} {z:15.8f}\n")

def arr(atoms):
    return np.array([[x,y,z] for _,x,y,z in atoms])

def kabsch(P,Q,X):
    pc=P.mean(0)
    qc=Q.mean(0)
    H=(P-pc).T@(Q-qc)
    U,s,Vt=np.linalg.svd(H)
    Rm=Vt.T@U.T
    if np.linalg.det(Rm)<0:
        Vt[-1]*=-1
        Rm=Vt.T@U.T
    return (X-pc)@Rm.T+qc

startB=readxyz(S/"model_B_CODH_NiFeS"/"model_B_CODH_NiFeS.xyz")
optB=readxyz(A/"model_B_xtb_optimized.xyz")
optC=readxyz(A/"model_C_xtb_optimized.xyz")
optA=readxyz(A/"model_A_xtb_optimized.xyz")

co2=[]
for line in (OUT/"3B52.pdb").read_text().splitlines():
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

co2=sorted(co2,key=lambda x: 0 if x[0].upper()=="C" else 1)

if len(co2)!=3:
    raise RuntimeError(f"Expected 3 CO2 atoms from 3B52, found {len(co2)}")

idx=[i for i,a in enumerate(startB) if a[0] in ("Fe","Ni","S")]
P=arr(startB)[idx]
X=np.array([[a[1],a[2],a[3]] for a in co2])

for name,target in [("B",optB),("C",optC)]:
    Q=arr(target)[idx]
    T=kabsch(P,Q,X)
    new=target+[[co2[i][0],*T[i]] for i in range(3)]
    writexyz(
        OUT/f"{name}_CO2_start.xyz",
        new,
        f"Model {name} + CO2; pose transferred from PDB 3B52"
    )

xyz=arr(optA)
met=[i for i,a in enumerate(optA) if a[0] in ("Fe","S")]
fe=next(i for i,a in enumerate(optA) if a[0]=="Fe")

center=xyz[met].mean(0)
u=xyz[fe]-center
u=u/np.linalg.norm(u)

ref=np.array([1.,0.,0.])
if abs(np.dot(u,ref))>0.9:
    ref=np.array([0.,1.,0.])

v=np.cross(u,ref)
v=v/np.linalg.norm(v)

C=xyz[fe]+2.05*u
O1=C+1.16*v
O2=C-1.16*v

new=optA+[["C",*C],["O",*O1],["O",*O2]]

writexyz(
    OUT/"A_CO2_start.xyz",
    new,
    "Model A + CO2; terminal Fe-C screening pose"
)

print("Created A_CO2_start.xyz, B_CO2_start.xyz, C_CO2_start.xyz")
