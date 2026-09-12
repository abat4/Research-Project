#!/usr/bin/env python3
"""Build Step-1 structures for Fe/S and CODH C-cluster models.

Open-source dependencies:
    RDKit (used only to add standard-valence ligand hydrogens)

The script accepts either the full PDB entry 3B52 or the bundled verbatim
active-site record subset. No geometry optimization, DFT, or spin assignment
is performed.
"""
from __future__ import annotations

import argparse
import math
from dataclasses import dataclass, replace
from pathlib import Path
from collections import defaultdict, Counter

try:
    from rdkit import Chem
    from rdkit.Geometry import Point3D
except ImportError as exc:
    raise SystemExit("RDKit is required (open-source). Install e.g. conda-forge::rdkit.") from exc


@dataclass
class Atom:
    element: str
    name: str
    x: float
    y: float
    z: float
    resname: str = "MOL"
    chain: str = "X"
    resseq: int = 1
    record: str = "HETATM"
    altloc: str = ""
    occupancy: float = 1.0
    source: str = ""

    @property
    def xyz(self):
        return (self.x, self.y, self.z)


def distance(a: Atom, b: Atom) -> float:
    return math.dist(a.xyz, b.xyz)


def parse_pdb(path: Path) -> list[Atom]:
    atoms = []
    for line in path.read_text().splitlines():
        if not (line.startswith("ATOM  ") or line.startswith("HETATM")):
            continue
        # Standard PDB fixed columns. Bundled source lines are copied from 3B52.
        atoms.append(Atom(
            element=(line[76:78].strip() or ''.join(c for c in line[12:16] if c.isalpha())[:2]).title(),
            name=line[12:16].strip(),
            altloc=line[16:17].strip(),
            resname=line[17:20].strip(),
            chain=line[21:22].strip() or "X",
            resseq=int(line[22:26]),
            x=float(line[30:38]), y=float(line[38:46]), z=float(line[46:54]),
            occupancy=float(line[54:60]),
            record=line[:6].strip(),
            source=line,
        ))
    return atoms


def choose_highest_occupancy(records: list[Atom], context: str) -> Atom:
    if not records:
        raise ValueError(f"Missing required atom: {context}")
    # Deterministic, data-driven altloc policy. Ties are not silently resolved.
    max_occ = max(a.occupancy for a in records)
    top = [a for a in records if abs(a.occupancy - max_occ) < 1e-9]
    if len(top) > 1:
        coords = {tuple(round(v, 4) for v in a.xyz) for a in top}
        if len(coords) > 1:
            raise ValueError(f"Unresolved equal-occupancy alternate coordinates for {context}: {top}")
    return top[0]


def selected_residue_heavy(source: list[Atom], resname: str, resseq: int, names: list[str]) -> tuple[list[Atom], list[str]]:
    chosen = []
    notes = []
    for name in names:
        cand = [a for a in source if a.chain == "X" and a.resname == resname and a.resseq == resseq and a.name == name]
        if len(cand) > 1:
            desc = ", ".join(f"altloc {a.altloc or '-'} occ={a.occupancy:.2f}" for a in cand)
            notes.append(f"{resname}{resseq} {name}: alternates present ({desc}); highest occupancy selected.")
        a = choose_highest_occupancy(cand, f"{resname}{resseq} {name}")
        chosen.append(replace(a, altloc="", occupancy=1.0, record="HETATM"))
    return chosen, notes


def rdkit_fragment(smiles: str, heavy_atoms: list[Atom], heavy_names: list[str], resname: str, resseq: int) -> list[Atom]:
    """Add H coordinates without moving deposited heavy atoms."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumAtoms() != len(heavy_atoms):
        raise ValueError(f"SMILES/heavy-atom mismatch for {resname}{resseq}")
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, a in enumerate(heavy_atoms):
        conf.SetAtomPosition(i, Point3D(a.x, a.y, a.z))
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)
    mol_h = Chem.AddHs(mol, addCoords=True)
    conf_h = mol_h.GetConformer()

    out = []
    for i, a in enumerate(heavy_atoms):
        p = conf_h.GetAtomPosition(i)
        # RDKit does not move heavy atoms; assert rather than assume.
        if math.dist((p.x, p.y, p.z), a.xyz) > 1e-6:
            raise ValueError(f"RDKit moved heavy atom {resname}{resseq} {heavy_names[i]}")
        out.append(replace(a, name=heavy_names[i], record="HETATM", altloc="", occupancy=1.0))

    h_counts = defaultdict(int)
    for i in range(len(heavy_atoms), mol_h.GetNumAtoms()):
        atom = mol_h.GetAtomWithIdx(i)
        if atom.GetSymbol() != "H":
            continue
        nbr = atom.GetNeighbors()[0].GetIdx()
        parent = heavy_names[nbr]
        h_counts[parent] += 1
        hname = f"H{parent}{h_counts[parent]}"[:4]
        p = conf_h.GetAtomPosition(i)
        out.append(Atom("H", hname, p.x, p.y, p.z, resname, "X", resseq, "HETATM", source="RDKit AddHs"))
    return out


def build_model_a() -> tuple[list[Atom], list[tuple[int,int]], str]:
    # Two interpenetrating tetrahedra. Fe-Fe target 2.70 A; nearest Fe-S target 2.30 A.
    d_ff = 2.70
    d_fs = 2.30
    r_fe = d_ff * math.sqrt(6.0) / 4.0
    # Solve rS^2 -(2/3)rFe*rS + rFe^2 - dFS^2 = 0 (positive physical root).
    B = -2.0 * r_fe / 3.0
    C = r_fe*r_fe - d_fs*d_fs
    r_s = (-B + math.sqrt(B*B - 4*C)) / 2.0
    invsqrt3 = 1.0 / math.sqrt(3.0)
    fe_dirs = [(1,1,1), (1,-1,-1), (-1,1,-1), (-1,-1,1)]
    s_dirs  = [(-1,1,1), (1,-1,1), (1,1,-1), (-1,-1,-1)]

    atoms: list[Atom] = []
    bonds: list[tuple[int,int]] = []
    for i, d in enumerate(fe_dirs, 1):
        atoms.append(Atom("Fe", f"FE{i}", *(r_fe*c*invsqrt3 for c in d), "F4S", "A", 1))
    for i, d in enumerate(s_dirs, 1):
        atoms.append(Atom("S", f"S{i}", *(r_s*c*invsqrt3 for c in d), "F4S", "A", 1))
    # Core topology = all Fe/S pairs shorter than 2.6 A.
    for i in range(4):
        for j in range(4,8):
            if distance(atoms[i], atoms[j]) < 2.6:
                bonds.append((i,j))

    # One terminal methylthiolate per Fe, along outward Fe-center radial direction.
    for k in range(4):
        fe = atoms[k]
        norm = math.sqrt(fe.x**2 + fe.y**2 + fe.z**2)
        u = (fe.x/norm, fe.y/norm, fe.z/norm)
        sidx = len(atoms)
        sx, sy, sz = (fe.x + 2.28*u[0], fe.y + 2.28*u[1], fe.z + 2.28*u[2])
        atoms.append(Atom("S", f"SL{k+1}", sx, sy, sz, "MTS", "A", k+2))
        cidx = len(atoms)
        cx, cy, cz = (sx + 1.82*u[0], sy + 1.82*u[1], sz + 1.82*u[2])
        atoms.append(Atom("C", f"C{k+1}", cx, cy, cz, "MTS", "A", k+2))
        bonds.extend([(k,sidx),(sidx,cidx)])

        # Tetrahedral CH3 H directions around C-S bond.
        # v points from C toward S. H directions have dot(h,v)=-1/3.
        v = (-u[0], -u[1], -u[2])
        ref = (0.0,0.0,1.0) if abs(v[2]) < 0.9 else (0.0,1.0,0.0)
        e1raw = (v[1]*ref[2]-v[2]*ref[1], v[2]*ref[0]-v[0]*ref[2], v[0]*ref[1]-v[1]*ref[0])
        n1 = math.sqrt(sum(q*q for q in e1raw)); e1=tuple(q/n1 for q in e1raw)
        e2 = (v[1]*e1[2]-v[2]*e1[1], v[2]*e1[0]-v[0]*e1[2], v[0]*e1[1]-v[1]*e1[0])
        for h,phi in enumerate((0.0, 2*math.pi/3, 4*math.pi/3),1):
            perp = tuple(math.cos(phi)*e1[q] + math.sin(phi)*e2[q] for q in range(3))
            direction = tuple((-1/3)*v[q] + math.sqrt(8/9)*perp[q] for q in range(3))
            hidx=len(atoms)
            atoms.append(Atom("H", f"H{k+1}{h}", cx+1.09*direction[0], cy+1.09*direction[1], cz+1.09*direction[2], "MTS", "A", k+2))
            bonds.append((cidx,hidx))
    assumptions = (
        "Idealized cubane-like starting geometry generated from interpenetrating Fe/S tetrahedra; "
        "target Fe-Fe 2.70 A and nearest Fe-S 2.30 A. Four terminal Fe-S bonds are 2.28 A, "
        "S-C bonds 1.82 A, and C-H bonds 1.09 A. No optimization was performed."
    )
    return atoms, bonds, assumptions


def build_model_b(source: list[Atom]) -> tuple[list[Atom], list[str]]:
    notes: list[str] = []
    ligands: list[Atom] = []

    # Keep CA and side chain; delete backbone N/C/O. CA becomes a methyl cap after AddHs.
    # Cys -> CH3-CH2-S(-); His -> capped neutral histidine side chain, ND1-H / NE2 donor.
    for resseq in (295,333,446,476,526):
        heavy, n = selected_residue_heavy(source, "CYS", resseq, ["CA","CB","SG"])
        notes.extend(n)
        ligands.extend(rdkit_fragment("CC[S-]", heavy, ["CA","CB","SG"], "CYS", resseq))
    heavy, n = selected_residue_heavy(source, "HIS", 261, ["CA","CB","CG","ND1","CE1","NE2","CD2"])
    notes.extend(n)
    ligands.extend(rdkit_fragment("CCc1[nH]cnc1", heavy, ["CA","CB","CG","ND1","CE1","NE2","CD2"], "HIS", 261))

    # Current 3B52 XCC cluster. Select highest occupancy for each unique atom name.
    cluster = []
    for name in ("FE1","FE2","FE3","FE4","S1","S2","S3","S4","NI"):
        cand = [a for a in source if a.chain == "X" and a.resname == "XCC" and a.resseq == 1006 and a.name == name]
        if len(cand)>1:
            desc = ", ".join(f"altloc {a.altloc or '-'} occ={a.occupancy:.2f} at ({a.x:.3f},{a.y:.3f},{a.z:.3f})" for a in cand)
            notes.append(f"XCC {name}: alternates present ({desc}); highest occupancy selected.")
        a = choose_highest_occupancy(cand, f"XCC1006 {name}")
        cluster.append(replace(a, altloc="", occupancy=1.0, record="HETATM"))

    notes.append("CO2 X1005 is present in the deposited structure but intentionally excluded because the requested Model B composition lists the metal cluster plus six protein ligands, not substrate CO2.")
    notes.append("Cys295 is independently highest-occupancy altloc A (0.80), whereas XCC FE2 is highest-occupancy altloc B (0.70). This mixed major-occupancy selection is explicit; no alternate coordinates were averaged.")
    notes.append("The deposited LINK records associate Cys295(A)-FE2(A) and His261-FE2(B). Therefore no new Cys295(A)-FE2(B) bond record is invented; its starting distance is reported only as a geometric contact.")
    notes.append("Five cysteines are represented as deprotonated thiolates; His261 is neutral with ND1-H and NE2 unprotonated. This protonation choice is a starting-model assumption, not an experimentally resolved hydrogen assignment.")
    return ligands + cluster, notes


def build_model_c(model_b: list[Atom]) -> list[Atom]:
    out=[]
    changed=0
    for a in model_b:
        if a.resname == "XCC" and a.name == "NI" and a.element.lower() == "ni":
            out.append(replace(a, element="Fe", name="FE5"))
            changed += 1
        else:
            out.append(replace(a))
    if changed != 1:
        raise ValueError(f"Expected to replace exactly one Ni atom; replaced {changed}")
    return out


def write_xyz(path: Path, atoms: list[Atom], comment: str):
    with path.open("w") as f:
        f.write(f"{len(atoms)}\n{comment}\n")
        for a in atoms:
            f.write(f"{a.element:2s} {a.x:14.8f} {a.y:14.8f} {a.z:14.8f}\n")


def pdb_atom_line(serial: int, a: Atom) -> str:
    rec = (a.record or "HETATM")[:6]
    # Right-justify short atom names for readable PDB; element kept explicit.
    atomname = a.name[:4].rjust(4)
    return (f"{rec:<6}{serial:5d} {atomname} {a.resname[:3]:>3s} {a.chain[:1] or 'X'}{a.resseq:4d}    "
            f"{a.x:8.3f}{a.y:8.3f}{a.z:8.3f}{1.00:6.2f}{0.00:6.2f}          {a.element.upper():>2s}")


def write_pdb(path: Path, atoms: list[Atom], remarks: list[str], bonds: list[tuple[int,int]] | None=None):
    with path.open("w") as f:
        for r in remarks:
            # Wrap simply to PDB REMARK width.
            while r:
                f.write("REMARK 999 " + r[:67] + "\n")
                r=r[67:]
        for i,a in enumerate(atoms,1):
            f.write(pdb_atom_line(i,a)+"\n")
        if bonds:
            for i,j in bonds:
                f.write(f"CONECT{i+1:5d}{j+1:5d}\n")
        f.write("END\n")


def element_counts(atoms: list[Atom]) -> Counter:
    return Counter(a.element.capitalize() for a in atoms)


def validate(atoms: list[Atom], model: str, bonds: list[tuple[int,int]] | None=None) -> list[str]:
    report=[]
    # Duplicate/near-duplicate check.
    closest=(1e9,None,None)
    for i in range(len(atoms)):
        for j in range(i+1,len(atoms)):
            d=distance(atoms[i],atoms[j])
            if d<closest[0]: closest=(d,i,j)
    report.append(f"Closest atom-atom separation: {closest[0]:.3f} A ({atoms[closest[1]].name}-{atoms[closest[2]].name}).")
    if closest[0] < 0.55:
        raise ValueError(f"{model}: probable duplicate atoms ({closest[0]:.3f} A)")

    counts=element_counts(atoms)
    if model=="A":
        if counts["Fe"]!=4 or counts["S"]!=8: raise ValueError(f"Model A element count wrong: {counts}")
    elif model=="B":
        if counts["Fe"]!=4 or counts["Ni"]!=1 or counts["S"]!=9: raise ValueError(f"Model B element count wrong: {counts}")
    elif model=="C":
        if counts["Fe"]!=5 or counts["Ni"]!=0 or counts["S"]!=9: raise ValueError(f"Model C element count wrong: {counts}")
    report.append(f"Element counts: {dict(sorted(counts.items()))}")

    # Bonded H length checks when explicit graph is known, else infer nearest heavy atom.
    if bonds:
        h_lengths=[]
        for i,j in bonds:
            if atoms[i].element=="H" or atoms[j].element=="H": h_lengths.append(distance(atoms[i],atoms[j]))
        if h_lengths:
            report.append(f"Explicit X-H bond range: {min(h_lengths):.3f}-{max(h_lengths):.3f} A.")
            if min(h_lengths)<0.80 or max(h_lengths)>1.25: raise ValueError(f"{model}: unrealistic explicit H bond length")
    else:
        h_nearest=[]
        heavy=[a for a in atoms if a.element!="H"]
        for h in [a for a in atoms if a.element=="H"]:
            h_nearest.append(min(distance(h,x) for x in heavy))
        if h_nearest:
            report.append(f"Nearest-heavy distance for added H atoms: {min(h_nearest):.3f}-{max(h_nearest):.3f} A.")
            if min(h_nearest)<0.75 or max(h_nearest)>1.25: raise ValueError(f"{model}: implausible added-H nearest-heavy distance")
        second_heavy=[]
        for h in [a for a in atoms if a.element=="H"]:
            ds=sorted(distance(h,x) for x in heavy)
            if len(ds)>1: second_heavy.append(ds[1])
        if second_heavy:
            report.append(f"Closest H--nonparent-heavy separation (second-nearest heavy): {min(second_heavy):.3f} A.")
            if min(second_heavy)<1.35: raise ValueError(f"{model}: possible H/nonparent-heavy clash")

    return report


def named(atoms, resseq=None, resname=None, name=None, element=None):
    out=[]
    for a in atoms:
        if resseq is not None and a.resseq!=resseq: continue
        if resname is not None and a.resname!=resname: continue
        if name is not None and a.name!=name: continue
        if element is not None and a.element.lower()!=element.lower(): continue
        out.append(a)
    if len(out)!=1:
        raise ValueError(f"Expected one atom for {resname}{resseq} {name}/{element}; got {len(out)}")
    return out[0]


def model_b_distance_report(atoms: list[Atom]) -> list[str]:
    pairs=[
        ("His261 NE2--FE2(B selected)", named(atoms,261,"HIS","NE2"), named(atoms,1006,"XCC","FE2")),
        ("Cys295(A selected) SG--FE2(B selected) [geometric contact; not deposited LINK]", named(atoms,295,"CYS","SG"), named(atoms,1006,"XCC","FE2")),
        ("Cys333 SG--FE1", named(atoms,333,"CYS","SG"), named(atoms,1006,"XCC","FE1")),
        ("Cys446 SG--FE3", named(atoms,446,"CYS","SG"), named(atoms,1006,"XCC","FE3")),
        ("Cys476 SG--FE4", named(atoms,476,"CYS","SG"), named(atoms,1006,"XCC","FE4")),
        ("Cys526 SG--Ni", named(atoms,526,"CYS","SG"), named(atoms,1006,"XCC","NI")),
    ]
    out=[f"{label}: {distance(a,b):.3f} A" for label,a,b in pairs]
    metals=[a for a in atoms if a.resname=="XCC" and a.element in ("Fe","Ni")]
    core_s=[a for a in atoms if a.resname=="XCC" and a.element=="S"]
    contacts=[]
    for m in metals:
        for s in core_s:
            d=distance(m,s)
            if d < 2.75:
                contacts.append((m.name,s.name,d))
    out.append("Metal--core-S contacts under 2.75 A: " + ", ".join(f"{m}-{s} {d:.2f}" for m,s,d in contacts))
    return out


def notes_text(title, atoms, metals, ligands, intended_charge, assumptions, validation):
    lines=[
        title,
        "="*len(title),
        f"Number of atoms: {len(atoms)}",
        f"Metal atoms present: {metals}",
        f"Ligand atoms/fragments: {ligands}",
        f"Overall intended charge: {intended_charge}",
        "",
        "Assumptions / modifications:",
    ]
    for x in assumptions: lines.append(f"- {x}")
    lines += ["", "Validation checks:"]
    for x in validation: lines.append(f"- {x}")
    return "\n".join(lines)+"\n"


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--pdb", type=Path, default=Path(__file__).parent/"source_data"/"3B52_active_site_records.pdb",
                    help="Full 3B52 PDB file or bundled active-site record subset")
    ap.add_argument("--out", type=Path, default=Path(__file__).parent)
    args=ap.parse_args()
    source=parse_pdb(args.pdb)
    root=args.out
    A_dir=root/"model_A_Fe4S4"; B_dir=root/"model_B_CODH_NiFeS"; C_dir=root/"model_C_CODH_Fe_control"
    for d in (A_dir,B_dir,C_dir): d.mkdir(parents=True,exist_ok=True)

    A, A_bonds, A_assumption=build_model_a()
    vA=validate(A,"A",A_bonds)
    # More structure-specific checks.
    coreFe=A[:4]; coreS=A[4:8]
    fs=[distance(f,s) for f in coreFe for s in coreS if distance(f,s)<2.6]
    ff=[distance(coreFe[i],coreFe[j]) for i in range(4) for j in range(i+1,4)]
    vA.append(f"Cubane Fe-S edge range: {min(fs):.3f}-{max(fs):.3f} A ({len(fs)} edges).")
    vA.append(f"Cubane Fe-Fe range: {min(ff):.3f}-{max(ff):.3f} A.")
    write_xyz(A_dir/"model_A_Fe4S4.xyz",A,"[Fe4S4(SCH3)4]2- idealized starting structure; no optimization")
    write_pdb(A_dir/"model_A_Fe4S4.pdb",A,["Model A: idealized [Fe4S4(SCH3)4]2- starting structure."],A_bonds)
    (A_dir/"model_A_info.txt").write_text(notes_text(
        "Model A: Simple [4Fe-4S] cluster", A,
        "4 Fe", "4 bridging sulfides + 4 methylthiolate S/C fragments (SCH3)", "-2",
        [A_assumption,"Formal charge follows the requested [Fe4S4(SCH3)4]2- composition; no spin state assigned."],vA))

    B, bnotes=build_model_b(source)
    vB=validate(B,"B")
    vB.extend(model_b_distance_report(B))
    write_xyz(B_dir/"model_B_CODH_NiFeS.xyz",B,"3B52 CODH C-cluster first-shell model; major-occupancy coordinate selection; no optimization")
    write_pdb(B_dir/"model_B_CODH_NiFeS.pdb",B,[
        "Model B: truncated 3B52 C-cluster first shell.",
        "Major-occupancy coordinate selection; no coordinate averaging.",
        "No CONECT/LINK invented for ambiguous Cys295(A)-FE2(B).",
    ])
    (B_dir/"model_B_info.txt").write_text(notes_text(
        "Model B: CODH Ni-Fe-S active-site model",B,
        "4 Fe + 1 Ni", "Cys295/Cys333/Cys446/Cys476/Cys526 truncated at CA; His261 truncated at CA; 4 inorganic cluster sulfides",
        "Not assigned in Step 1 (depends on C-cluster redox state/protonation; five cysteine fragments carry -1 each in the ligand graph)",
        bnotes + ["Backbone N/C/O atoms were removed. The retained CA atoms were saturated with H by RDKit, making each CA a methyl cap; side-chain heavy-atom coordinates remain exactly experimental."],vB))

    C=build_model_c(B)
    vC=validate(C,"C")
    # Confirm coordinate identity to B atom-by-atom.
    maxdelta=max(math.dist(a.xyz,c.xyz) for a,c in zip(B,C))
    vC.append(f"Maximum coordinate displacement relative to Model B: {maxdelta:.6f} A (must be 0).")
    fe5=named(C,1006,"XCC","FE5")
    c526=named(C,526,"CYS","SG")
    vC.append(f"Artificial FE5--Cys526 SG starting distance: {distance(fe5,c526):.3f} A; this short distance is intentionally inherited from the experimental Ni--S coordinate and was not relaxed.")
    write_xyz(C_dir/"model_C_CODH_Fe_control.xyz",C,"Artificial control: Model B with Ni->Fe at identical coordinates; no optimization")
    write_pdb(C_dir/"model_C_CODH_Fe_control.pdb",C,[
        "Model C: artificial Fe-substituted control.",
        "Exact coordinates copied from Model B; XCC NI renamed FE5 and element set to Fe.",
    ])
    (C_dir/"model_C_info.txt").write_text(notes_text(
        "Model C: Fe-substituted CODH control",C,
        "5 Fe (Model B Ni site replaced by Fe)", "Identical ligand atoms and coordinates to Model B",
        "Not assigned in Step 1; substitution changes element identity only, not a formally specified redox/protonation state",
        ["Exact copy of Model B coordinates followed by one atom-identity substitution: Ni -> Fe (named FE5).","No geometry relaxation, spin assignment, or DFT calculation was performed."],vC))

    # Machine-readable validation summary.
    report=["STRUCTURE VALIDATION SUMMARY","", "MODEL A"]+vA+["","MODEL B"]+vB+["","MODEL C"]+vC
    (root/"validation_report.txt").write_text("\n".join(report)+"\n")

    print("Built:")
    print(f"  Model A: {len(A)} atoms")
    print(f"  Model B: {len(B)} atoms")
    print(f"  Model C: {len(C)} atoms")
    print("Validation completed without fatal geometry/duplicate/H-placement errors.")

if __name__ == "__main__":
    main()
