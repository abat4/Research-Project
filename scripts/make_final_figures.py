from pathlib import Path
import matplotlib.pyplot as plt

R = Path.home()/"Research-Project"
F = R/"analysis"/"figures"
F.mkdir(parents=True,exist_ok=True)

# ---------- Figure 1: baseline metal-CO2 binding ----------
labels=["Ni model","Fe control"]

dist=[1.8362,3.7643]
wbo=[0.6792,0.0]

fig,ax=plt.subplots(figsize=(6,4))
ax.bar(labels,dist)
ax.set_ylabel("Metal–C distance (Å)")
ax.set_title("Optimized baseline CO₂ coordination")
fig.tight_layout()
fig.savefig(F/"01_metal_C_distance.png",dpi=300)
plt.close(fig)

fig,ax=plt.subplots(figsize=(6,4))
ax.bar(labels,wbo)
ax.set_ylabel("Metal–C Wiberg bond order")
ax.set_title("Baseline metal–CO₂ bonding")
fig.tight_layout()
fig.savefig(F/"02_metal_C_WBO.png",dpi=300)
plt.close(fig)

# ---------- Figure 2: charge transfer ----------
baseline=[-0.6803,-0.7633]
reduced=[-0.8492,-0.8391]

x=[0,1]
width=0.35

fig,ax=plt.subplots(figsize=(6,4))
ax.bar([i-width/2 for i in x],baseline,width,label="Baseline (-3)")
ax.bar([i+width/2 for i in x],reduced,width,label="Reduced (-4)")
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.set_ylabel("Charge on CO₂ (e)")
ax.set_title("Redox response of CO₂")
ax.legend()
fig.tight_layout()
fig.savefig(F/"03_CO2_charge_redox.png",dpi=300)
plt.close(fig)

# ---------- Figure 3: frontier separation ----------
baseline_gap=[0.2954198,0.1634122]
reduced_gap=[0.1641049,0.1865364]

fig,ax=plt.subplots(figsize=(6,4))
ax.bar([i-width/2 for i in x],baseline_gap,width,label="Baseline (-3)")
ax.bar([i+width/2 for i in x],reduced_gap,width,label="Reduced (-4)")
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.set_ylabel("Frontier-orbital separation (eV)")
ax.set_title("Redox sensitivity")
ax.legend()
fig.tight_layout()
fig.savefig(F/"04_frontier_gap.png",dpi=300)
plt.close(fig)

# ---------- Figure 4: proton approach ----------
r=[3.00,2.50,2.00,1.60,1.30,1.00]
B=[0.00,15.88,51.42,118.61,132.14,16.38]
C=[0.00,16.01,51.80,117.66,125.36,5.11]

fig,ax=plt.subplots(figsize=(6,4))
ax.plot(r,B,marker="o",label="Ni model")
ax.plot(r,C,marker="o",label="Fe control")
ax.invert_xaxis()
ax.set_xlabel("O–H approach distance (Å)")
ax.set_ylabel("Relative electronic energy (kJ/mol)")
ax.set_title("Unrelaxed proton-approach screening")
ax.legend()
fig.tight_layout()
fig.savefig(F/"05_proton_approach.png",dpi=300)
plt.close(fig)

print("Figures written to:",F)
for p in sorted(F.glob("*.png")):
    print(p.name)
