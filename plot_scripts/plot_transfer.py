"""Verification plot for the field transfer (BOUT++ elm-pb P -> PCMS 3D index-space field -> 2x finer mesh).

Slice shown: iy_file = 32 (outboard midplane, where |P| peaks), axes (ix, iz).
Data written by the transfer test (run inside the BOUT++ pause at output 1, t = 1):
  transfer_slice.csv        ix,iz,P_bout            the SOURCE slice, 64 x 16 nodes, from BOUT++ memory
  transfer_fine_slice.csv   ix_real,iz_real,P_pcms  the TRANSFERRED field on the 2x mesh, same row (127 x 33 nodes; iz = 16 is the periodic copy of iz = 0)
  transfer_targets.csv      xi_ix,zeta_iz,P_bout,P_pcms,...   1024 fractional fine nodes of that row: BOUT++'s XZBilinear of the original field, and the PCMS value
"""
import plots_dir  # data in, plots out: see plots_dir.py
import csv
import numpy as np
import matplotlib.pyplot as plt

BASE = plots_dir.DATA
rd = lambda path: list(csv.DictReader(open(path)))
col = lambda rows, k: np.array([float(r[k]) for r in rows])

sl = rd(f"{BASE}/transfer_slice.csv")
ix, iz, P = col(sl, "ix").astype(int), col(sl, "iz").astype(int), col(sl, "P_bout")
fn = rd(f"{BASE}/transfer_fine_slice.csv")
fx, fz, fp = col(fn, "ix_real"), col(fn, "iz_real"), col(fn, "P_pcms")
tg = rd(f"{BASE}/transfer_targets.csv")
pb, pp = col(tg, "P_bout"), col(tg, "P_pcms")
rt = rd(f"{BASE}/transfer_random_targets.csv")   # random real (ix, iz) in the same row: general fractions
rb, rp = col(rt, "P_bout"), col(rt, "P_pcms")

m = np.abs(P).max()                       # normalisation: max |P| over the source slice
print(f"source slice: n={len(P)}  min {P.min():.4e}  max {P.max():.4e}  max|P| = {m:.4e}")
print(f"fine slice:   n={len(fp)}  min {fp.min()/m:+.4f}  max {fp.max()/m:+.4f} (normalised)")

# the source slice as a 2D array, rows = iz, columns = ix - 2
src2d = np.full((iz.max() + 1, ix.max() - ix.min() + 1), np.nan)
src2d[iz, ix - ix.min()] = P

# fine nodes that coincide with source nodes: compare with the stored BOUT++ value (the iz = 16 plane is the periodic copy of iz = 0)
coin = (np.abs(fx - np.round(fx)) < 1e-12) & (np.abs(fz - np.round(fz)) < 1e-12)
orig = np.array([src2d[int(round(z)) % 16, int(round(x)) - ix.min()] for x, z in zip(fx[coin], fz[coin])])
dc = (fp[coin] - orig) / m
print(f"coincident nodes: n={int(coin.sum())}  PCMS identical to the original BOUT++ value at {int((dc == 0).sum())};  max|PCMS - BOUT++|/max|P| = {np.abs(dc).max():.3e}")

# fractional targets: BOUT++ XZBilinear of the original field
dt = (pp - pb) / m
print(f"fractional targets: n={len(pb)}  PCMS identical to BOUT++ XZBilinear at {int((pp == pb).sum())};  max|PCMS - BOUT++|/max|P| = {np.abs(dt).max():.3e}")

# the fine slice as a 2D array, rows = iz (step 0.5), columns = ix (step 0.5)
fine2d = np.full((33, 127), np.nan)
fine2d[np.round(fz / 0.5).astype(int), np.round((fx - ix.min()) / 0.5).astype(int)] = fp / m

fig, axs = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
ext_s = [ix.min() - 0.5, ix.max() + 0.5, -0.5, iz.max() + 0.5]
ext_f = [ix.min() - 0.25, ix.max() + 0.25, -0.25, 16.25]
im = axs[0, 0].imshow(src2d / m, origin="lower", aspect="auto", extent=ext_s, vmin=-1, vmax=1)
axs[0, 0].set_title("source: elm-pb P on the slice nodes (ix, iz), t = 1", fontsize=10)
axs[0, 0].set_xlabel("ix"); axs[0, 0].set_ylabel("iz")
axs[0, 1].imshow(fine2d, origin="lower", aspect="auto", extent=ext_f, vmin=-1, vmax=1)
axs[0, 1].set_title("PCMS transfer, 2x resolution", fontsize=10)
axs[0, 1].set_xlabel("ix"); axs[0, 1].set_ylabel("iz")
fig.colorbar(im, ax=[axs[0, 0], axs[0, 1]], label=f"P / max|P|   (max|P| = {m:.3e})", shrink=0.8)

axs[1, 0].plot(dc, ".", ms=3, label="PCMS Source− original BOUT++ value")
axs[1, 0].set_ylim(-4e-16, 4e-16)
axs[1, 0].set_title(f"{int(coin.sum())} fine nodes that coincide with source nodes", fontsize=10)
axs[1, 0].set_xlabel("node index"); axs[1, 0].set_ylabel("difference / max|P|"); axs[1, 0].legend(loc="upper right")
dr = (rp - rb) / m
print(f"random real fractions: n={len(rb)}  PCMS identical to BOUT++ XZBilinear at {int((rp == rb).sum())};  max|PCMS - BOUT++|/max|P| = {np.abs(dr).max():.3e}")
axs[1, 1].plot(np.arange(len(rb)), dr, ".", ms=3, label="PCMS Target − BOUT++ XZBilinear")
axs[1, 1].set_ylim(-4e-16, 4e-16)
axs[1, 1].set_title(f"{len(pb)} target nodes", fontsize=10)
axs[1, 1].set_xlabel("target index"); axs[1, 1].set_ylabel("difference / max|P|"); axs[1, 1].legend(loc="upper right")
fig.savefig("transfer_verification.png", dpi=130, bbox_inches="tight")
print("wrote transfer_verification.png")
