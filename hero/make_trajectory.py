"""Build the header-animation data from an ASE-readable trajectory.

    python make_trajectory.py traj.xyz --stop 300 --stride 5 --zcut 18 --bond O H > traj_blob.json

then paste the output over the `const T = {...};` line in index.html.

--stride N   keep every Nth frame          --start/--stop   frame range (inclusive stop)
--zcut Z     drop atoms below z=Z (Å) so only the top layers are drawn (0 = keep all)
--bond A B   track the distance between the A atom and the B atom that end up furthest apart
             among the adsorbate atoms (used for the "bond being broken" readout); omit to skip
The adsorbate (first atoms of species B, e.g. H) is centred in the cell with periodic wrapping.
Species are stored by first letter, except Ir -> "I"; add a RAD/COL entry in index.html for new ones.
"""
import argparse, json, sys
import numpy as np
from ase.io import read

ap = argparse.ArgumentParser()
ap.add_argument("path"); ap.add_argument("--start", type=int, default=0); ap.add_argument("--stop", type=int, default=None)
ap.add_argument("--stride", type=int, default=5); ap.add_argument("--zcut", type=float, default=0.0)
ap.add_argument("--bond", nargs=2, default=None)
ap.add_argument("--model", default=None); ap.add_argument("--walltime", default=None)
a = ap.parse_args()

frames = read(a.path, index=":")
stop = len(frames) - 1 if a.stop is None else min(a.stop, len(frames) - 1)
steps = list(range(a.start, stop + 1, a.stride))
f0 = frames[0]; cell = f0.cell; A, B = cell[0, 0], cell[1, 1]
syms = f0.get_chemical_symbols(); z0 = f0.positions[:, 2]
o = hl = None; ads = []
if a.bond:
    sa, sb = a.bond
    bs = [i for i, s in enumerate(syms) if s == sb]
    o = min([i for i, s in enumerate(syms) if s == sa], key=lambda i: np.linalg.norm(f0.positions[i] - f0.positions[bs[0]]))
    pe = frames[stop].positions; hl = max(bs, key=lambda i: np.linalg.norm(pe[o] - pe[i]))
    shift = np.array([A/2 - f0.positions[o, 0], B/2 - f0.positions[o, 1], 0.0])
    ads = [o] + bs                       # drawn as full balls; everything else as the lattice
else:
    shift = np.zeros(3)
slab = [i for i in range(len(f0)) if z0[i] > a.zcut and i not in ads]
# drop slab atoms that have no bond partner among the drawn slab atoms (they would float unattached)
P = f0.positions
slab = [i for i in slab if any(syms[i] != syms[j] and np.linalg.norm(P[i] - P[j]) < 2.35 for j in slab if j != i)]
keep = slab + ads   # adsorbate last
try:
    free = ~np.array([c.index.tolist() for c in f0.constraints][0] if f0.constraints else [], dtype=int).any() if False else None
except Exception: free = None
from ase.constraints import FixAtoms
fixed = set(); [fixed.update(c.index.tolist()) for c in f0.constraints if isinstance(c, FixAtoms)]
freemask = np.array([i not in fixed for i in range(len(f0))])
zmid = (z0[keep].max() + z0[keep].min()) / 2
ref = f0.positions[keep] + shift; ref[:, 0] %= A; ref[:, 1] %= B     # wrap frame 0 once ...

E, D, FM, out = [], [], [], []
for s in steps:
    fr = frames[s]; p = fr.positions
    try: E.append(round(float(fr.get_potential_energy()), 3))
    except Exception: E.append(0.0)
    try: FM.append(round(float(np.linalg.norm(fr.get_forces(), axis=1)[freemask].max()), 3))
    except Exception: FM.append(0.0)
    if o is not None: D.append(round(float(np.linalg.norm(p[o] - p[hl])), 3))
    q = p[keep] + shift                                              # ... then keep each atom in the image nearest its frame-0 position
    q[:, 0] -= np.round((q[:, 0] - ref[:, 0]) / A) * A; q[:, 1] -= np.round((q[:, 1] - ref[:, 1]) / B) * B
    q[:, 0] -= A/2; q[:, 1] -= B/2; q[:, 2] -= zmid
    out.append([round(float(v), 2) for v in q.ravel()])

blob = {"n": len(keep), "nads": len(ads), "fmax": FM, "sym": "".join("I" if syms[i] == "Ir" else syms[i][0] for i in keep),
        "cell": [round(float(A), 3), round(float(B), 3)], "steps": steps, "E": E, "frames": out,
        "total_atoms": len(f0), "source": f"{a.path}, steps {a.start}-{stop} every {a.stride}"}
if D: blob["OH"] = D
if a.model: blob["model"] = a.model
if a.walltime: blob["walltime"] = a.walltime
json.dump(blob, sys.stdout, separators=(",", ":"))
