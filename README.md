# OxideFlow — results page

A single-file, dependency-free results page for **OxideFlow**: an MLIP-only pipeline that turns a
bulk oxide crystal into relaxed surfaces with vacancies and adsorbates, plus a 15-model benchmark
of those energies against five published DFT studies.

Everything is in `index.html` — the charts are hand-rolled inline SVG built from the run data at
load time, so there is no build step, no bundler and no runtime dependency. The only external
request is the Google Fonts stylesheet (Archivo / IBM Plex Sans / IBM Plex Mono); the page falls
back to system fonts cleanly if it is blocked.

## Deploy to GitHub Pages

```bash
# from inside this folder
git init -b main
git add .
git commit -m "OxideFlow results page"

# create the public repo and push in one step (the URL below assumes "oxideflow-page")
gh repo create oxideflow-page --public --source=. --remote=origin --push

# turn on Pages from main, root folder
gh api -X POST repos/{owner}/oxideflow-page/pages \
  -f 'source[branch]=main' -f 'source[path]=/'
```

Live a minute or two later at `https://<your-username>.github.io/oxideflow-page/`.

If you prefer the UI: **Settings → Pages → Source: Deploy from a branch → main / (root)**.

To serve it from `https://<your-username>.github.io/` directly, name the repo
`<your-username>.github.io` instead.

## Local preview

```bash
python3 -m http.server 8000   # then open http://localhost:8000
```

## Editing the numbers

All figure data lives in one JSON object assigned to `const D = {...}` in the first `<script>`
block near the bottom of `index.html`:

| key | drives |
|---|---|
| `scatter` | trend-vs-accuracy scatter and the sortable 15-checkpoint table (per-model O*/OH* MAE and Spearman ρ, params, OVFE MAE, s/relax) |
| `elements` | periodic convergence map (per-oxide converged counts out of 30, for O* and OH*) |
| `walltime` | cost chart (s/relax, s/step, steps, n) |
| `ovfe` | oxygen-vacancy chart (per-model E_vac at O2c and O3c, MAE) |
| `gas` | gas-reference spread chart (per-species range across models, mean abs. deviation from median) |

Values in the prose, the "best model per benchmark" table, the UMA pair table and the config-tuning
table are written inline in the HTML rather than driven by `D`.

## Header animation

The animation is an OxideFlow relaxation, not a video: H2O on rutile IrO2(110), 195 atoms, 600
steps, of which steps 0–300 (every step) are embedded in the page as the `const T = {...}` JSON
and drawn as a ball model on a `<canvas>` — no WebGL, no library. Only atoms above z = 20 Å (the top layer) are drawn; periodic images are chosen per atom relative
to frame 0 so nothing jumps across the cell boundary; the slab is drawn as muted ball-and-stick
and the adsorbate (the last `nads` atoms in the blob) as full-size balls. `--model` / `--walltime` (or `"model"` / `"walltime"` keys in the blob) add the run line to the caption.
The embedded run is SevenNet-omni-omat24 on IrO2(110), Ir5c site, 1.0 Å nudge, 600 steps, 222 s. The right-hand panel reads `T.steps`, `T.E`, `T.fmax` (largest force on a free atom) and `T.OH`
(the O–H distance of the bond that breaks). To rebuild it from a new trajectory:

```bash
python hero/make_trajectory.py path/to/traj.xyz --stop 300 --stride 1 --zcut 20 --bond O H --model "SevenNet-omni-omat24" --walltime "3 min 42 s" > hero/traj_blob.json
# paste hero/traj_blob.json over the `const T = {...};` line in index.html
```

Element radii and colours live in the `RAD` / `COL` tables at the top of the hero script
(Ir is stored as "I"; other elements by first letter). Camera: `yaw` and `pitch` on the same lines.

**To use a rendered video instead** (OVITO, Blender): put it at `hero/opening.mp4` and replace the
`<canvas id="atoms" …>` element with
`<video autoplay muted playsinline loop poster="hero/poster.png"><source src="hero/opening.mp4" type="video/mp4"></video>`.
Delete the right-hand panel if the video carries its own annotation.

## Files

```
index.html                       the whole page
hero/traj_blob.json              the embedded animation data (steps 0–300 of the IrO2 / H2O run)
hero/make_trajectory.py          regenerate it from any ASE-readable trajectory
.nojekyll                        stops GitHub Pages running the content through Jekyll
README.md                        this file
```

## Credit

Sean Song · NSERC USRA S2026 · UBC Materials Engineering — Sustainable Materials Lab + UBC
Catalysis Lab. Supervisor: Karthik Akkiraju. All studies ran on UBC ARC Sockeye
(allocation `st-akkiraju-1`). Pipeline source: https://github.com/dehydrogenated/oxideflow
