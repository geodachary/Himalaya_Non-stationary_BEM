# Himalaya Non-stationary Boundary Element Method

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22116240.svg)](https://doi.org/10.5281/zenodo.22116240)

Boundary-element models of interseismic locking and creep on the Main Himalayan Thrust
(MHT), using a **non-stationary** asperity formulation in which the locked–creeping
boundary is a free parameter and stress accumulated on a process zone ("ring") around
each asperity is carried explicitly.

The code discretizes the megathrust into triangular dislocation elements, compresses the
element-to-element stress-interaction matrix with an H-matrix approximation (HMMVP),
solves the quasi-static creep problem for a given set of locking boundaries, and inverts
GPS-derived observations for those boundaries with a Metropolis–Hastings MCMC sampler.

Both an **elastic half-space** setup and a **layered viscoelastic** setup (Green's
functions precomputed with VISCO3D) are included.

---

## Table of contents

- [Requirements](#requirements)
- [Quick start](#quick-start)
- [Workflow](#workflow)
- [Model configuration](#model-configuration)
- [Repository Structure](#repository-structure)
- [Input data](#input-data)
- [Results and output formats](#results-and-output-formats)
- [Large data files (Zenodo)](#large-data-files-zenodo)
- [Third-party software and data](#third-party-software-and-data)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)

---

## Requirements

| Requirement | Notes |
|---|---|
| MATLAB R2020b or newer | Developed and tested on macOS |
| Mapping Toolbox | `shaperead`, used by `tools/plot_coast_xy.m` to draw coastlines |
| Signal Processing Toolbox | `medfilt1`, used only for the slip-rate profile plots |
| C++ compiler | Only if the bundled HMMVP MEX files must be rebuilt (see `hmmvp0.16/make.m`) |
| Disk space | ~1.7 GB for the H-matrices and viscoelastic Green's functions |

Prebuilt HMMVP MEX binaries for macOS (Intel `.mexmaci64` and Apple Silicon `.mexmaca64`)
ship with the repository. On Linux or Windows, run `make.m` inside `hmmvp0.16/` first.

The Natural Earth 1:10m coastline shapefile used by the plotting helpers is bundled in
`ne_10m_coastline/`, so no external download is needed.

---

## Quick start

**Elastic model, inverting baseline length-change rates:**

```matlab
% 1. Mesh + H-matrix + Green's functions
setup_elastic                              % edit block_id / invert_vel at the top

% 2a. One forward model for a prescribed locking geometry
forward_propagation_onestep_boundary_nodes

% 2b. ...or sample the posterior
mcmc_inversion_boundary                    % edit folder_name at the top

% 3. Plot the chain and the posterior means
plot_MCMC_inversion_results
```

**Viscoelastic model:**

```matlab
setup_visco                                % edit block_id / gf_name at the top
mcmc_inversion_boundary
plot_MCMC_inversion_results
```

The setup scripts leave the mesh (`el`, `nd`), the convergence field (`rakes`, `rates`),
the Green's functions (`Ge`/`Gn`/`Gu` or `Gbase`), the observations, and the H-matrix
handle (`gamb.hmat.id`) in the base workspace. Every downstream script reads them from
there, so **run a setup script first in the same MATLAB session** — do not `clear` between
steps.

---

## Workflow

The complete workflow is four steps.

### 1. Setup — `setup_elastic.m` / `setup_visco.m`

Defines the discretized fault and locates or builds the hierarchical stress-interaction
matrix.

- `setup_elastic.m` — elastic half-space version. Builds the mesh via
  `build_mesh/make_mesh_elastic.m`, computes elastic displacement Green's functions with
  triangular dislocations (`make_dispG_triangular` / `Get_Gs_baselines`), and locates or
  builds the H-matrix through `tools/setup_Hmat.m`.
- `setup_visco.m` — viscoelastic version. Loads a precomputed mesh from `visco_mesh/` and
  the matching VISCO3D interseismic velocity Green's functions from
  `Green_function_visco/`, converts them to baseline Green's functions where needed
  (`Get_Gs_baselines_visco`), and locates the H-matrix for the same block geometry.

Both scripts then load the GPS data, optionally convert station velocities to baseline
length-change rates, plot the data, and compute the per-element self-stiffness
(`self_rate`) used to scale the ring stress.

### 2. Forward model — `forward_propagation_onestep_boundary_nodes.m`

Defines the locking boundaries, solves for interseismic stress accumulation, calculates
creep rates and the locked–creeping patch distribution, and plots the results.

For a prescribed pair of along-strike boundary curves (`locked_depths_U`,
`locked_depths_L` at `num_Dpts` nodes) and a prescribed ring stress (`ring_tau`), it:

1. rotates the mesh so the mean strike aligns with the *y* axis and interpolates the
   boundary depths onto element centroids;
2. flags locked elements (`i_locked`) and the surrounding process-zone ring (`i_ring`, of
   width `D`);
3. accumulates back-slip stress on the asperities over the time step `dT` and adds the
   stiffness-scaled ring stress;
4. solves the resulting system for creep on the unlocked elements with GMRES, using
   H-matrix matrix–vector products;
5. plots creep rate, slip-deficit rate, ring stress, predicted vs. observed data, and
   slip-rate profiles across the interface.

Set `load_inversion = true` to start from the last sample of a previous MCMC run instead
of the hard-coded boundaries.

### 3. Inversion — `mcmc_inversion_boundary.m`

Runs the same forward calculation inside a Metropolis–Hastings sampler over the locking
boundaries and ring stresses, and streams the accepted models to disk.

- **Parameters** (`3 × num_Dpts`, default `num_Dpts = 20`, so 60 unknowns):
  upper boundary depth, lower boundary depth, and ring stress at each along-strike node.
- **Proposal:** one parameter perturbed per iteration with a symmetric random step of size
  `stepsize_*`; the chain is written once per full parameter sweep (thinning).
- **Bounds:** `0 ≤ locked_depths_U ≤ locked_depths_L ≤ 40 km`, `0 ≤ ring_tau ≤ 5 MPa`.
- **Likelihood:** Gaussian, on baseline elongation rates (`invert_vel = false`) or on
  station velocities with the long-term block field added (`invert_vel = true`), weighted
  by the reported GPS uncertainties.
- **Output:** written continuously to `folder_name`; set `continuing = true` to append to
  an existing chain, `false` to overwrite it.

Watch the printed acceptance rate and retune `stepsize_*` if it drifts far from a healthy
range.

### 4. Plotting — `plot_MCMC_inversion_results.m`

The same script is used for both the elastic and the viscoelastic runs. It creates
overview plots of the obtained results: the log-likelihood trace, posterior mean
and spread of the locking boundaries, mean creep rate, probability of locking, mean ring
stress, and observed vs. predicted data. Set `folder_name` to the chain directory and
`discard` to the number of burn-in samples to drop.

`moment_rates_plot.m` is an optional add-on that converts a posterior ensemble into
moment-deficit rates and equivalent magnitudes.

---

## Model configuration

All switches live in the `BEGIN SETUP` block at the top of each script.

### Block model (convergence field and mesh)

`block_id` selects which plate-convergence field, mesh, and H-matrix are used. Rakes and
rates were derived from the Euler poles published in the corresponding studies.

| Elastic `block_id` | Viscoelastic `block_id` | Convergence field | H-matrix | Source |
|---|---|---|---|---|
| `lindsey_2018` | `himalaya_visco` | `lindsey_2018_rake_rate.txt` | `hmat/hmat_lindsey` | Lindsey et al. (2018) |
| `Block_I` | `visco_block_I` | `Block_I_rake_rate.txt` | `hmat/hmat_I` | Panda & Lindsey (2024) |
| `Block_II` | `visco_block_II` | `Block_II_rake_rate.txt` | `hmat/hmat_II` | Panda & Lindsey (2024) |
| `Block_III` | `visco_block_III` | `Block_III_rake_rate.txt` | `hmat/hmat_III` | Panda & Lindsey (2024) |

### Viscosity structure (viscoelastic runs only)

`gf_name` in `setup_visco.m` is built as `<block_id>_interseismic_vel_<viscosity>`:

| Tag | Mantle | Lower crust |
|---|---|---|
| `10eta` | 10¹⁹ Pa·s | 10¹⁹ Pa·s |
| `100eta` | 10¹⁹ Pa·s | 10²⁰ Pa·s |
| `1000eta` | 10¹⁹ Pa·s | 10²¹ Pa·s |
| `100_1000eta` | 10²⁰ Pa·s | 10²¹ Pa·s |

### Observation type

`invert_vel = false` inverts **baseline length-change rates** between station pairs (the
configuration used for the published runs); `invert_vel = true` inverts **station
velocities**, in which case the long-term block velocity field
(`Himalaya_block_velocities.txt`) is added to the model prediction.

### Stationary vs. non-stationary

The ring stress is what makes the model non-stationary. To run a **stationary** model,
zero it out: set `ring_tau = 0` in the forward script, and
`stepsize_ring_tau = 0e6*ones(num_Dpts,1)` together with a zero starting value in the MCMC
script, so the ring parameters never move.

### Other physical parameters

| Variable | Default | Meaning |
|---|---|---|
| `num_Dpts` | 20 | along-strike boundary nodes |
| `D` | 3 km | process-zone (ring) width |
| `dT` | 10 yr (MCMC), 100 yr (forward) | stress-accumulation time step |
| `mu` | 3 × 10¹⁰ Pa | shear modulus |
| `intv` | 12 km | nominal triangular patch side length |
| `origin` | [83°E, 28°N] | origin of the local Cartesian frame |
| `Niterations` | 10⁶ | MCMC iterations |

---

## Repository Structure

A typical layout (names as used in the scripts):

```
.
├── setup_elastic.m                              # step 1 — elastic setup
├── setup_visco.m                                # step 1 — viscoelastic setup
├── forward_propagation_onestep_boundary_nodes.m # step 2 — forward model
├── mcmc_inversion_boundary.m                    # step 3 — MCMC inversion
├── plot_MCMC_inversion_results.m                # step 4 — plots (elastic + viscoelastic)
├── moment_rates_plot.m                          # optional moment-rate analysis
│
├── build_mesh/            # fault geometry, convergence field, mesh generation
├── visco_mesh/            # precomputed meshes for the viscoelastic runs
├── Green_function_visco/  # VISCO3D interseismic velocity Green's functions
├── hmat/                  # H-matrix stress kernels (not in git — see Zenodo)
├── vel_field_lindsey/     # GPS velocities, baselines, long-term block velocities
├── ne_10m_coastline/      # Natural Earth coastline shapefile, used by the plots
├── MCMC_results/          # posterior-mean results for the published runs
│
├── tools/                 # geometry, Green's-function, solver, plotting, and
│                          #   input-data helpers
└── hmmvp0.16/             # HMMVP hierarchical-matrix software (third party)
```

### `build_mesh/`

```
him_contours_45.txt          make_mesh_elastic.m
lindsey_2018_rake_rate.txt   make_mesh_visco.m
Block_I_rake_rate.txt
Block_II_rake_rate.txt
Block_III_rake_rate.txt
```

**`him_contours_45.txt`** — depth contours of the MHT (`lon, lat, depth in km, negative
downward`), the geometry the mesh is built from. It was constructed from the Slab2 model
of Hayes et al. (2018) combined with the structural model of Hubbard et al. (2016):

- In the **Nepal section**, Slab2 was modified so the interface reproduces the ramp–flat
  structure of the MHT reported by Hubbard et al. (2016).
- In the **remaining sections**, Slab2 and Hubbard et al. (2016) were combined.
- Hubbard et al. (2016) report the geometry only to **30 km** depth, so depths were
  **linearly extrapolated to 45 km** using the dip measured between 25 and 30 km. The same
  procedure was extended to **75 km** for the viscoelastic setup.

**`*_rake_rate.txt`** — plate convergence at the element centroids
(`lon, lat, rake in degrees, rate in mm/yr`), one file per block model. Rake is the motion
direction of the overriding plate relative to the subducting plate, right-hand rule (thumb
along strike, index finger down dip, palm down). These fields were obtained by taking the
**Euler poles published in Lindsey et al. (2018) and Panda & Lindsey (2024) and converting
them to rake and rate** at each patch (see `tools/get_rake_rate_from_EulerPole.m`).

**`make_mesh_elastic.m` / `make_mesh_visco.m`** — reconstruct the triangular mesh for the
elastic and viscoelastic setups, interpolate rake and rate onto patch centroids, and plot
the mesh, the convergence field, and the strike/dip vectors. The mesh algorithm follows
code by Yo Fukushima (2010). The current configuration produces **9,474 triangular
elements**.

### `visco_mesh/`

Precomputed meshes used by the viscoelastic runs, one per block geometry, matching the
VISCO3D Green's functions element for element:

```
himalaya_visco_mesh_inversion.mat    visco_block_II_mesh_inversion.mat
visco_block_I_mesh_inversion.mat     visco_block_III_mesh_inversion.mat
```

### `Green_function_visco/`

Interseismic surface velocity Green's functions computed with the spectral-element code
VISCO3D, named `<block_id>_interseismic_vel_<viscosity>.mat` — 4 block models × 4
viscosity structures = 16 files, ~46 MB each. Each file provides `Ge_inter`, `Gn_inter`,
and `Gu_inter` (east, north, up velocity per unit slip rate on each element).

### `hmat/`

H-matrix approximations of the element-to-element stress-interaction kernel, one directory
per block geometry:

```
hmat_lindsey/    hmat_I/    hmat_II/    hmat_III/
```

Each holds a single `Hmat_ne9474hmat1rerr-3.00.dat` (~229 MB), built at an element-wise
relative accuracy of 10⁻³. The naming encodes the element count, the H-matrix flag, and
log₁₀ of the tolerance, so it is regenerated automatically whenever the mesh changes.
**This directory is excluded from git** — see [Large data files](#large-data-files-zenodo).

### `tools/` and `hmmvp0.16/`

- `tools/` — geometry, Green's-function, solver, plotting, and input-data helpers.
- `hmmvp0.16/` — HMMVP hierarchical-matrix software (third party, see below).

---

## Input data

`vel_field_lindsey/` holds the observations. The top-level copies of
`observed_vel_subset.txt` and `Himalaya_block_velocities.txt` are the ones the scripts
load by default.

| File | Contents |
|---|---|
| `Observed_velocities.txt` | Full GPS velocity field of Lindsey et al. (2018): `lon, lat, Ve, Vn, σe, σn, corr, …, site` (mm/yr) |
| `observed_vel_subset.txt` | **The subset used in this study** — 215 sites, `lon, lat, Ve, Vn, σe, σn`. The Shillong block was not modelled here (only the Himalayan arc), so those sites are excluded |
| `Himalaya_block_velocities.txt` | Long-term block velocities at the same sites, in the same row order: `lon, lat, Ve, Vn` (mm/yr). Used only when `invert_vel = true` |
| `vbase.txt` | Observed baseline length-change rates between GPS station pairs (583 baselines, mm/yr) |

Uncertainties below 1 mm/yr are doubled in the setup scripts as an error floor. Baselines
that cross elements breaking the surface are handled by integrating strain rates across
them (`tools/make_baseline_rate_changes.m`).

---

## Results and output formats

### Raw MCMC chains

`mcmc_inversion_boundary.m` writes one row per retained sample to `folder_name`:

| File | Columns |
|---|---|
| `M_locked_depths_U.txt` | upper locking boundary depth at each of the `num_Dpts` nodes (km) |
| `M_locked_depths_L.txt` | lower locking boundary depth at each node (km) |
| `M_ring_tau.txt` | ring stress at each node (Pa) |
| `locked_index.txt` | locked (1) / creeping (0) flag per element |
| `creep_rates.txt` | creep rate per element (mm/yr) |
| `dhat.txt` | predicted data (baseline rates or velocities, mm/yr) |
| `logrho.txt` | log-likelihood of the sample |

### Posterior-mean results

`MCMC_results/` contains the posterior means of the published inversion runs, organized
as:

```
MCMC_results/
├── Elastic_mean_results/     Bline_lindsey_2018, Bline_Block_I/II/III
├── Visco_mean_results/       <block_id>_<viscosity>  (16 combinations)
├── Stationary_mean_results/  stationary (zero ring stress) reference run
└── Ensemble_mean_results/    mean across the model ensemble
```

Each directory holds plain text files with these columns:

| File | Columns | Units |
|---|---|---|
| `mean_creep_rates.txt` | `lon, lat, depth, mean creep rate` — one row per element centroid | depth km (negative down), rate mm/yr |
| `mean_locked_index.txt` | `lon, lat, depth, mean probability of locking` | probability 0–1 |
| `mean_ring_tau.txt` | `lon, lat, depth, mean ring stress` | — |
| `mean_stressing_rates.txt` | `lon, lat, depth, mean stressing rate` (ensemble folder) | — |
| `mean_dhat.txt` | `lon1, lat1, lon2, lat2, mean length-change rate` — one row per baseline, giving the coordinates of the two GPS stations and the rate of change of the distance between them | mm/yr |

Element files have 9,474 rows (one per triangular patch); baseline files have 583 rows.

---

## Large data files (Zenodo)

The `hmat/` directory contains four H-matrix stress-interaction kernels used by the elastic and viscoelastic boundary-element calculations. Each file is approximately 229 MB (~1 GB total). Because GitHub rejects regular Git objects larger than 100 MB and is not well suited for distributing files of this size, **`hmat/` is excluded from the repository via `.gitignore` and archived separately on Zenodo.**

The archived dataset is:

> Acharya, D. (2026). *H-matrix Stress Kernels for Himalaya Non-stationary Boundary Element Method* [Dataset]. Zenodo.
> https://doi.org/10.5281/zenodo.22114918

Download the archived H-matrix dataset, extract it, and place the `hmat/` directory at the repository root before running the models. Preserve the following directory structure:

```text
hmat/
├── hmat_lindsey/
│   └── Hmat_ne9474hmat1rerr-3.00.dat
├── hmat_I/
│   └── Hmat_ne9474hmat1rerr-3.00.dat
├── hmat_II/
│   └── Hmat_ne9474hmat1rerr-3.00.dat
└── hmat_III/
    └── Hmat_ne9474hmat1rerr-3.00.dat
```

### SHA-256 checksums

The integrity of the downloaded H-matrix files can be verified using the following SHA-256 checksums:

| File                                              | SHA-256                                                            |
| ------------------------------------------------- | ------------------------------------------------------------------ |
| `hmat/hmat_lindsey/Hmat_ne9474hmat1rerr-3.00.dat` | `b4592f23e128faa007a1512ebaef556b547de920ebe14dab167c0428a2e60641` |
| `hmat/hmat_I/Hmat_ne9474hmat1rerr-3.00.dat`       | `05da266e3e06adae9879374d67518922ea74e2aa59de248d9e3e36823b2ba5a0` |
| `hmat/hmat_II/Hmat_ne9474hmat1rerr-3.00.dat`      | `2be73ad932c3f52f3ba89836f8a324591aa90c53e821eff23b22dde7c976e7ac` |
| `hmat/hmat_III/Hmat_ne9474hmat1rerr-3.00.dat`     | `ebd898d8560e707fb7f76becbca13baef978c9fdb1e8d407048b7641344833bc` |

The compressed archive deposited on Zenodo can additionally be verified using:

```text
Himalaya_Non-stationary_BEM_hmat_v1.0.0.tar.gz
SHA-256: a6e548b46e301e50ee8de7f9ec7ef3dac4c7cb74e959a86d6c92c8dbbdeae151
```

Alternatively, the required H-matrix can be generated locally. Edit `setup_elastic.m` or `setup_visco.m` and set:

```matlab
compute_hmat = true;
```

If the H-matrix corresponding to the selected `block_id` is missing, the setup routine will generate it automatically. **This calculation is substantially more computationally expensive than the model calculations that follow**, so using the archived Zenodo files is recommended whenever possible.

After the required H-matrix has been generated, set:

```matlab
compute_hmat = false;
```

again for subsequent model runs.


---

## Third-party software and data

This repository bundles third-party components that are **not** covered by its own
license:

- **HMMVP 0.16** by Andrew M. Bradley (`hmmvp0.16/`), distributed under the **Eclipse
  Public License 1.0**; see `hmmvp0.16/readme.txt`.
- **Triangular dislocation elements** by Brendan Meade (`tools/CalcTriDisps.m`,
  `tools/tridisloc3d.m`), Copyright © 2006 Brendan Meade, MIT license (notice in the file
  headers) — Meade (2007), after Comninou & Dunders (1975). The variants
  `tools/CalcTriDisps_O.m` and `tools/CalcTriStrains_O.m` are modifications by
  Wen-Jeng Huang. **These are the routines used by the model.**
- **Triangular dislocation routines** by Mehdi Nikkhoo (`tools/Nikkhoo/`,
  `tools/TDstressHS.m`), Copyright © 2014 Mehdi Nikkhoo, MIT license (notice in the file
  headers) — Nikkhoo & Walter (2015). Bundled as an **optional alternative kernel only**;
  every call site is commented out in `tools/CalcG_rake.m`, so no published result depends
  on it.
- **Mesh generation** based on code by Yo Fukushima (2010), adapted in
  `build_mesh/make_mesh_*.m` and `tools/make_tri_mesh.m`.
- **Natural Earth** 1:10m Coastline shapefile, version 4.1.0 (`ne_10m_coastline/`) —
  public domain, <https://www.naturalearthdata.com>.
- **Observational data** redistributed from Lindsey et al. (2018) — see
  [Input data](#input-data).
- **Fault geometry** derived from Slab2 (Hayes et al., 2018) and Hubbard et al. (2016).
- **Viscoelastic Green's functions** (`Green_function_visco/`) computed with VISCO3D
  (Pollitz, 2014, 2025); the input files follow the method of Sherrill et al. (2026).
  See [Citation](#citation) for the full references.

Before public release, confirm and document the license, source, and version of every
bundled third-party component.

---

## Citation

This section records how the model was developed, including the software implementation,
large numerical datasets, methods, codes, and observational datasets on which it builds.

### Software

If you use this software, please cite:

> Acharya, D., Johnson, K. M., & Sherrill, E. M. (2026). *Himalaya Non-stationary BEM*
> (Version v1.0.0) [Computer software]. Zenodo.
> https://doi.org/10.5281/zenodo.22116240

### Associated publication

> *(Accompanying paper — citation to be added upon publication.)*

### H-matrix dataset

The H-matrix stress kernels required for the model are archived separately as:

> Acharya, D. (2026). *H-matrix Stress Kernels for Himalaya Non-stationary Boundary Element
> Method* [Dataset]. Zenodo. https://doi.org/10.5281/zenodo.22114918

### HMMVP

Its author asks that work using the H-matrix code cite:

> A. M. Bradley (2012). *HMMVP: Software to Compute Matrix-Vector Products with an
> H-Matrix.*

### Methods and data sources

The references below are the sources of the methods, codes, and data the model is built
from; each is credited in the section it belongs to.

The **non-stationary asperity approach** and the model setup are adapted from:

> Johnson, K., & Sherrill, E. (2026, March 9). *Eroding asperities imply larger locked
> regions on subduction megathrusts.* PREPRINT (Version 1), Research Square.
> https://doi.org/10.21203/rs.3.rs-8714153/v1

The **viscoelastic Green's functions** were computed with the spectral-element method
(VISCO3D) following:

> Pollitz, F. F. (2014). Post-earthquake relaxation using a spectral element method: 2.5-D
> case. *Geophysical Journal International*, 198(1), 308–326.
> https://doi.org/10.1093/gji/ggu114

> Pollitz, F. F. (2025). 3D Viscoelastic Models of Slip-Deficit Rate Along the Cascadia
> Subduction Zone. *Journal of Geophysical Research: Solid Earth*, 130(1), e2024JB029847.
> https://doi.org/10.1029/2024JB029847

The **VISCO3D input files** were created by adopting the method used in:

> Sherrill, E. M., Johnson, K. M., & Pollitz, F. F. (2026). Competing effects of elastic
> heterogeneity and viscous flow on interseismic coupling at Cascadia. *Geophysical
> Research Letters*, 53, e2026GL124108. https://doi.org/10.1029/2026GL124108

The **fault geometry** distributed here (`build_mesh/him_contours_45.txt`) is a new
composite surface built for this study. It combines the two models below — modifying Slab2
to reproduce the ramp–flat structure of the MHT reported by Hubbard et al. (2016) in the
Nepal section, merging the two models elsewhere, and extrapolating the interface from the
30 km limit of Hubbard et al. (2016) down to 45 km (and to 75 km for the viscoelastic
setup) using the dip measured between 25 and 30 km. See
[`build_mesh/`](#build_mesh) for the full construction.

If you use this fault geometry, please cite the software release:

> Acharya, D., Johnson, K. M., & Sherrill, E. M. (2026). *Himalaya Non-stationary BEM*
> (Version v1.0.0) [Computer software]. Zenodo.
> https://doi.org/10.5281/zenodo.22116240

together with the source models it is derived from:

> Hayes, G. P., Moore, G. L., Portner, D. E., Hearne, M., Flamme, H., Furtney, M., &
> Smoczyk, G. M. (2018). Slab2, a comprehensive subduction zone geometry model. *Science*,
> 362(6410), 58–61. https://doi.org/10.1126/science.aat4723

> Hubbard, J., Almeida, R., Foster, A., Sapkota, S. N., Bürgi, P., & Tapponnier, P. (2016).
> Structural segmentation controlled the 2015 Mw 7.8 Gorkha earthquake rupture in Nepal.
> *Geology*, 44(8), 639–642. https://doi.org/10.1130/G38077.1

The **GPS data, block models, and Euler poles** come from:

> Lindsey, E. O., Almeida, R., Mallick, R., Hubbard, J., Bradley, K., Tsang, L. L. H., …
> Hill, E. M. (2018). Structural Control on Downdip Locking Extent of the Himalayan
> Megathrust. *Journal of Geophysical Research: Solid Earth*, 123(6), 5265–5278.
> https://doi.org/10.1029/2018JB015868

> Panda, D., & Lindsey, E. O. (2024). Overriding Plate Deformation Controls Inferences of
> Interseismic Coupling Along the Himalayan Megathrust. *Journal of Geophysical Research:
> Solid Earth*, 129(9), e2024JB029819. https://doi.org/10.1029/2024JB029819

The **triangular dislocation elements** used for the stress kernel and the elastic Green's
functions are those of:

> Meade, B. J. (2007). Algorithms for the calculation of exact displacements, strains, and
> stresses for triangular dislocation elements in a uniform elastic half space. *Computers
> & Geosciences*, 33(8), 1064–1075. https://doi.org/10.1016/j.cageo.2006.12.003

---

## License

The **BSD 3-Clause License** in `LICENSE` applies to the original Himalaya Non-stationary
BEM code authored by Durga Acharya. It **does not** replace the terms governing bundled
third-party software or data — see
[Third-party software and data](#third-party-software-and-data).

---

## Contact

Durga Acharya
