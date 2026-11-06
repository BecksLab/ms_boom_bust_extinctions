# ms_boom_bust_extinctions

Code and data for comparing how food-web reconstruction and dynamical burn-in
affect network structure and extinction robustness. The main pipeline is
implemented in Julia; downstream statistical analyses and figures are
implemented in R.

## Project layout

| Path | Contents |
|---|---|
| `code/data/pre/`, `code/data/post/` | PRE and POST community interaction data |
| `code/data/plankton/` | Plankton data used when generating communities |
| `code/data/feeding_rules.csv` | Feeding rules used by network-generation scripts |
| `code/01_buildNetworks.jl`–`code/04_getTSS.jl` | Main Julia workflow |
| `code/10_S1Bodysize.jl` | Separate, PRE-only body-size sensitivity analysis |
| `code/libs/` | Shared Julia and R helpers |
| `code/outputs/` | Intermediate and final tabular/network outputs |
| `code/*.R` | Statistical analyses and figure-generation scripts |
| `figures/`, `tables/` | Figures and tables used by the manuscript |
| `index.qmd`, `supplement/` | Quarto manuscript and supplement |

## Requirements and setup

- Julia 1.12.6 (the version recorded in `code/Manifest.toml`)
- R and the packages imported by the R scripts

Instantiate the Julia environment once:

```sh
cd code
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

The Julia project and manifest are in `code/`. R dependencies are not managed
by an `renv.lock`; check the `library()` calls in the R script you plan to run.
The R analyses use packages including `tidyverse`, `genzplyr`, `lme4`,
`emmeans`, `multcomp`, `vegan`, `FactoMineR`, `patchwork`, and `ggdendro`.

**Run the Julia and R scripts from `code/`.** Their input and output paths are
relative to the current working directory. Running a script from the repository
root, or another directory, will not resolve those paths as intended.

## Main workflow

Run the four Julia scripts in order. Network generation and extinction
simulations can both be computationally long; the later topology and TSS steps
reuse their outputs.

```sh
cd code

# 1. Generate paired PRE and POST networks
julia --project=. 01_buildNetworks.jl

# 2. Run topological and dynamic extinction simulations
JULIA_NUM_THREADS=4 julia --project=. 02_extinctionSims.jl

# 3. Export PRE network topology and per-species degrees
julia --project=. 03_getTopology.jl

# 4. Calculate species-matched TSS scores
julia --project=. 04_getTSS.jl
```

Choose `JULIA_NUM_THREADS` when starting Julia; the second script parallelizes
independent dynamic-extinction scenarios. `4` is an example, not a required
setting. Adjust it to available CPU and memory resources. Threading is not
guaranteed to make every workload faster.

### Julia steps and outputs

| Script | Reads | Writes |
|---|---|---|
| `01_buildNetworks.jl` | PRE/POST community files, plankton files, feeding rules | `outputs/pre/` and `outputs/post/`, each with `realised_network_archive.jld2`, `adjacency_matrices.jld2`, `network_summary.csv`, `species_metadata.csv`, and `downsample_metadata.csv` |
| `02_extinctionSims.jl` | PRE and POST realised-network archives | `outputs/extinction_summary.csv`, `outputs/extinction_curves.csv`, `outputs/tss_networks/pre_at_post_richness.jld2`, and `outputs/tss_networks/species.csv` |
| `03_getTopology.jl` | PRE adjacency matrices | `outputs/paleo_topology.csv` and `outputs/paleo_spp_degrees.csv` |
| `04_getTSS.jl` | PRE/POST archives and the PRE simulations saved by step 2 | `outputs/tss_results.csv`, `outputs/tss_scores.csv`, and `outputs/tss_networks/tss_scores.jld2` |

The extinction summary is a wide table with one row per network and robustness
metrics in columns. The extinction-curves table stores the per-step curves.
Step 2 saves its results after the simulations finish, so the final CSV/JLD2
outputs will not appear as completed results while that script is still
running. Step 1 likewise writes its archives after network generation.

### Species identity and TSS

Species IDs are retained alongside adjacency matrices in the same row/column
order. The archive metadata and `species_metadata.csv` provide the network
identity; `outputs/tss_networks/species.csv` additionally maps each simulated
matrix position to its `species_id`.

TSS matching uses species IDs, not matrix positions. PRE and POST networks
undergo burn-in independently, so their realised survivor lists need not be
identical. The TSS step uses the PRE generation species universe, the realised
POST species IDs, and the explicit species IDs for each simulated network.
When using or joining network matrices downstream, preserve these identifiers
and do not assume that equal matrix positions represent the same species.

## R analyses

Run the R scripts from `code/` after the Julia outputs they consume exist.
These scripts write figures under `../figures/` and, where applicable, tables
under `../tables/` or derived data under `outputs/derived/`.

```sh
cd code
Rscript 04_downsamplingEffect.R
Rscript 05_robustnessTopoVDynamic.R
Rscript 06_robustnessDownsampling.R
Rscript 07_robustnessContrast.R
Rscript 08_robustnessClustering.R
Rscript 09_robustnessComposite.R
Rscript 11_tss.R
```

The analyses use these main inputs:

| Script | Main inputs | Purpose |
|---|---|---|
| `04_downsamplingEffect.R` | `outputs/paleo_topology.csv` | Topology ordination and burn-in displacement |
| `05_robustnessTopoVDynamic.R` | `outputs/extinction_summary.csv` | Compare topological and dynamic robustness; writes a manuscript table |
| `06_robustnessDownsampling.R` | Extinction summary and PRE topology | Test how downsampling relates to robustness |
| `07_robustnessContrast.R` | Extinction summary | Contrasts among network types; writes `outputs/derived/robustnessEmmTbl.rds` |
| `08_robustnessClustering.R` | The derived RDS from script 07 | Cluster network types using estimated robustness |
| `09_robustnessComposite.R` | Extinction summary | Composite robustness figure |
| `11_tss.R` | `outputs/tss_results.csv` | Plot TSS scores |

Scripts 07 and 08 must be run in that order. The remaining scripts consume the
Julia exports and can otherwise be run independently.

## PRE-only body-size sensitivity analysis

This is a separate branch; it does not read or generate POST networks. It uses
the Türkiye PRE community and plankton inputs and compares empirical with
internally specified body sizes.

```sh
cd code
julia --project=. 10_S1Bodysize.jl
Rscript S1_bodysize.R
```

The Julia script writes its archive and tables to `outputs/S1_bodysize/`,
including `network_archive.jld2`, `species_metadata.csv`,
`network_summary.csv`, `adjacency_matrices.jld2`, `topology.csv`,
`spp_degrees.csv`, and `motif_breakdown.csv`. The R script uses those outputs
to create the body-size sensitivity figures under `figures/bodysize/`.