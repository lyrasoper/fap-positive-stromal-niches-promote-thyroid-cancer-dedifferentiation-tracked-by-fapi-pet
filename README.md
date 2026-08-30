# FAP-positive stromal niches and thyroid cancer dedifferentiation — figure code

Bioinformatics figure code for the study of FAP-positive stromal niches, epithelial
thyroid-lineage attenuation, and [18F]FAPI-42 PET.
Guo, Ye, Pang, Xie, *et al.* (Sun Yat-sen University).

---

## Scope — read this first

This repository contains **the bioinformatics figure code only**, and it covers
**Figures 1–3 and a subset of the supplementary figures**. It is not a complete
reproduction package for the paper.

**What is here**

- Human single-cell atlas analyses and figure builders (Figure 1)
- Visium spatial-transcriptomics analyses and figure builders (Figure 2)
- Ligand–receptor nomination, deconvolution/niche analyses, the ITGA2 axis, and
  the full render of the submitted figure (Figure 3)
- The multi-omic supplementary analyses that accompany Figure 4
  (Supplementary Figs. S14 and S15)
- The ssGSEA recomputation and panel builder for Figure 5I
- The mouse source-characteristics and bulk-validation supplement
  (Supplementary Fig. S21)

**What is deliberately not here**

- Wet-lab quantification: immunoblot, IHC, multiplex immunofluorescence,
  in vivo mouse pharmacology (Figure 4 main panels, parts of Figure 2)
- The main single-cell panels of the *Braf*-driven mouse progression model
  (Figure 5A–H)
- Small-animal and clinical PET/CT analyses (Figures 6 and 7)

**What is read rather than re-run.** Several upstream analyses were run outside
this repository and their result tables are read as inputs: LIANA cell–cell
communication aggregation, RCTD deconvolution, Space Ranger and Cell Ranger
processing, and the bulk RNA/DIA-proteomic differential testing. The scripts here
run NicheNet, decoupleR ULM scoring, GSEA, the per-cell/per-spot statistics, and
all figure rendering.

---

## Figure correspondence

| Manuscript figure | Code in this repository |
|---|---|
| Figure 1 — human single-cell atlas | `scripts/fig1_human_atlas/` |
| Figure 2 — spatial transcriptomics | `scripts/fig2_spatial/` (spatial panels only) |
| Figure 3 — spatial multi-omic nomination of ITGA2 | `scripts/fig3_spatial_mechanism/` |
| Figure 4 — host *Fap* loss and α2β1–FAK–ERK signaling | main panels not included; supplementary multi-omics in `scripts/fig4_fap_multiomics/` |
| Figure 5 — *Braf*-driven mouse progression | panel I and Supplementary Fig. S21 in `scripts/fig5_mouse_braf/`; panels A–H not included |
| Figure 6 — longitudinal [18F]FAPI-42 PET in mice | not included |
| Figure 7 — paired [18F]FAPI-42 / [18F]FDG PET/CT | not included |

### A note on figure numbering inside the scripts

The manuscript was renumbered during revision. Some **output directory names and
inline comments** inside the scripts still use the earlier numbering, most often:

- `fig7_*` (for example `outputs/fig7_redesign`, `Fig7A_...`) → **current Figure 3**
- `fig8_*` (for example `fig8_supp_panel`) → **current Figure 4**
- internal `supp1`–`supp6`, `supp11_15`, `Supp Fig 19` names → this repository's own
  build order, not submission numbering

These are working paths, not scientific claims, and they were left unchanged so the
scripts keep working against existing intermediate files. Every script carries a
two-line header stating which manuscript figure it builds; trust the header.

---

## Layout

```
scripts/
├── _shared/                  paper_A house style, sourced by the figure scripts
├── fig1_human_atlas/         Figure 1 — human scRNA atlas + bulk TCGA/GEO
├── fig2_spatial/             Figure 2 — Visium spatial transcriptomics
├── fig3_spatial_mechanism/   Figure 3 — LIANA/NicheNet nomination, RCTD niche,
│                             ITGA2 axis, PAX8 distance decay, final render
├── fig4_fap_multiomics/      Supplementary Figs. S14, S15 — RNA/DIA proteomics
│                             from Fap-deficient vs wild-type hosts
└── fig5_mouse_braf/          Figure 5I ssGSEA recomputation and panel, plus
                              Supplementary Fig. S21 — mouse single-cell source
                              characteristics and bulk-cohort validation
```

Within each folder, numbered scripts run in order. Unnumbered `.py` files in
`fig3_spatial_mechanism/` are **modules, not steps** — they are imported by the
numbered entry points.

### The Figure 3 render pipeline

The submitted Figure 3 is produced by one entry point,
`13_render_fig3_final.py`, which composes five modules:

```
build_fig3_panelA_paired.py   panel A, paired sender -> receiver nomination
fig3_full_render.py           full-figure layout and panels B-H
fig3_moderate_variant.py      layout variant used by the beautify passes
fig1_fig3_revision_helpers.py loaded for its matplotlib rcParams side effects
fig3_beautify.py              the eight post-processing passes
        |
13_render_fig3_final.py       monkeypatches panel A to color couplings by
                              receiver and italicize the non-integrin receivers
                              (CD44, SDC1, DDR1), then renders
```

Run `13_render_fig3_final.py`; the earlier numbered scripts in that folder
produce the upstream tables it consumes.

---

## A rebuild that was not merged

**Figure 4C–E.** A re-derivation from the raw RNA and DIA matrices exists but
reports Hedges *g* with per-sample points, whereas the submitted legend
describes group-mean standardized module differences. It was therefore not
merged; the panels in the paper come from the earlier build.

## Data

No data are included. The scripts locate their inputs through two environment
variables (see **Configuration**). Accession numbers and access terms for the
underlying data are given in the manuscript's Data Availability statement.

Two inputs live outside the repository by design: the NicheNet prior model cache
and the bulk-mTC helper module `integrate_thyroid_markers.py`, both under
`EXTERNAL_DATA`. Scripts that need them fail with an explicit message naming the
missing path.

---

## Requirements

**R (≥ 4.2)**
`Seurat`, `Matrix`, `data.table`, `dplyr`, `tidyr`, `tibble`, `readr`, `readxl`,
`stringr`, `ggplot2`, `ggpubr`, `ggrepel`, `ggsci`, `ggtext`, `patchwork`,
`cowplot`, `scales`, `gridExtra`, `png`, `ragg`, `ComplexHeatmap`, `circlize`,
`clustree`, `Nebulosa`, `MASS`, `openxlsx`, `rstatix`, `lme4`, `lmerTest`,
`GSVA`, `decoupleR`, `nichenetr`

**Python (≥ 3.10)**
`numpy`, `pandas`, `scipy`, `matplotlib`, `seaborn`, `Pillow`, `gseapy`

---

## Configuration

The scripts read environment variables, with angle-bracket placeholders as
fall-backs:

| Variable | Meaning |
|---|---|
| `PROJECT_ROOT` | this repository's working directory (outputs are written here) |
| `EXTERNAL_DATA` | where the input objects and upstream result tables live |
| `PYTHON_BIN` | a Python with `gseapy` / `scanpy` installed (used by the GSEA step) |

For R, copy `.Renviron.example` to `.Renviron`; R reads it at startup. For Python
and shell, export the variables directly. `config.example.yaml` documents the same
keys for reference.

---

## License

See `LICENSE`.
