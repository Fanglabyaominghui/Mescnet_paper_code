# Mescnet paper code

Analysis code accompanying the manuscript describing **Mescnet** (*maximum-entropy
single-cell gene network*), a network-based framework for identifying gene programs
from single-cell and spatial transcriptomic data.

The R analysis functions themselves are distributed as a separate package:
**<https://github.com/Fanglabyaominghui/Mescnet>**.

## Repository layout

```
figure1/    fig1_code.R   fig1.svg
figure2/    fig2_code.R   fig2.svg
figure3/    fig3_code.R   fig3.svg
figure4/    fig4.R        fig4.svg
figure5/    fig5_code.R   fig5_code_python.html   fig5_code2_python.html   fig5.svg
figure6/    fig6_code.R   fig6_code_python.html   fig6.svg
supplement/
  code/     sup1_code.R  sup2_code.R  sup4.code.R  sup5code.R
            sup6_code.R  sup7.code.R  sup8_code.R
  figures/  sup1.svg ... sup8.svg
  tables/   table1.csv  table2.csv  table3.csv
```

## Figure / supplement legend

| Path | Content |
| --- | --- |
| `figure1/fig1_code.R` | Mescnet conceptual framework and analysis workflow: network trajectory (UMAP + network) representation of the Waddington-like landscape and the three analysis scales of the method. |
| `figure2/fig2_code.R` | Main benchmark. Mescnet module recovery versus pathway-guided / existing network methods across PBMC and 10 *hotspot* benchmark datasets; EGAD and PPI validation of the reconstructed networks. |
| `figure3/fig3_code.R` | Colorectal cancer epithelial atlas. Builds an `INTepi` Seurat object with `Mescnet::CreaterMEscnet()`, detects modules with `ComputeMEscnetModules()`, and visualises PHATE trajectories and module activities along the proCSC -> CSC -> revCSC axis. |
| `figure4/fig4.R` | Spatial transcriptomics. SPARK-X spatially variable gene analysis of Visium HD sections, module score projection, and cross-sample / cross-species (human-mouse) correlation matrices. |
| `figure5/fig5_code.R` (+ `fig5_code_python.html`, `fig5_code2_python.html`) | Spatial adenocarcinoma-to-squamous transition (AST) program in lung cancer: RCTD deconvolution of Visium HD (R) with the reference-building / CellOracle virtual-perturbation steps in Python. |
| `figure6/fig6_code.R` (+ `fig6_code_python.html`) | Spatial microenvironmental remodelling across the LUAD-LUSC interface: Module 1 / Module 3 spatial gradients, AST score, random-forest TF prioritisation, and CellChat ligand-receptor analysis. |
| `supplement/code/sup1_code.R` | PBMC UCell benchmark of Mescnet modules versus reference gene-program methods. |
| `supplement/code/sup2_code.R` | Merged EGAD + PPI benchmark script: 10 *hotspot* datasets x 8 network methods, aggregation, ranking and per-dataset violin plots. |
| `supplement/code/sup4.code.R` | Fixed module colour library and supplementary module-level visualisations. |
| `supplement/code/sup5code.R` | Colorectal cancer Slingshot trajectory inference (proCSC -> CSC -> revCSC) with per-module UCell scoring and GO enrichment. |
| `supplement/code/sup6_code.R` | Spatial module activity maps per Mescnet module (rasterised scatter, q70-q95 thresholding). |
| `supplement/code/sup7.code.R` | TCGA-LUAD / TCGA-LUSC download and Scissor-based patient stratification. |
| `supplement/code/sup8_code.R` | CellChat network centrality and collagen-pathway communication analysis, plus TCGA LUAD M1/M3 Cox survival analysis. |
| `supplement/figures/` | Final rendered supplementary figures (SVG). `sup3.svg` has no standalone script in this release. |
| `supplement/tables/` | Supplementary tables (CSV). |

## Requirements

### R

The scripts were developed under R 4.3 / 4.4 and Seurat v5. Core dependencies:

```r
install.packages(c(
  "Seurat", "SeuratObject", "tidyverse", "dplyr", "ggplot2", "patchwork",
  "cowplot", "grid", "igraph", "ggraph", "tidygraph", "uwot", "reshape2",
  "viridis", "scales", "ggrastr", "mascarade", "ggpubr", "Matrix",
  "SingleCellExperiment", "zellkonverter", "spacexr", "SPARK", "schard",
  "UCell", "randomForest", "tibble", "slingshot", "CellChat", "Scissor",
  "ppcor", "EGAD", "AnnotationDbi", "org.Hs.eg.db", "BiocManager"
))

# Bioconductor
BiocManager::install(c("SingleCellExperiment", "zellkonverter", "AnnotationDbi", "org.Hs.eg.db"))

# GitHub
devtools::install_github("smorabit/hdWGCNA")
devtools::install_github("Fanglabyaominghui/Mescnet")
devtools::install_github("dmcable/spacexr")
```

### Python

`figure5/fig5_code_python.html`, `figure5/fig5_code2_python.html` and
`figure6/fig6_code_python.html` are exported Jupyter notebooks (scverse /
`scanpy`, `anndata`, `CellOracle`, `scipy`). Python 3.10+ was used.

## Running the code

Each script is self-contained and runs top to bottom. Two things must be checked
before running:

1. **Working directories.** The scripts keep the authors' original absolute paths
   (e.g. `E:/deskup/...` on Windows and `/home/user/Fanglab1/...` on Linux). Edit
   the `setwd()` calls and the input/output path variables at the top of each
   script to point at your local copies of the data.
2. **Input data.** Inputs are intermediate `.Rdata` / `.rds` / `.h5ad` objects
   produced by earlier steps of the pipeline (for example `plot.Rdata`,
   `INTepi.rds`, `merged_data.Rdata`, `TOMs.Rdata`, `go.Rdata`,
   `pbmc.RDS`, `benchmark.Rdata`). They are not shipped with this repository;
   see the data availability statement of the manuscript.

Several scripts (notably `fig2_code.R`, `sup2_code.R`, `sup4.code.R`) preserve the
intermediate objects they need and can be re-entered from any section, so a single
figure can be regenerated without re-running the whole benchmarking pipeline.

## Data availability

Public datasets used for benchmarking and validation are available from their
original sources: GSE178429 (PBMC), the *hotspot* benchmark collection,
TCGA-LUAD / TCGA-LUSC (GDC), and the spatial transcriptomics samples described in
the manuscript. Processed objects are available from the authors on request.

## License

The code in this repository is released under the MIT License - see `LICENSE`.

## Contact

Fang lab - <https://github.com/Fanglabyaominghui>
