# Epilepsy

## Overview

This repository provides the R code used in the study:

"Uncovering the key circuit FOSL2/FOS/EGR3/EGR1 contributing to hyperexcitability of excitatory neurons in epileptic temporal cortex and hippocampus"

Epilepsy is mainly characterized by spontaneous seizures caused by hyperactive neural circuits. To delineate the celltype-specific mechanisms underlying neuronal hyperexcitability, we resolve the hyperexcitability of excitatory neurons across epileptic human brain transfoci at single-cell resolution to identify key drivers and diagnostic signatures.

------------------------------------------------------------------------

## Main Results

-   Dissecting transregional cellular composition independent and joint-triggering epileptic effects
-   Excitatory neurons exhibit highly active state in separate epileptic regions
-   Excitatory neurons exert hyperactive influence in epileptic human temporal cortex and hippocampus
-   Glial cells mediate cellular junction assembly and synaptic organization functions along the hippocampal anterior and posterior axis
-   Revealing circuit FOSL2/FOS/EGR3/EGR1 transregional crosstalk promoting excitatory neuronal activation


------------------------------------------------------------------------

## Code Organization

The scripts are organized by results to facilitate reproducibility. Each folder contains the specific analysis code used to generate the corresponding results presented in the manuscript.

For example:
- `code/Result1/` contains scripts for the complete analysis and visualizations of the Result "Dissecting transregional cellular composition independent and joint-triggering epileptic effects".
- `code/Result2/` contains scripts for the complete analysis and visualizations of the Result "Excitatory neurons exhibit highly active state in separate epileptic regions".
- Subsequent folders correspond to the remaining results.

Users can follow these scripts to reproduce the results step by step.

------------------------------------------------------------------------

## Repository Structure

Epilepsy/
├── README.md
├── LICENSE
├── code/
│   ├── Result1/
│   │   ├── Data Processing_Hippocampus.R
│   │   ├── Data Processing_Temporal cortex.R
│   │   ├── DEGs_Figure1G_Hippocampus.R
│   │   ├── DEGs_Figure1G_Temporal cortex.R
│   │   ├── Figure1G_Temp_Hippo.R
│   │   ├── scCODA.ipynb
│   │   └── sccoda_data.R
│   │
│   ├── Result2/
│   │   ├── Cell subpopulation re-annotation_Hoppocampus.R
│   │   ├── Cell subpopulation re-annotation_Temporal cortex.R
│   │   ├── Hotspot_addmodulescore.R
│   │   ├── Hotspot_Hippocampus.ipynb
│   │   ├── Hotspot_module_enrichment.R
│   │   ├── hotspot_Temporal cortex.ipynb
│   │   ├── Neuronal_activation_gene_expression.R
│   │   ├── Neuronal_activation_Score.R
│   │   └── Transcriptomic similarity.R
│   │
│   ├── Result3/
│   │   ├── monocle3_Hippocampus.R
│   │   └── monocle3_Temporal cortex.R
│   │
│   ├── Result4/
│   │   ├── Glial_subtype.R
│   │   └── Glial_enrichment.R
│   │
│   ├── Result5/
│   │   ├── pySCENIC_Data preparation_Visualization.R
│   │   ├── scenic.bash.R
│   │   ├── Neuronal_activation.R
│   │   └── Bulk RNA data processing.R
│   │
│
└── sessionInfo.txt

------------------------------------------------------------------------

## Data Availability

The datasets used in this study are publicly available:

- The data that support the findings of this study are openly available in GEO at GSE160189 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160189) and GitHub (https://github.com/khodosevichlab/Epilepsy19). 
- Bulk RNA-seq datasets analyzed in this study were sourced from the Gene Expression Omnibus (GEO) database, including GSE256068 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE256068), GSE139914 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139914), and GSE140393 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140393).

Due to file size limitations, raw datasets are not included in this repository.

------------------------------------------------------------------------

## Requirements

R = 4.2.2

Packages: Seurat, AUCell, CBNplot, ggsankey, GOplot, ClusterGVis, clusterProfiler, jjAnno

------------------------------------------------------------------------

## Contact

Yunpeng Zhang: zhangyp@hrbmu.edu.cn\
Xia Li: lixia@hrbmu.edu.cn

------------------------------------------------------------------------

## License

MIT License
