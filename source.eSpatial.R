# omics
library(Signac)
library(Seurat)
library(SeuratDisk)
library(SummarizedExperiment)
library(GenomeInfoDb)
# library(EnsDb.Mmusculus.v79)
# library(EnsDb.Hsapiens.v86)
library(BSgenome.Mmusculus.UCSC.mm10)
library(JASPAR2020)

#library(DropletUtils) #
library(harmony)
#library(gprofiler2)
library(loomR)
library(SeuratWrappers)
library(monocle3)
library(ArchR)
addArchRGenome("mm10")
addArchRThreads(threads = 1) 

#spatial
library(STutility)
library(parallel)
library(ica)
library(spdep)
library(jpeg)

#plots
library(ggplot2)
library(ggrepel)
library(dplyr)
library(ggsci)
library(patchwork)
library(scales)
library(RColorBrewer)
library(viridis)
library(ggpubr)
library(gplots)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)

#other
library(dplyr)
library(purrr)
library(Matrix)
library(data.table)
library(future)
library(pbapply)


source("../functions.R")
source("../as_matrix.R")