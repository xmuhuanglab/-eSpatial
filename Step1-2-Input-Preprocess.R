load("/cluster/huanglab/dhong/project/eSpatial/Work/R3/1.gr_embryos.Rdata")
load("/cluster/huanglab/dhong/project/Genomes/EnsDb/EnsDb.Mmusculus.v79.Rdata")

fragment = "/GSM6043255_ME11_20um_fragments.tsv.gz"
imagefile =  "/ME11_20um_spatial/"

#-----------------
#---Input file for eSpatial
#-----------------
sample = 'E11'
fragment = fragments[sample]
atac.frags = Signac::CreateFragmentObject(path = fragment)
counts = Signac::FeatureMatrix(fragments = atac.frags, features = gr)
atac.assay <- Signac::CreateChromatinAssay(
  counts = counts,
  min.features = 0,
  fragments = atac.frags
)
object = Seurat::CreateSeuratObject(counts = atac.assay, assay = "peaks")
object$sample <- sample
Annotation(object) = annotations

# add image 
image = Read10X_Image(image.dir = imagefile,filter.matrix = T)
coord = image@coordinates
rownames(coord) = paste0(rownames(coord),"-1") # RF data need
image@coordinates = coord
cells = rownames(coord)
image = image[Cells(x = object)]
DefaultAssay(object = image) <- "Spatial"
object[["slice1"]] <- image
cod = object@images$slice1@coordinates[,c("row","col")]
colnames(cod) = c("cod1","cod2")
drobj = CreateDimReducObject(
  embeddings = as.matrix(cod),
  key = "cod"
)
object@reductions[["cod"]] = drobj

# calculate gene activity
DefaultAssay(object) <- "peaks"
Annotation(object) <- annotations
gene.activities <- GeneActivity(object)
genes <- rownames(gene.activities)
genes_filter <- c(grep("Pcdh", genes), grep("PCDH", genes),
                  grep("UGT", genes), grep("Ugt", genes),
                  grep("Gm", genes), grep("Rik", genes))
gene.activities_filtered <- gene.activities[-genes_filter,]
object[['RNA']] <- CreateAssayObject(counts = gene.activities_filtered)

object = subset(object,cells = cells)
#---RNA clustering 
DefaultAssay(object) = "RNA"
object <- NormalizeData(object) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:30) %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters()
#object = FindClusters(object,resolution = 0.8) 
object$seurat_clusters_rna <- paste0("R",object$seurat_clusters)
object$seurat_clusters_rna <- factor(object$seurat_clusters_rna,
                                     levels = paste0("R",0:length(object$seurat_clusters)))

pdf(paste0("2.",sample,"_Dimplot_seurat_clusters_rna.pdf"),width = 8,height = 4)
DimPlot(object,reduction = "cod",cols = cols,group.by = "seurat_clusters_rna",pt.size = 0.5)+
  DimPlot(object,cols = cols,group.by = "seurat_clusters_rna")
dev.off()

#---ATAC clustering 
DefaultAssay(object) = "peaks"
object <- FindTopFeatures(object, min.cutoff = 10)
object <- RunTFIDF(object)
object <- RunSVD(object)
DepthCor(object,n = 30)
object <- RunUMAP(object, reduction = 'lsi', dims = 2:10) %>% # 11-2 2:10
  FindNeighbors(reduction = 'lsi', dims = 2:10) %>%
  FindClusters(algorithm = 3) 

object = FindClusters(object,resolution = 0.5) # 11-2 :0.5 11-1: 1
object$seurat_clusters_atac <- paste0("A",object$seurat_clusters)
object$seurat_clusters_atac <- factor(object$seurat_clusters_atac,
                                      levels = paste0("A",0:length(unique(object$seurat_clusters_atac))))

pdf(paste0("2.",sample,"_Dimplot_seurat_clusters_atac.pdf"),width = 8,height = 4)
DimPlot(object,reduction = "cod",cols = cols,group.by = "seurat_clusters_atac",pt.size = 0.5)+
  DimPlot(object,cols = cols,group.by = "seurat_clusters_atac")
dev.off()

# build a joint neighbor graph using both assays
object <- FindMultiModalNeighbors(
  object = object,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:10, 2:10),
  modality.weight.name = "peak.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
object <- RunUMAP(
  object = object,
  nn.name = "weighted.nn",
  assay = "peaks",
  verbose = TRUE
)

object <-  FindClusters(object,graph.name = "wsnn", algorithm = 3,resolution = 0.5)
object$seurat_clusters_joint <- paste0("J",object$seurat_clusters)
object$seurat_clusters_joint <- factor(object$seurat_clusters_joint,
                                       levels = paste0("J",0:length(object$seurat_clusters_joint)))

pdf(paste0("2.",sample,"_Dimplot_seurat_clusters_joint.pdf"),width = 8,height = 4)
DimPlot(object,reduction = "cod",cols = cols,group.by = "seurat_clusters_joint",pt.size = 0.5)+
  DimPlot(object,cols = cols,group.by = "seurat_clusters_joint")
dev.off()
saveRDS(object, file = paste0("1.Obj_RNA_cluster_",sample,".rds"))

#-----------------
#---Preprocess
#-----------------

#---RNA
DefaultAssay(object) = "RNA"
sub = DietSeurat(object,assays = "RNA")
meta = sub@meta.data
meta = meta[,c("nCount_peaks","nFeature_peaks","seurat_clusters_rna","seurat_clusters_atac","seurat_clusters_joint")]
cood = object@reductions$cod@cell.embeddings
meta = cbind(meta,cood)
sub@meta.data = meta
SaveH5Seurat(sub,filename = paste0(sample,".h5Seurat"),overwrite = T)
Convert(paste0(sample,".h5Seurat"),"h5ad",overwrite = T)

#---LSI
mat = object@reductions$lsi@cell.embeddings[,2:50]
assay = CreateAssayObject(counts = mat)
sub = CreateSeuratObject(counts = t(mat),
                         meta.data = meta)
SaveH5Seurat(sub,filename = paste0(sample,"_LSI_-1.h5Seurat"),overwrite = T)
Convert(paste0(sample,"_LSI_-1.h5Seurat"),"h5ad",overwrite = T)

#---joint
mat1 = object@reductions$lsi@cell.embeddings[,2:50]
mat2 = object@reductions$pca@cell.embeddings[,1:50]
mat = cbind(mat1,mat2)
assay = CreateAssayObject(counts = mat)
sub = CreateSeuratObject(counts = t(mat),
                         meta.data = meta)
SaveH5Seurat(sub,filename = paste0(sample,"_PCA_LSI.h5Seurat"),overwrite = T)
Convert(paste0(sample,"_PCA_LSI.h5Seurat"),"h5ad",overwrite = T)



