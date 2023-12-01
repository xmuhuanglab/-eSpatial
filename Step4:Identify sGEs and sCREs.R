
object <- readRDS("../1.object_RNA_cluster_E11.rds")
#---add cell neigbors
coord = object@reductions$cod@cell.embeddings
nCells = nearestCells(coord = coord,
                      n=6)

#---add stagate joint
meta = read.csv("../STAGATE/E11_PCA_LSI_meta_pl_a0.8r0.7.csv")
rownames(meta) = meta$X
object$stagate_joint = paste0("C",meta[,"louvain"])
object$stagate_joint <- factor(object$stagate_joint,
                               levels = paste0("C",0:length(unique(object$stagate_joint))))
#---smooth
sm = smoothBydist(nCells = nCells,
                  metadata = object@meta.data,
                  smoothBy = "stagate_joint")
object$stagate_joint_sm = sm$round1
object$stagate_joint_sm <- factor(object$stagate_joint_sm,
                                  levels = paste0("C",0:length(unique(object$stagate_joint_sm))))
#---add stagate umap
umap = read.csv("../STAGATE/E11_PCA_LSI_umap_pl_a0.8r0.7.csv")
umap = data.frame(row.names = Cells(object),umap[,2:3])
colnames(umap) = c("umap1","umap2")
drobject = CreateDimReducobjectect(
  embeddings = as.matrix(umap),
  key = "umap"
)
object@reductions[["atac_umap"]] = drobject

pdf("2.E11_Dimplot_stagate_joint.pdf",width = 12,height = 4)
DimPlot(object,reduction = "cod",cols = cols,group.by = "stagate_joint",pt.size = 0.5)+
  DimPlot(object,reduction = "cod",cols = cols,group.by = "stagate_joint_sm",pt.size = 0.5)+
  DimPlot(object,reduction = "atac_umap",cols = cols,group.by = "stagate_joint_sm",pt.size = 0.5)
dev.off()

#---add stagate ATAC 
meta = read.csv("../STAGATE/E11_LSI_meta_pl_a0.8r0.7.csv")
rownames(meta) = meta$X
object$stagate_atac = paste0("A",meta[,"louvain"])
object$stagate_atac <- factor(object$stagate_atac,
                              levels = paste0("A",0:length(unique(object$stagate_atac))))
#---smooth
sm = smoothBydist(nCells = nCells,
                  metadata = object@meta.data,
                  smoothBy = "stagate_atac")
object$stagate_atac_sm = sm$round1
object$stagate_atac_sm <- factor(object$stagate_atac_sm,
                                 levels = paste0("A",0:length(unique(object$stagate_atac_sm))))
#---add stagate umap
umap = read.csv("../STAGATE/E11_LSI_umap_pl_a0.8r0.7.csv")
umap = data.frame(row.names = Cells(object),umap[,2:3])
colnames(umap) = c("umap1","umap2")
drobject = CreateDimReducobjectect(
  embeddings = as.matrix(umap),
  key = "umap"
)
object@reductions[["atac_umap"]] = drobject

pdf("2.E11_Dimplot_stagate_atac.pdf",width = 12,height = 4)
DimPlot(object,reduction = "cod",cols = cols,group.by = "stagate_atac",pt.size = 0.5)+
  DimPlot(object,reduction = "cod",cols = cols,group.by = "stagate_atac_sm",pt.size = 0.5)+
  DimPlot(object,reduction = "atac_umap",cols = cols,group.by = "stagate_atac_sm",pt.size = 0.5)
dev.off()

#---add stagate RNA 
meta = read.csv("../STAGATE/E11_rna_meta_pl_a0.8r0.7.csv")
rownames(meta) = meta$X
object$stagate_rna = paste0("R",meta[,"louvain"])
object$stagate_rna <- factor(object$stagate_rna,
                             levels = paste0("R",0:length(unique(object$stagate_rna))))
#---smooth
sm = smoothBydist(nCells = nCells,
                  metadata = object@meta.data,
                  smoothBy = "stagate_rna")
object$stagate_rna_sm = sm$round1
object$stagate_rna_sm <- factor(object$stagate_rna_sm,
                                levels = paste0("R",0:length(unique(object$stagate_rna_sm))))
#---add stagate umap
umap = read.csv("../STAGATE/E11_rna_umap_pl_a0.8r0.7.csv")
umap = data.frame(row.names = Cells(object),umap[,2:3])
colnames(umap) = c("umap1","umap2")
drobject = CreateDimReducobjectect(
  embeddings = as.matrix(umap),
  key = "umap"
)
object@reductions[["atac_umap"]] = drobject

pdf("2.E11_Dimplot_stagate_rna.pdf",width = 12,height = 4)
DimPlot(object,reduction = "cod",cols = cols,group.by = "stagate_rna",pt.size = 0.5)+
  DimPlot(object,reduction = "cod",cols = cols,group.by = "stagate_rna_sm",pt.size = 0.5)+
  DimPlot(object,reduction = "atac_umap",cols = cols,group.by = "stagate_rna_sm",pt.size = 0.5)
dev.off()


#---rename idents by class 
Idents(object) = object$stagate_joint_sm
object = RenameIdents(object, 
                      #cns
                      "C3" = "C0",
                      "C4" = "C1",
                      #spine
                      "C6" = "C2",
                      "C0" = 'C3',
                      
                      "C2" = 'C4',
                      "C1" = "C5",
                      "C5" = "C6"
)

object$spatial_domain = Idents(objectect = object)

pdf("2.E11_Dimplot_spd.pdf",width = 12,height = 4)
DimPlot(object,reduction = "cod",cols = cols,group.by = "stagate_atac",pt.size = 0.5)+
  DimPlot(object,reduction = "cod",cols = cols,group.by = "spatial_domain",pt.size = 0.5)+
  DimPlot(object,reduction = "atac_umap",cols = cols,group.by = "spatial_domain",pt.size = 0.5)
dev.off()

pdf("2.E11_Dimplot_spd_split.pdf",width = 9,height = 9)
DimPlot(object,reduction = "cod",cols = cols,split.by = "spatial_domain",
        pt.size = 0.5,ncol = 3,group.by = "spatial_domain")
dev.off()
save(object,file = "1.object_cluster_spd_E11.Rdata")
#---find markers
Idents(object) = object$spatial_domain
DefaultAssay(object) = "RNA"
markers <- FindAllMarkers(object, only.pos = T, logfc.threshold = 0.1, min.pct = 0.05) 
markers = markers %>% filter(p_val_adj < 0.05)
markers = markers %>% group_by(gene) %>% filter(avg_log2FC==max(avg_log2FC))
save(markers, file = "1.markers_spd_E11.Rdata")

markers = markers %>% filter(avg_log2FC > 0.2)
genes = unique(markers$gene)
save(markers, file = "1.markers_spd_E11_2306.Rdata")

#---Foldchange
clusters = as.character(unique(object$spatial_domain))
#---genes fc
DefaultAssay(object) = "RNA"
genes = rownames(object)
exp.avg = data.frame(matrix(nrow = length(genes),ncol = length(clusters)),row.names = genes)
colnames(exp.avg) = clusters
for (i in clusters) {
  avgfc = FoldChange(object, ident.1 = i)
  exp.avg[,i] = avgfc$avg_log2FC 
}
save(exp.avg,file ="1.E11_avgFC_spd_genes.Rdata")

#---peaks fc
DefaultAssay(object) = "peaks"
peaks = rownames(object)
peaks.avg = data.frame(matrix(nrow = length(peaks),ncol = length(clusters)),row.names = peaks)
colnames(peaks.avg) = clusters
for (i in clusters) {
  avgfc = FoldChange(object, ident.1 = i)
  peaks.avg[,i] = avgfc$avg_log2FC 
}
save(peaks.avg,file ="1.E11_avgFC_spd_peaks.Rdata")
