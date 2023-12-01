#-----------------
#---Denoised the spatial gene expression and chromatin accessibility
#-----------------

#--- prepare RNA matrix For DCA
clusters = as.character(unique(obj$spatial_domain))
i = "C11"
for (i in clusters) {
  cells = Cells(obj)[which(obj$spatial_domain == i)]
  mtx = as_matrix(obj@assays$RNA@counts[,cells])
  write.csv(mtx, paste0("./RNA1/",i,"_genes.csv"))
}


#--- prepare peak matrix For DCA
clusters = as.character(unique(obj$spatial_domain))
i = "C11"
for (i in clusters) {
  cells = Cells(obj)[which(obj$spatial_domain == i)]
  mtx = as_matrix(obj@assays$peaks@data[,cells])
  write.csv(mtx, paste0("./ATAC1/",i,"_peaks.csv"))
}

# for i in `ls *genes.csv`; do
# dca $i ${i%%.*} --threads 3 --saveweights
# done

# for i in `ls *_peaks.csv`; do
# dca $i ${i%%.*} --threads 8 --nosizefactors --nonorminput --nologinput --nocheckcounts --saveweights
# done
