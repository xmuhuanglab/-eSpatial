#-----------------
#---Define gene expression modules
#-----------------
genes = as.character(unique(spatialLinksFilt$Gene))
exp.mat = as.matrix(exp.avg)
cre.mat = as.matrix(peaks.avg)
#---gene kmeans 
genes = intersect(rownames(exp.mat),spatialLinksFilt$Gene)
exp.mat = exp.mat[genes,]
exp.Z = t(scale(t(exp.mat),center = T))
k = 7
exp.km = kmeans(exp.Z, k,nstart = 25)
i = 1
cluster = exp.km$cluster
kmeansInfo = data.frame(row.names = names(cluster),
                        peaks = names(cluster),
                        kmeans = paste0("k",cluster))
kDomain = lapply(X = seq_along(1:k),FUN = function(i){
  tmp = cluster[cluster == i]
  genes = names(tmp)
  means = colMeans(exp.mat[genes,])
  df = matrix(data = means,nrow = 1)
}) %>% do.call(what = rbind)

rownames(kDomain) = paste0("k",1:k)
colnames(kDomain) = colnames(exp.mat)

#---order enhancer modules
modules = list("M1" = c("C0","C1","C2"),
               "M2" = c("C5","C6"),
               "M3" = "C4",
               "M4" = 'C3')
clusters = unlist(modules)

kDomain = kDomain[,clusters]

pheatmap(kDomain,
         scale = "row",
         name = "Enrichment",
         cluster_rows = T,
         cluster_cols = F,
         show_rownames = T,
         show_colnames = T,
         color = ArchRPalettes$blueYellow,
         breaks = c(seq(-2,2,length=300)),
         legend = T,
         legend_labels = NA,
         annotation_legend = T,
         #annotation_row = anno_row,
         # annotation_col = anno_cols,
         # annotation_colors = list(module = (module_cols)),
         fontsize = 18)

kmeansorder = c("k1","k4","k5","k7","k2","k3","k6")
#kmeansorder = c("k1","k3","k2","k6","k7","k9","k5","k4","k8") #hrg
kDF = lapply(X = seq_along(1:k),FUN = function(i){
  tmp = kmeansInfo %>% filter(kmeans == kmeansorder[i])
}) %>% do.call(what = rbind)


b = data.frame(colMeans(exp.mat))
exp.adj = exp.mat
#exp.adj[,c("C12","C9")] = exp.adj[,c("C12","C9")]*0.8
# exp.adj[,"C7"] = exp.adj[,"C7"]*0.5
# exp.adj[,"C4"] = exp.adj[,"C4"]*0.6
exphm = exp.adj[kDF$peaks,clusters]
ArchRPalettes$solarExtra
anno_row = data.frame(row.names = rownames(kDF),module = kDF$kmeans)
p = pheatmap(exphm,
             scale = "row",
             name = "Enrichment",
             cluster_rows = F,
             cluster_cols = F,
             show_rownames = F,
             show_colnames = T,
             color = ArchRPalettes$blueYellow,
             #breaks = c(seq(-2,2,length=300)),
             legend = T,
             legend_labels = NA,
             annotation_legend = T,
             annotation_row = anno_row,
             # annotation_col = anno_cols,
             # annotation_colors = list(module = (module_cols)),
             fontsize = 18)

kDomain = t(kDomain)
kDomain = kDomain[,kmeansorder]
kDomain = kDomain[clusters,]

#-----------------
#---Define spatial patterns of enhancers
#-----------------
modules = list("M1" = c("C0","C1","C2"),
               "M2" = c("C5","C6"),
               "M3" = "C4",
               "M4" = 'C3')
clusters = unlist(modules)
cre.mat = cre.mat[,clusters]


kmeansorder = c("k1","k4","k5","k7","k2","k3","k6")
peakod = data.frame()
geneod = data.frame()

m = 1

cutoff = mean(cre.mat) + sd(cre.mat)
cre.hm = matrix()
for (m in 1:length(kmeansorder)) {
  module = kmeansorder[m]
  domain = kDomain[,module]
  domain = names(domain)[which(domain > 0.03)]
  rmdomain = rownames(kDomain)[-which(rownames(kDomain) %in% domain)]
  kgene = kDF %>% filter(kmeans == module)
  kgene = kgene$peaks
  kgene = intersect(kgene,rownames(exp.mat))
  
  #---peak module
  kpeak = spatialLinksFilt %>% filter(Gene %in% kgene)
  kpeak = kpeak$Enhancer
  kpeak = intersect(rownames(cre.mat),kpeak)
  cre.sub = cre.adj[kpeak,domain]
  
  cre.bin = cre.sub
  cre.bin[cre.bin >= cutoff] = 1
  cre.bin[cre.bin < cutoff] = 0
  if(length(domain) > 1){
    #---gene module
    comb = lapply(X = kpeak, FUN = function(x){
      com = paste0(cre.bin[x,],collapse = "")
      com = data.frame(peaks = x,expDomainNum = length(domain),geneModule = paste0("m",m),
                       type = com, enhDomainNum = sum(cre.bin[x,]))
    }) %>% do.call(what = rbind)
    
    #---filtered ramdom combination
    ncom = round(nrow(comb)*0.005)
    comInfo = data.frame(table(comb$type))
    radom = comInfo$Var1[which(comInfo$Freq < ncom)]
    comb$type[which(comb$type %in% radom)] = "0-radom"
    combType = unique(comb$type)
    combType = combType[order(combType,decreasing = T)]
    
    combType = combType[-length(combType)]
    
    pod = lapply(X = seq_along(combType), FUN = function(x){
      pp = comb %>% filter(type == combType[x])
      pp$type = paste0("M",x)
      return(pp)
    })%>% do.call(what = rbind)
    
    pod.sub = pod %>% filter(enhDomainNum > 0)
    cre.bin.hm = cre.bin[pod.sub$peaks,domain]
    #adjust
    # pp = pod.sub$peaks[pod.sub$type == "M1"]
    # cre.bin.hm[pp,"C2"] = 1
    anno_row = data.frame(row.names = pod.sub$peaks,type = pod.sub$type)
    annotation_cols = cols[1:length(unique(pod.sub$type))]
    names(annotation_cols) = unique(pod.sub$type)
    
    p = pheatmap(cre.bin.hm,
                 #border_color = "NA",
                 #scale = "row",
                 name = "Enrichment",
                 cluster_rows = F,
                 cluster_cols = F,
                 show_rownames = F,
                 show_colnames = T,
                 color = c("white","#D51F26"),
                 #breaks = c(seq(-3,3,length=300)),
                 legend = T,
                 legend_labels = NA,
                 annotation_legend = T,
                 annotation_row = anno_row,
                 # annotation_col = anno_cols,
                 annotation_colors = list(type = annotation_cols),
                 fontsize = 18)
    pdf(paste0("./",module,"_EMhm.pdf"),height = 8,width = 6)
    print(p)
    dev.off()
  }
  
}




