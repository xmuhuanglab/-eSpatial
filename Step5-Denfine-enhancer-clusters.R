source("../eNet/R/MainFunc.R")
source("../eNet/R/utils.R")
load("/cluster/huanglab/dhong/project/eSpatial/Work/R3/E11/Signac/1.Obj_cluster_spd_E11-2.Rdata")
exp.mat = obj@assays$RNA@data
cre.mat = obj@assays$peaks@data
GPTab = GPCor(cre.mat = cre.mat,
              exp.mat = exp.mat,
              genome = "mm10")
save(GPTab,file = "1.GPTab_all.Rdata")
GPTabFilt = FindNode(GPTab = GPTab,
                     genome = "mm10")
save(GPTabFilt,file = "1.GPTabFilt.Rdata")