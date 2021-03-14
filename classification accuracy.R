#classification accuracy
randomcell<-read.csv("sp.random-celllist50.csv",head=T,row.names = 1)
lastcell<-read.csv("sp.random-lastcelllist50.csv",head=T,row.names = 1)
geneuse<-read.csv("sp.random-genelist50.csv",head=T,row.names = 1)
genes.use<-geneuse$scmap

sp.integrated=RunPCA(sp.integrated, npcs=30, verbose=FALSE, features=genes.use)
sp.integrated=RunUMAP(sp.integrated, reduction="pca", dims=1:8)
sp.integrated=FindNeighbors(sp.integrated, dims=1:10)
sp.integrated=FindClusters(object=sp.integrated, resolution=0.6, verbose=FALSE)

ref<-as.data.frame(sp.integrated@assays$RNA@counts[rownames(sp.integrated@assays$RNA@counts)%in%geneuse$scmap,colnames(sp.integrated@assays$RNA@counts)%in%randomcell$x])
ref<-t(ref)
ref<-as.data.frame(ref)
ref$cluster<-sp.integrated@meta.data$seurat_clusters[rownames(sp.integrated@meta.data)%in%randomcell$x]
colnames(ref)<-gsub("-","_",colnames(ref))
ref[,2001]<-factor(ref[,2001])
test<-as.data.frame(sp.integrated@assays$RNA@counts[rownames(sp.integrated@assays$RNA@counts)%in%geneuse$scmap,colnames(sp.integrated@assays$RNA@counts)%in%lastcell$x])
test<-t(test)
test<-as.data.frame(test)
test$cluster<-sp.integrated@meta.data$seurat_clusters[rownames(sp.integrated@meta.data)%in%lastcell$x]
colnames(test)<-gsub("-","_",colnames(test))
test[,2001]<-factor(test[,2001])


set.seed(3)
data.ref<-randomForest(cluster~.,data=ref,mtry=2000,importance = T,proximity=T)
print(data.ref)
wp<-predict(data.ref,test)
res<-data.frame(name=character(175),predict=numeric(175),cluster=numeric(175))
res$name<-rownames(test)
res$predict<-wp
res$cluster<-test$cluster
res2<-table(test[,2001],wp)




#names
write.csv(res,"scmap-predict50.csv")
write.csv(res2,"scmap-cp-predict50.csv")



miss<-sum(res$predict!=res$cluster)
E<-list()
E[1]<-miss/nrow(res)

genes.use<-geneuse$scran
sp.integrated=RunPCA(sp.integrated, npcs=30, verbose=FALSE, features=genes.use)
sp.integrated=RunUMAP(sp.integrated, reduction="pca", dims=1:8)
sp.integrated=FindNeighbors(sp.integrated, dims=1:10)
sp.integrated=FindClusters(object=sp.integrated, resolution=0.6, verbose=FALSE)

ref<-as.data.frame(sp.integrated@assays$RNA@counts[rownames(sp.integrated@assays$RNA@counts)%in%geneuse$scran,colnames(sp.integrated@assays$RNA@counts)%in%randomcell$x])
ref<-t(ref)
ref<-as.data.frame(ref)
ref$cluster<-sp.integrated@meta.data$seurat_clusters[rownames(sp.integrated@meta.data)%in%randomcell$x]
colnames(ref)<-gsub("-","_",colnames(ref))
ref[,2001]<-factor(ref[,2001])
test<-as.data.frame(sp.integrated@assays$RNA@counts[rownames(sp.integrated@assays$RNA@counts)%in%geneuse$scran,colnames(sp.integrated@assays$RNA@counts)%in%lastcell$x])
test<-t(test)
test<-as.data.frame(test)
test$cluster<-sp.integrated@meta.data$seurat_clusters[rownames(sp.integrated@meta.data)%in%lastcell$x]
colnames(test)<-gsub("-","_",colnames(test))
test[,2001]<-factor(test[,2001])


set.seed(3)
data.ref<-randomForest(cluster~.,data=ref,mtry=2000,importance = T,proximity=T)
print(data.ref)
wp<-predict(data.ref,test)
res<-data.frame(name=character(175),predict=numeric(175),cluster=numeric(175))
res$name<-rownames(test)
res$predict<-wp
res$cluster<-test$cluster
res2<-table(test[,2001],wp)

#names
write.csv(res,"scran-predict50.csv")
write.csv(res2,"scran-cp-predict50.csv")

miss<-sum(res$predict!=res$cluster)

E[2]<-miss/nrow(res)

genes.use<-geneuse$hvg
sp.integrated=RunPCA(sp.integrated, npcs=30, verbose=FALSE, features=genes.use)
sp.integrated=RunUMAP(sp.integrated, reduction="pca", dims=1:8)
sp.integrated=FindNeighbors(sp.integrated, dims=1:10)
sp.integrated=FindClusters(object=sp.integrated, resolution=0.6, verbose=FALSE)

ref<-as.data.frame(sp.integrated@assays$RNA@counts[rownames(sp.integrated@assays$RNA@counts)%in%geneuse$hvg,colnames(sp.integrated@assays$RNA@counts)%in%randomcell$x])
ref<-t(ref)
ref<-as.data.frame(ref)
ref$cluster<-sp.integrated@meta.data$seurat_clusters[rownames(sp.integrated@meta.data)%in%randomcell$x]
colnames(ref)<-gsub("-","_",colnames(ref))
ref[,2001]<-factor(ref[,2001])
test<-as.data.frame(sp.integrated@assays$RNA@counts[rownames(sp.integrated@assays$RNA@counts)%in%geneuse$hvg,colnames(sp.integrated@assays$RNA@counts)%in%lastcell$x])
test<-t(test)
test<-as.data.frame(test)
test$cluster<-sp.integrated@meta.data$seurat_clusters[rownames(sp.integrated@meta.data)%in%lastcell$x]
colnames(test)<-gsub("-","_",colnames(test))
test[,2001]<-factor(test[,2001])

set.seed(3)
data.ref<-randomForest(cluster~.,data=ref,mtry=2000,importance = T,proximity=T)
print(data.ref)
wp<-predict(data.ref,test)
res<-data.frame(name=character(175),predict=numeric(175),cluster=numeric(175))
res$name<-rownames(test)
res$predict<-wp
res$cluster<-test$cluster
res2<-table(test[,2001],wp)

#names
write.csv(res,"hvg-predict50.csv")
write.csv(res2,"hvg-cp-predict50.csv")

miss<-sum(res$predict!=res$cluster)

E[3]<-miss/nrow(res)

genes.use<-geneuse$sct
sp.integrated=RunPCA(sp.integrated, npcs=30, verbose=FALSE, features=genes.use)
sp.integrated=RunUMAP(sp.integrated, reduction="pca", dims=1:8)
sp.integrated=FindNeighbors(sp.integrated, dims=1:10)
sp.integrated=FindClusters(object=sp.integrated, resolution=0.6, verbose=FALSE)

ref<-as.data.frame(sp.integrated@assays$RNA@data[rownames(sp.integrated@assays$RNA@data)%in%geneuse$sct,colnames(sp.integrated@assays$RNA@data)%in%randomcell$x])
ref<-t(ref)
ref<-as.data.frame(ref)
ref$cluster<-sp.integrated@meta.data$seurat_clusters[rownames(sp.integrated@meta.data)%in%randomcell$x]
colnames(ref)<-gsub("-","_",colnames(ref))
ref[,2001]<-factor(ref[,2001])
test<-as.data.frame(sp.integrated@assays$RNA@data[rownames(sp.integrated@assays$RNA@data)%in%geneuse$sct,colnames(sp.integrated@assays$RNA@data)%in%lastcell$x])
test<-t(test)
test<-as.data.frame(test)
test$cluster<-sp.integrated@meta.data$seurat_clusters[rownames(sp.integrated@meta.data)%in%lastcell$x]
colnames(test)<-gsub("-","_",colnames(test))
test[,2001]<-factor(test[,2001])


set.seed(3)
data.ref<-randomForest(cluster~.,data=ref,mtry=2000,importance = T,proximity=T)
print(data.ref)
wp<-predict(data.ref,test)
res<-data.frame(name=character(175),predict=numeric(175),cluster=numeric(175))
res$name<-rownames(test)
res$predict<-wp
res$cluster<-test$cluster
res2<-table(test[,2001],wp)

#names
write.csv(res,"sct-predict50.csv")
write.csv(res2,"sct-cp-predict50.csv")

miss<-sum(res$predict!=res$cluster)

E[4]<-miss/nrow(res)

genes.use<-geneuse$rogue_norm
sp.integrated=RunPCA(sp.integrated, npcs=30, verbose=FALSE, features=genes.use)
sp.integrated=RunUMAP(sp.integrated, reduction="pca", dims=1:8)
sp.integrated=FindNeighbors(sp.integrated, dims=1:10)
sp.integrated=FindClusters(object=sp.integrated, resolution=0.6, verbose=FALSE)

ref<-as.data.frame(sp.integrated@assays$RNA@data[rownames(sp.integrated@assays$RNA@data)%in%geneuse$rogue_norm,colnames(sp.integrated@assays$RNA@data)%in%randomcell$x])
ref<-t(ref)
ref<-as.data.frame(ref)
ref$cluster<-sp.integrated@meta.data$seurat_clusters[rownames(sp.integrated@meta.data)%in%randomcell$x]
colnames(ref)<-gsub("-","_",colnames(ref))
ref[,2001]<-factor(ref[,2001])
test<-as.data.frame(sp.integrated@assays$RNA@data[rownames(sp.integrated@assays$RNA@data)%in%geneuse$rogue_norm,colnames(sp.integrated@assays$RNA@data)%in%lastcell$x])
test<-t(test)
test<-as.data.frame(test)
test$cluster<-sp.integrated@meta.data$seurat_clusters[rownames(sp.integrated@meta.data)%in%lastcell$x]
colnames(test)<-gsub("-","_",colnames(test))
test[,2001]<-factor(test[,2001])


set.seed(3)
data.ref<-randomForest(cluster~.,data=ref,mtry=2000,importance = T,proximity=T)
print(data.ref)
wp<-predict(data.ref,test)
res<-data.frame(name=character(175),predict=numeric(175),cluster=numeric(175))
res$name<-rownames(test)
res$predict<-wp
res$cluster<-test$cluster
res2<-table(test[,2001],wp)

#names
write.csv(res,"rogue-norm-predict50.csv")
write.csv(res2,"rogue-norm-cp-predict50.csv")

miss<-sum(res$predict!=res$cluster)

E[5]<-miss/nrow(res)

genes.use<-geneuse$rogue
sp.integrated=RunPCA(sp.integrated, npcs=30, verbose=FALSE, features=genes.use)
sp.integrated=RunUMAP(sp.integrated, reduction="pca", dims=1:8)
sp.integrated=FindNeighbors(sp.integrated, dims=1:10)
sp.integrated=FindClusters(object=sp.integrated, resolution=0.6, verbose=FALSE)

ref<-as.data.frame(sp.integrated@assays$RNA@counts[rownames(sp.integrated@assays$RNA@counts)%in%geneuse$rogue,colnames(sp.integrated@assays$RNA@counts)%in%randomcell$x])
ref<-t(ref)
ref<-as.data.frame(ref)
ref$cluster<-sp.integrated@meta.data$seurat_clusters[rownames(sp.integrated@meta.data)%in%randomcell$x]
colnames(ref)<-gsub("-","_",colnames(ref))
ref[,2001]<-factor(ref[,2001])
test<-as.data.frame(sp.integrated@assays$RNA@counts[rownames(sp.integrated@assays$RNA@counts)%in%geneuse$rogue,colnames(sp.integrated@assays$RNA@counts)%in%lastcell$x])
test<-t(test)
test<-as.data.frame(test)
test$cluster<-sp.integrated@meta.data$seurat_clusters[rownames(sp.integrated@meta.data)%in%lastcell$x]
colnames(test)<-gsub("-","_",colnames(test))
test[,2001]<-factor(test[,2001])


set.seed(3)
data.ref<-randomForest(cluster~.,data=ref,mtry=2000,importance = T,proximity=T)
print(data.ref)
wp<-predict(data.ref,test)
res<-data.frame(name=character(175),predict=numeric(175),cluster=numeric(175))
res$name<-rownames(test)
res$predict<-wp
res$cluster<-test$cluster
res2<-table(test[,2001],wp)

#names
write.csv(res,"rogue-predict50.csv")
write.csv(res2,"rogue-cp-predict50.csv")

miss<-sum(res$predict!=res$cluster)

E[6]<-miss/nrow(res)

genes.use<-geneuse$vst
sp.integrated=RunPCA(sp.integrated, npcs=30, verbose=FALSE, features=genes.use)
sp.integrated=RunUMAP(sp.integrated, reduction="pca", dims=1:8)
sp.integrated=FindNeighbors(sp.integrated, dims=1:10)
sp.integrated=FindClusters(object=sp.integrated, resolution=0.6, verbose=FALSE)

ref<-as.data.frame(sp.integrated@assays$RNA@data[rownames(sp.integrated@assays$RNA@data)%in%geneuse$vst,colnames(sp.integrated@assays$RNA@data)%in%randomcell$x])
ref<-t(ref)
ref<-as.data.frame(ref)
ref$cluster<-sp.integrated@meta.data$seurat_clusters[rownames(sp.integrated@meta.data)%in%randomcell$x]
colnames(ref)<-gsub("-","_",colnames(ref))
ref[,2001]<-factor(ref[,2001])
test<-as.data.frame(sp.integrated@assays$RNA@data[rownames(sp.integrated@assays$RNA@data)%in%geneuse$vst,colnames(sp.integrated@assays$RNA@data)%in%lastcell$x])
test<-t(test)
test<-as.data.frame(test)
test$cluster<-sp.integrated@meta.data$seurat_clusters[rownames(sp.integrated@meta.data)%in%lastcell$x]
colnames(test)<-gsub("-","_",colnames(test))
test[,2001]<-factor(test[,2001])


set.seed(3)
data.ref<-randomForest(cluster~.,data=ref,mtry=2000,importance = T,proximity=T)
print(data.ref)
wp<-predict(data.ref,test)
res<-data.frame(name=character(175),predict=numeric(175),cluster=numeric(175))
res$name<-rownames(test)
res$predict<-wp
res$cluster<-test$cluster
res2<-table(test[,2001],wp)

#names
write.csv(res,"vst-predict50.csv")
write.csv(res2,"vst-cp-predict50.csv")

miss<-sum(res$predict!=res$cluster)

E[7]<-miss/nrow(res)

genes.use<-geneuse$disp
sp.integrated=RunPCA(sp.integrated, npcs=30, verbose=FALSE, features=genes.use)
sp.integrated=RunUMAP(sp.integrated, reduction="pca", dims=1:8)
sp.integrated=FindNeighbors(sp.integrated, dims=1:10)
sp.integrated=FindClusters(object=sp.integrated, resolution=0.6, verbose=FALSE)

ref<-as.data.frame(sp.integrated@assays$RNA@data[rownames(sp.integrated@assays$RNA@data)%in%geneuse$disp,colnames(sp.integrated@assays$RNA@data)%in%randomcell$x])
ref<-t(ref)
ref<-as.data.frame(ref)
ref$cluster<-sp.integrated@meta.data$seurat_clusters[rownames(sp.integrated@meta.data)%in%randomcell$x]
colnames(ref)<-gsub("-","_",colnames(ref))
ref[,2001]<-factor(ref[,2001])
test<-as.data.frame(sp.integrated@assays$RNA@data[rownames(sp.integrated@assays$RNA@data)%in%geneuse$disp,colnames(sp.integrated@assays$RNA@data)%in%lastcell$x])
test<-t(test)
test<-as.data.frame(test)
test$cluster<-sp.integrated@meta.data$seurat_clusters[rownames(sp.integrated@meta.data)%in%lastcell$x]
colnames(test)<-gsub("-","_",colnames(test))
test[,2001]<-factor(test[,2001])


set.seed(3)
data.ref<-randomForest(cluster~.,data=ref,mtry=2000,importance = T,proximity=T)
print(data.ref)
wp<-predict(data.ref,test)
res<-data.frame(name=character(175),predict=numeric(175),cluster=numeric(175))
res$name<-rownames(test)
res$predict<-wp
res$cluster<-test$cluster
res2<-table(test[,2001],wp)

#names
write.csv(res,"disp-predict50.csv")
write.csv(res2,"disp-cp-predict50.csv")

miss<-sum(res$predict!=res$cluster)

E[8]<-miss/nrow(res)

genes.use<-geneuse$schs
sp.integrated=RunPCA(sp.integrated, npcs=30, verbose=FALSE, features=genes.use)
sp.integrated=RunUMAP(sp.integrated, reduction="pca", dims=1:8)
sp.integrated=FindNeighbors(sp.integrated, dims=1:10)
sp.integrated=FindClusters(object=sp.integrated, resolution=0.6, verbose=FALSE)

ref<-as.data.frame(sp.integrated@assays$RNA@data[rownames(sp.integrated@assays$RNA@data)%in%geneuse$schs,colnames(sp.integrated@assays$RNA@data)%in%randomcell$x])
ref<-t(ref)
ref<-as.data.frame(ref)
ref$cluster<-sp.integrated@meta.data$seurat_clusters[rownames(sp.integrated@meta.data)%in%randomcell$x]
colnames(ref)<-gsub("-","_",colnames(ref))
ref[,2001]<-factor(ref[,2001])
test<-as.data.frame(sp.integrated@assays$RNA@data[rownames(sp.integrated@assays$RNA@data)%in%geneuse$schs,colnames(sp.integrated@assays$RNA@data)%in%lastcell$x])
test<-t(test)
test<-as.data.frame(test)
test$cluster<-sp.integrated@meta.data$seurat_clusters[rownames(sp.integrated@meta.data)%in%lastcell$x]
colnames(test)<-gsub("-","_",colnames(test))
test[,2001]<-factor(test[,2001])


set.seed(3)
data.ref<-randomForest(cluster~.,data=ref,mtry=2000,importance = T,proximity=T)
print(data.ref)
wp<-predict(data.ref,test)
res<-data.frame(name=character(175),predict=numeric(175),cluster=numeric(175))
res$name<-rownames(test)
res$predict<-wp
res$cluster<-test$cluster
res2<-table(test[,2001],wp)

#names
write.csv(res,"schs-predict50.csv")
write.csv(res2,"schs-cp-predict50.csv")

miss<-sum(res$predict!=res$cluster)

E[9]<-miss/nrow(res)

error.rates<-unlist(E)
write.csv(error.ra*tes,"error.rates.predict50.csv")

setwd("C:/Users/zyn05/Desktop/gene-select/spBM1 predict/")
library(ggplot2)
library(reshape2)
error<-read.csv("spBM1-predict-2000-errorrates.csv",head=T,row.names=1)
error<-error[,c(1,3,5,7,9,11,13,15,17)]
colnames(error)<-c("scmap","scran","hvg","sct","rogue_norm","rogue","vst","disp","schs")
error<-melt(error)
colnames(error)<-c("method","rates")
error$rates<-1-error$rates
error$method<-factor(error$method,levels=c("hvg","vst","sct","disp","rogue","rogue_norm","scmap","scran","schs"))
ggplot(error,aes(method,rates))+geom_boxplot()
ggplot(error,aes(method,rates))+geom_violin()+theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
ggplot(error,aes(method,rates)) +    geom_violin(aes(fill=method)) +      geom_jitter(shape=".",position = position_jitter(width = 0.4))+     theme_bw()+theme(legend.position = "none")+theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))

