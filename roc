setwd("C:/Users/zyn05/Desktop/gene-select")
##
library(Seurat)
library(dplyr)
library(Matrix)
library(reticulate)
library(ggplot2)
library(ROCR)



##1. ???荼??????蠊菇?seurat????
abc.data=Read10X(data.dir="C:/Users/zyn05/Desktop/gene-select/")
abc=CreateSeuratObject(counts=abc.data, project="abc", min.cells=50, min.features=1, names.field=1, names.delim="_")
#????indiv,tissue??lin??息
W=read.table("01_meta.txt", sep="\t", header=T, row.names=1)
abc@meta.data$indiv=W$indiv
abc@meta.data$tissue=W$tissue
abc@meta.data$lin=W$lin

##2. ????一系?械?QC
abc=subset(abc, subset=nFeature_RNA >= 1000 & nFeature_RNA <= 8000)
abc=subset(abc, subset=nCount_RNA >= 2000 & nCount_RNA <= 500000)

##3. SubsetData----indiv.list[[indiv]]??spBM1,spBM2...
X = rownames(abc@meta.data)
abc.hspc=subset(abc,cells=grep(pattern="spBM1", X, value=T))

##4. NormalizeData + FindVariableFeatures
abc.hspc=NormalizeData(abc.hspc, verbose=FALSE)
sp.integrated=abc.hspc
sp.integrated=ScaleData(sp.integrated, verbose=FALSE)

setwd("C:/Users/zyn05/Desktop/gene-select/spBM1 predict/")


filename <- sapply(c(1:50),function(x){
  paste("sp.random-celllist",x, ".csv", sep='')})
path <- getwd()
filePath <- sapply(filename, function(x){
  paste(path,x,sep='/')})
r.cell <- lapply(filePath, function(x){
  read.csv(x,row.names = 1)})

filename2 <- sapply(c(1:50),function(x){
  paste("sp.random-lastcelllist",x, ".csv", sep='')})
path <- getwd()
filePath2 <- sapply(filename2, function(x){
  paste(path,x,sep='/')})
l.cell <- lapply(filePath2, function(x){
  read.csv(x,row.names = 1)})

filename3 <- sapply(c(1:50),function(x){
  paste("sp.random-genelist",x, ".csv", sep='')})
path <- getwd()
filePath3 <- sapply(filename3, function(x){
  paste(path,x,sep='/')})
gene <- lapply(filePath3, function(x){
  read.csv(x,row.names = 1)})

outPath <- getwd()##????路??
out_fileName <- sapply(c(1:50),function(x){
  paste("vst-prob-",x, ".csv", sep='')}) ##csv??式
out_filePath  <- sapply(out_fileName, function(x){
  paste(outPath ,x,sep='/')}) ##????路????

out_fileName2 <- sapply(c(1:50),function(x){
  paste("vst-auc-",x, ".csv", sep='')}) ##csv??式
out_filePath2  <- sapply(out_fileName2, function(x){
  paste(outPath ,x,sep='/')}) ##????路????

E<-list()

for (k in 10:50){

genes.use<-gene[[k]][,7]
abc.hspc=NormalizeData(abc.hspc, verbose=FALSE)
sp.integrated=abc.hspc
sp.integrated=ScaleData(sp.integrated, verbose=FALSE)
sp.integrated=RunPCA(sp.integrated, npcs=30, verbose=FALSE, features=genes.use)
sp.integrated=RunUMAP(sp.integrated, reduction="pca", dims=1:8)
sp.integrated=FindNeighbors(sp.integrated, dims=1:10)
sp.integrated=FindClusters(object=sp.integrated, resolution=0.6, verbose=FALSE)

 randomcell<-r.cell[[k]]$x
  lastcell<-l.cell[[k]]$x
  ref<-as.data.frame(sp.integrated@assays$RNA@counts[rownames(sp.integrated@assays$RNA@counts)%in%genes.use,colnames(sp.integrated@assays$RNA@counts)%in%randomcell])
  ref<-t(ref)
  ref<-as.data.frame(ref)
  ref$cluster<-sp.integrated@meta.data$seurat_clusters[rownames(sp.integrated@meta.data)%in%randomcell]
  colnames(ref)<-gsub("-","_",colnames(ref))
  ref$cluster<-factor(ref$cluster)
  test<-as.data.frame(sp.integrated@assays$RNA@counts[rownames(sp.integrated@assays$RNA@counts)%in%genes.use,colnames(sp.integrated@assays$RNA@counts)%in%lastcell])
  test<-t(test)
  test<-as.data.frame(test)
  test$cluster<-sp.integrated@meta.data$seurat_clusters[rownames(sp.integrated@meta.data)%in%lastcell]
  colnames(test)<-gsub("-","_",colnames(test))
  test$cluster<-factor(test$cluster)
  library(rpart)
  library(randomForest)
  set.seed(3)
  data.ref<-randomForest(cluster~.,data=ref,mtry=length(genes.use),importance = T,proximity=T)
  print(data.ref)
  wp<-predict(data.ref,test,type='prob')
  res<-data.frame(name=character(length(lastcell)),predict=numeric(length(lastcell)),cluster=numeric(length(lastcell)))
  res$name<-rownames(test)
  res$predict<-wp
  res$cluster<-test$cluster

  wp3<-predict(data.ref,test)
  res3<-data.frame(name=character(length(lastcell)),predict=numeric(length(lastcell)),cluster=numeric(length(lastcell)))
  res3$name<-rownames(test)
  res3$predict<-wp3
  res3$cluster<-test$cluster
  result<-cbind(res,as.vector(res3$predict))
  auclist<-matrix(nrow=20,ncol=1)

  for (i in 1:ncol(res$predict)){
    cluster0<-cbind(res$predict[,i],as.vector(res3$cluster),as.vector(res3$predict))
    cluster<-cluster0[cluster0[,3]==i-1,]
    cluster<-as.data.frame(cluster)
    colnames(cluster)<-c("score","raw","predict")
    if (nrow(cluster)>0){
      cluster$label<-"0"
    for (j in 1:nrow(cluster)){
      if (cluster[j,2]==cluster[j,3]){
        cluster[j,4]=1}else{
          cluster[j,4]=0
        }
    }
    cluster[nrow(cluster)+1,]<-"0"
    pred <-prediction(as.numeric(cluster$score),as.numeric(cluster$label))
    perf <- performance(pred,"tpr","fpr")
    plot(perf, col='blue',lty=2)
    auc <- performance(pred,'auc')
    auc = unlist(slot(auc,"y.values"))
    picname<-paste("vst-roc-rep",k,"c",i-1,".pdf", sep='')
    picPath  <- paste(outPath,picname,sep='/')
    pdf(file=picPath)
    plot(perf,
         xlim=c(0,1), ylim=c(0,1),col='red',
         main=paste("ROC curve (", "AUC = ",auc,")"),
         lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
    abline(0,1)
    dev.off()
    auclist[i,]<-auc
    }else{
    auclist[i,]<-"0"  
    }
    

    
  }
  #names
  write.csv(result, file=out_filePath[k])
  write.csv(auclist, file=out_filePath2[k])


  E[k]<-mean(na.omit(auclist))

}
er<-unlist(E)
write.csv(er, "vst.auc.mean-2.csv")



