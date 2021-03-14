#reproducibility
filename<-list.files(pattern = "random-genelist*")
path <- getwd()
filePath <- sapply(filename, function(x){
  paste(path,x,sep='/')})
data <- lapply(filePath, function(x){
  read.csv(x)})
a<-list()
for (i in 1:9){
  a[[i]]<-matrix(nrow=2000,ncol=50)
  a[[i]]<-as.data.frame(a[[i]])
}

for (i in 1:9){
  for (j in 1:50){
    a[[i]][,j]<-data[[j]][,i+1]
  }
}


b<-list()
for (i in 1:9){
  b[[i]]<-matrix(nrow=50,ncol=50)
  b[[i]]<-as.data.frame(b[[i]])
}
for  (k in 1:9) {
  for (i in 1:50){
    for (j in 1:50){
      b[[k]][i,j]<-sum(a[[k]][,i]%in%a[[k]][,j])
    }
  }
}

repro<-list()
for (i in 1:9){
  repro[[i]]<-matrix(nrow=50,ncol=50)
  repro[[i]]<-as.data.frame(repro[[i]])
}
for  (k in 1:9) {
  for (i in 1:49){
    for (j in (i+1):50){
      repro[[k]][i,j]<-sum(a[[k]][,i]%in%a[[k]][,j])
    }
  }
}
reproducity<-matrix(ncol=9,nrow=25*49)
reproducity<-as.data.frame(reproducity)
for (i in 1:9){
  reproducity[,i]<-na.omit(unlist(repro[[i]]))
}
colnames(reproducity)<-c("Scmap","Scran","Seurat-SCT","SCHS","Rogue","Rogue_n","Seurat-vst","Seurat-disp","M3Drop")
reproducity<-reproducity/2000
library(reshape2)
repro_plot<-melt(reproducity)
colnames(repro_plot)<-c("method","reproducity")
repro_plot$method<-factor(repro_plot$method,levels = c("M3Drop","Seurat-vst","Seurat-SCT","Rogue","Seurat-disp","Rogue_n","Scmap","Scran","SCHS"))
library(ggplot2)
ggplot(repro_plot,aes(method,reproducity))+geom_boxplot()+labs(title="HSPC")
