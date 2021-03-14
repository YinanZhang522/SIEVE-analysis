#purity
library(ROGUE)
purity_vst<-rogue(gene_expr_norm,label=label$vst,platform = "full-length",sample)
purity_rogue<-rogue(gene_expr,label=label$rogue,platform = "full-length",sample)
purity_sct<-rogue(gene_expr_norm,label=label$sct,platform = "full-length",sample)
purity_rogue_norm<-rogue(gene_expr_norm,label=label$rogue_norm,platform = "full-length",sample)
purity_hvg<-rogue(gene_expr,label=label$hvg,platform = "full-length",sample)
purity_scmap<-rogue(gene_expr,label=label$scmap,platform = "full-length",sample)
purity_scran<-rogue(gene_expr,label=label$scran,platform = "full-length",sample)
purity_schs<-rogue(gene_expr_norm,label=label$schs,platform = "full-length",sample)
purity_disp<-rogue(gene_expr_norm,label=label$disp,platform = "full-length",sample)

purity<-cbind(purity_disp,purity_hvg,purity_rogue,purity_rogue_norm,purity_schs,purity_scmap,purity_scran,purity_sct,purity_vst)
purity<-as.vector(as.matrix(purity))
purity<-data.frame(method=character(64),score=purity)
purity$method=rep(c("disp","hvg","rogue","rogue_norm","schs","scmap","scran","sct","vst"),c(7,8,7,7,6,7,7,8,7))
purity$score[1:7]<-purity_disp
purity$score[8:15]<-purity_hvg
purity$score[16:22]<-purity_rogue
purity$score[23:29]<-purity_rogue_norm
purity$score[30:35]<-purity_schs
purity$score[36:42]<-purity_scmap
purity$score[43:49]<-purity_scran
purity$score[50:57]<-purity_sct
purity$score[58:64]<-purity_vst
#colnames(purity)<-c("method","score")
purity$score<-as.numeric(purity$score)
purity$method<-factor(purity$method,levels=c("hvg","vst","sct","disp","rogue","rogue_norm","scmap","scran","schs"))
pic=ggplot(purity,aes(method,score))+geom_violin()+theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
pic
ggsave("spBM1-purity.pdf", pic, width=12, height=10, units="cm")

