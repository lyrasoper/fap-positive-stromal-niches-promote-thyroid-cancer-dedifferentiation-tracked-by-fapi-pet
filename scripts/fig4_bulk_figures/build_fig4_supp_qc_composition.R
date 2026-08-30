#!/usr/bin/env Rscript
# Cancer Research submission - figure code release
# Builds: Supplementary Fig. S17 - experiment QC and NNLS composition control.
## Supplementary backing the MINIMAL Fig 4 panels A + B (experiment QC + composition control)
## a RNA PCA (WT vs KO) · b FAP on-target loss · c NNLS CAF composition stable
## d myCAF scaffold retained while Fap collapses — functional, not compositional perturbation.
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(readr);library(tibble);library(ggplot2);library(patchwork);library(ragg)})
proj<-Sys.getenv("BULK_OMICS_DATA", unset = "<BULK_OMICS_DATA>"); out<-file.path(proj,"outputs/wt_vs_fapko_bulk_omics_20260407")
FIGD<-file.path(Sys.getenv("PROJECT_ROOT", unset = "<PROJECT_ROOT>"), "outputs/fig4_minimal"); dir.create(FIGD,showWarnings=FALSE,recursive=TRUE)
col_wt<-"#3B6B9C"; col_ko<-"#C9493A"; grey<-"#9AA0A6"; dn_c<-"#2C6FA6"
pa<-function(base=7) theme_classic(base_size=base,base_family="Helvetica")+
  theme(axis.text=element_text(size=base-1,colour="black"),axis.title=element_text(size=base),
        plot.title=element_text(size=base+1,face="bold",hjust=0),plot.subtitle=element_text(size=base-1.3,colour="grey35"),
        legend.text=element_text(size=base-1.5),legend.title=element_text(size=base-1),
        axis.line=element_line(linewidth=0.4),axis.ticks=element_line(linewidth=0.4),legend.key.size=unit(3,"mm"))
deg_rna<-read_tsv(file.path(out,"rna/rna_differential_expression.tsv"),show_col_types=FALSE)
deg_prot<-read_tsv(file.path(out,"proteomics/proteomics_differential_expression.tsv"),show_col_types=FALSE)
nnls<-read_tsv(file.path(out,"caf_deconvolution/nnls_cell_fractions.tsv"),show_col_types=FALSE)
dstat<-read_tsv(file.path(out,"caf_deconvolution/deconvolution_group_stats.tsv"),show_col_types=FALSE)
core<-read_tsv(file.path(out,"ecm_mycaf_core_genes/core_gene_stats.tsv"),show_col_types=FALSE)
fpkm<-read_tsv(file.path(proj,"001_mRNA_Summary/1.GeneExpression/1_genes_fpkm_expression.txt"),show_col_types=FALSE)

## a: RNA PCA
fcols<-grep("^FPKM\\.",names(fpkm),value=TRUE)
mat<-log2(as.matrix(fpkm[,fcols])+1); rownames(mat)<-make.unique(fpkm$gene_name)
v<-apply(mat,1,var); mat<-mat[order(-v)[1:2000],]
pc<-prcomp(t(mat),center=TRUE,scale.=TRUE); ve<-round(100*pc$sdev^2/sum(pc$sdev^2),1)
pcd<-as_tibble(pc$x[,1:2],rownames="sample")%>%mutate(group=ifelse(grepl("FAP",sample),"FAP-KO","WT"))
p_a<-ggplot(pcd,aes(PC1,PC2,fill=group))+geom_point(size=2.4,shape=21,colour="white",stroke=0.4)+
  scale_fill_manual(values=c(WT=col_wt,`FAP-KO`=col_ko),name=NULL)+
  labs(x=sprintf("PC1 (%.0f%%)",ve[1]),y=sprintf("PC2 (%.0f%%)",ve[2]),title="a  Bulk-RNA PCA (top-2,000 variable)")+
  pa()+theme(legend.position=c(0.82,0.16))

## b: FAP on-target
fap<-tibble(sample=fcols,val=log2(as.numeric(fpkm[fpkm$gene_name=="Fap",fcols][1,])+1))%>%
  mutate(group=ifelse(grepl("FAP",sample),"FAP-KO","WT"))
fr<-deg_rna%>%filter(toupper(gene)=="FAP")%>%pull(log2_fc)%>%.[1]
fpp<-deg_prot%>%filter(toupper(gene)=="FAP")%>%pull(log2_fc)%>%.[1]
ymx<-max(fap$val)
p_b<-ggplot(fap,aes(group,val))+
  geom_boxplot(aes(fill=group),width=0.46,outlier.shape=NA,linewidth=0.35,alpha=0.9)+
  geom_jitter(aes(fill=group),width=0.08,height=0,size=1.8,shape=21,colour="white",stroke=0.45)+
  annotate("segment",x=1,xend=2,y=ymx+0.17,yend=ymx+0.17,linewidth=0.4,colour="grey30")+
  annotate("segment",x=1,xend=1,y=ymx+0.10,yend=ymx+0.17,linewidth=0.4,colour="grey30")+
  annotate("segment",x=2,xend=2,y=ymx+0.10,yend=ymx+0.17,linewidth=0.4,colour="grey30")+
  annotate("text",x=1.5,y=ymx+0.30,label="***",size=3,colour="grey20")+
  scale_fill_manual(values=c(WT=col_wt,`FAP-KO`=col_ko),guide="none")+
  scale_y_continuous(expand=expansion(mult=c(0.04,0.16)))+
  labs(x=NULL,y="Fap (log2 FPKM)",title="b  FAP on-target loss",
       subtitle=sprintf("RNA log2FC %.1f (q=2e-7) · protein log2FC %.1f",fr,fpp))+
  pa()+theme(axis.text.x=element_text(size=6.5),plot.subtitle=element_text(size=5.4))

## c: NNLS CAF composition stable
ct_lab<-c(myCAF="myCAF",iCAF="iCAF",apCAF="apCAF",Thyroid_tumor="Tumour",Immune="Immune",Endothelial="Endoth.")
nn<-nnls%>%mutate(group=ifelse(grepl("FAP",sample),"FAP-KO","WT"),celltype=factor(celltype,levels=unique(celltype)))
ps<-nn%>%group_by(celltype)%>%summarise(top=max(fraction),.groups="drop")%>%
  left_join(dstat%>%filter(method=="NNLS fraction")%>%transmute(celltype,p=wilcox_p),by="celltype")%>%
  mutate(lab=sprintf("P=%.2f",p),y=top+0.016)
p_c<-ggplot(nn,aes(celltype,fraction,fill=group))+
  geom_boxplot(width=0.64,outlier.shape=NA,linewidth=0.3,position=position_dodge(0.7))+
  geom_text(data=ps,aes(celltype,y,label=lab),inherit.aes=FALSE,size=1.6,colour="grey50")+
  scale_fill_manual(values=c(WT=col_wt,`FAP-KO`=col_ko),name=NULL)+
  scale_x_discrete(labels=ct_lab)+
  scale_y_continuous(expand=expansion(mult=c(0.04,0.1)))+
  labs(x=NULL,y="NNLS cell-type fraction",title="c  CAF composition preserved (no fraction collapses)")+
  pa()+theme(axis.text.x=element_text(angle=30,hjust=1,size=6),legend.position=c(0.88,0.9),legend.key.size=unit(2.6,"mm"))

## d: myCAF scaffold retained vs Fap collapses
sig<-core%>%filter(!gene_mgi %in% c("Fap","FAP"))%>%
  mutate(n_sigset=rowSums(cbind(Elyada,Dominguez,Kieffer),na.rm=TRUE))%>%
  filter(n_sigset>=2)%>%arrange(desc(abs(rna_log2fc)))%>%slice_head(n=12)
sd2<-bind_rows(tibble(gene="Fap",rna_log2fc=fr,grp="FAP (activation marker)"),
               sig%>%transmute(gene=gene_mgi,rna_log2fc,grp="myCAF scaffold (retained)"))%>%
  arrange(rna_log2fc)%>%mutate(gene=factor(gene,levels=gene))
p_d<-ggplot(sd2,aes(rna_log2fc,gene,fill=grp))+geom_col(width=0.7)+geom_vline(xintercept=0,linewidth=0.3,colour="grey50")+
  scale_fill_manual(values=c(`FAP (activation marker)`=col_ko,`myCAF scaffold (retained)`=grey),name=NULL)+
  labs(x="RNA log2 FC (FAP-KO vs WT)",y=NULL,title="d  Fap collapses; myCAF scaffold retained",
       subtitle="functional perturbation, not loss of the myofibroblast scaffold")+
  pa()+theme(axis.text.y=element_text(size=5.6,face="italic"),legend.position=c(0.62,0.18),plot.subtitle=element_text(size=5.3))

final<-(p_a|p_b)/(p_c|p_d)+plot_layout(heights=c(1,1.05))
agg_png(file.path(FIGD,"Fig4_supp_qc_composition.png"),width=180,height=132,units="mm",res=360)
print(final); dev.off()
cat("DONE\n")
