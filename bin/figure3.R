#!/usr/bin/env Rscript

## Copyright (C) 2020 by Julien Dorier, BioInformatics Competence Center, University of Lausanne, Switzerland.
## 
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

library(ggplot2)
library(cowplot)
library(data.table)
library(R.utils)
library(readxl)
set.seed(83712)

#############################################
##paths
#############################################
path.rnaseq.ctcf0.vs.wt="data/SupplementaryData4_RNAseq_CTCF0-WT.xlsx"
path.chipseq.ctcf.wt.vs.ctcf0="data/SupplementaryData1_ChIPseq_CTCF_WT-CTCF0.xlsx"
path.eigenvectors="data/SupplementaryData5_HiC_eigenvectors.xlsx"
path.blacklist="https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/dm6-blacklist.v2.bed.gz"
##blacklist:
##[1] Amemiya HM, Kundaje A, Boyle AP. The ENCODE blacklist: identification of problematic regions of the genome. Sci Rep. 2019 Dec; 9(1) 9354 DOI: 10.1038/s41598-019-45839-z
##[2] ENCODE Project Consortium. An integrated encyclopedia of DNA elements in the human genome. Nature. 2012 Sep 6;489(7414):57-74. doi: 10.1038/nature11247

dir.create("output/", showWarnings = FALSE,recursive=TRUE)

#############################################
##variables
#############################################
rnaseq.threshold.log2fc=log2(1.5)
rnaseq.threshold.padj=0.05
chrs=c("chr2L","chr2R","chr3L","chr3R","chr4","chrX")
plot.range=100000
binsize=2000 #bp

#############################################
##functions
#############################################
blacklisted=function(ch,position){                
    apply(data.table(ch,position),1,function(x){
        blacklist[,any(chr==x[1]&as.integer(x[2])>=start&as.integer(x[2])<=end)]
    })
}            

#############################################
##read data
#############################################
##CTCF peaks in WT
chipseq.ctcf.wt.vs.ctcf0=as.data.table(read_xlsx(path.chipseq.ctcf.wt.vs.ctcf0))
##keep only "up" peaks with positive fold change (filter out false positive)
chipseq.ctcf.wt.vs.ctcf0=chipseq.ctcf.wt.vs.ctcf0[direction=="up"&best.logFC>0]
##keep only selected chrs
chipseq.ctcf.wt.vs.ctcf0=chipseq.ctcf.wt.vs.ctcf0[chr%in%chrs]

##blacklisted regions
blacklist=fread(path.blacklist)[,1:3,with=FALSE]
setnames(blacklist,c("chr","start","end"))
setkey(blacklist,chr,start,end)

##Eigenvectors
eigenvectors=as.data.table(read_xlsx(path.eigenvectors))

##RNA seq
rnaseq.ctcf0.vs.wt=as.data.table(read_xlsx(path.rnaseq.ctcf0.vs.wt))
##convert Strand to numeric
rnaseq.ctcf0.vs.wt[,Strand:=ifelse(Strand=="+",+1,ifelse(Strand=="-",-1,0))]
##add TSS position
rnaseq.ctcf0.vs.wt[,TSS:=ifelse(Strand>0,Start,End)]
##keep only selected chrs
rnaseq.ctcf0.vs.wt=rnaseq.ctcf0.vs.wt[Chromosome%in%chrs]
##remove genes with TSS in blacklisted regions
rnaseq.ctcf0.vs.wt=rnaseq.ctcf0.vs.wt[!blacklisted(Chromosome,TSS)]
##differentially expressed genes
rnaseq.ctcf0.vs.wt.DE=rnaseq.ctcf0.vs.wt[abs(log2FoldChange)>=rnaseq.threshold.log2fc&padj<=rnaseq.threshold.padj]
##non differentially expressed genes
rnaseq.ctcf0.vs.wt.nonDE=rnaseq.ctcf0.vs.wt[!(abs(log2FoldChange)>=rnaseq.threshold.log2fc&padj<=rnaseq.threshold.padj)]

#############################################
## Sample same number of non-differentially
## expressed genes as differentially expressed
## genes with same distribution of mean expression
#############################################


##appoximate distribution of differentially expressed genes log average expression (baseMean=(w1118r1+w1118r2+w1118r3+ctcf0r1+ctcf0r2+ctcf0r3)/6.0)
DE.genes.dens=density(rnaseq.ctcf0.vs.wt.DE[,log(baseMean)])
DE.genes.dens_f=approxfun(DE.genes.dens$x,DE.genes.dens$y,yleft=0,yright=0)
nonDE.genes.dens=density(rnaseq.ctcf0.vs.wt.nonDE[,log(baseMean)])
nonDE.genes.dens_f=approxfun(nonDE.genes.dens$x,nonDE.genes.dens$y,yleft=0,yright=0)
prob_f=function(x){ifelse(nonDE.genes.dens_f(x)<1e-8,0,DE.genes.dens_f(x)/nonDE.genes.dens_f(x))}               
##Sample same number of nonDE genes as DE genes, with same distribution of log baseMean as DE.genes
rnaseq.ctcf0.vs.wt.nonDE.sample=rnaseq.ctcf0.vs.wt.nonDE[sample.int(nrow(rnaseq.ctcf0.vs.wt.nonDE),size=nrow(rnaseq.ctcf0.vs.wt.DE),replace=TRUE,prob=prob_f(rnaseq.ctcf0.vs.wt.nonDE[,log(baseMean)]))]
##Note: method used to generate random integers in sample.int() changed in R 3.6.0.

    
#############################################
## CTCF around differentially expressed genes
#############################################

## CTCF around DE gene TSSs 
ctcf.around.DEgenes=rbindlist(lapply(1:nrow(rnaseq.ctcf0.vs.wt.DE),function(i){
    strand1=rnaseq.ctcf0.vs.wt.DE[i,Strand]
    chr1=rnaseq.ctcf0.vs.wt.DE[i,Chromosome]
    position1=rnaseq.ctcf0.vs.wt.DE[i,TSS]
    tmp=chipseq.ctcf.wt.vs.ctcf0[chr==chr1&best.pos>=position1-plot.range-binsize&best.pos<=position1+plot.range+binsize] 
    tmp[,distance.binned:=strand1*round((best.pos-position1)/binsize)*binsize] ##distance measured in the direction of gene transcription
    if(nrow(tmp)==0)
        return(cbind(i,data.table(distance.binned=NaN)))
    cbind(i,tmp[,.(distance.binned)])
}))
nb.DEgenes=nrow(rnaseq.ctcf0.vs.wt.DE)
##fraction of genes with AT LEAST one ctcf at given distance (=> unique())
ctcf.around.DEgenes=unique(ctcf.around.DEgenes)[,.(fraction=.N/nb.DEgenes),by=distance.binned]
##add missing distances
distance.binned.all=data.table(distance.binned=seq(-ceiling(plot.range/binsize)*binsize,plot.range,by=binsize))
ctcf.around.DEgenes=ctcf.around.DEgenes[distance.binned.all,on="distance.binned"]
ctcf.around.DEgenes[is.na(fraction),fraction:=0]

## CTCF around non DE gene TSSs (same number as DE genes)
ctcf.around.nonDEgenes=rbindlist(lapply(1:nrow(rnaseq.ctcf0.vs.wt.nonDE.sample),function(i){
    strand1=rnaseq.ctcf0.vs.wt.nonDE.sample[i,Strand]
    chr1=rnaseq.ctcf0.vs.wt.nonDE.sample[i,Chromosome]
    position1=rnaseq.ctcf0.vs.wt.nonDE.sample[i,TSS]
    tmp=chipseq.ctcf.wt.vs.ctcf0[chr==chr1&best.pos>=position1-plot.range-binsize&best.pos<=position1+plot.range+binsize] 
    tmp[,distance.binned:=strand1*round((best.pos-position1)/binsize)*binsize] ##distance measured in the direction of gene transcription
    if(nrow(tmp)==0)
        return(cbind(i,data.table(distance.binned=NaN)))
    cbind(i,tmp[,.(distance.binned)])
}))
nb.nonDEgenes=nrow(rnaseq.ctcf0.vs.wt.nonDE.sample)
##fraction of genes with AT LEAST one ctcf at given distance (=> unique())
ctcf.around.nonDEgenes=unique(ctcf.around.nonDEgenes)[,.(fraction=.N/nb.nonDEgenes),by=distance.binned]
##add missing distances
distance.binned.all=data.table(distance.binned=seq(-ceiling(plot.range/binsize)*binsize,plot.range,by=binsize))
ctcf.around.nonDEgenes=ctcf.around.nonDEgenes[distance.binned.all,on="distance.binned"]
ctcf.around.nonDEgenes[is.na(fraction),fraction:=0]



##plot
pdf("output/figure_3f.pdf",10,5)
ctcf.around.DEgenes[,label:=paste0("DE in CTCF0\n(n=",nrow(rnaseq.ctcf0.vs.wt.DE)," genes)")]
ctcf.around.nonDEgenes[,label:=paste0("non-DE in CTCF0\n(n=",nrow(rnaseq.ctcf0.vs.wt.nonDE.sample)," genes)")]
alpha=c(0.3,1)
names(alpha)=c(paste0("non-DE in CTCF0\n(n=",nrow(rnaseq.ctcf0.vs.wt.nonDE.sample)," genes)"),paste0("DE in CTCF0\n(n=",nrow(rnaseq.ctcf0.vs.wt.DE)," genes)"))
fraction.DE.at.0=ctcf.around.DEgenes[distance.binned==0,fraction]
fraction.nonDE.average=ctcf.around.nonDEgenes[,mean(fraction)]
annotation=paste0(round(fraction.DE.at.0*100),"% (",round(fraction.DE.at.0/fraction.nonDE.average),"x enrichment)")
p=ggplot(aes(x=distance.binned/1000,y=fraction,alpha=label),data=rbind(ctcf.around.DEgenes,ctcf.around.nonDEgenes))
p=p+geom_line()
p=p+scale_alpha_manual("",values=alpha)
p=p+scale_x_continuous("Distance to gene TSS (kb)",limits=c(-plot.range/1000,plot.range/1000))
p=p+scale_y_continuous(paste0("% of genes with CTCF peak\nat a given distance (per ",binsize/1000,"kb)"),labels = scales::percent)
p=p+ggtitle("Enrichment of CTCF peaks around DE gene TSSs")
p=p+annotate("text",x=0,y=ctcf.around.DEgenes[,max(fraction,na.rm=TRUE)],label=annotation,hjust=-0.1,vjust=1)
p=p+theme_bw()
print(p)
dev.off()

##save data
fwrite(rbind(ctcf.around.DEgenes,ctcf.around.nonDEgenes)[,.(type=gsub("\n"," ",label),distance=format(distance.binned,scientific=FALSE,trim=TRUE),fraction)],file="output/figure_3f.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE,na="NA")

#############################################
## Differentially expressed genes around CTCF
#############################################

## DE gene TSSs around CTCF 
DEgenes.around.ctcf=rbindlist(lapply(1:nrow(chipseq.ctcf.wt.vs.ctcf0),function(i){
    chr1=chipseq.ctcf.wt.vs.ctcf0[i,chr]
    position1=chipseq.ctcf.wt.vs.ctcf0[i,best.pos]
    tmp=rnaseq.ctcf0.vs.wt.DE[Chromosome==chr1&TSS>=position1-plot.range-binsize&TSS<=position1+plot.range+binsize] 
    tmp[,distance.binned:=Strand*round((TSS-position1)/binsize)*binsize]
    if(nrow(tmp)==0)
        return(cbind(i,data.table(distance.binned=NaN)))
    cbind(i,tmp[,.(distance.binned)])
}))
nb.ctcf=nrow(chipseq.ctcf.wt.vs.ctcf0)
##fraction of CTCF with AT LEAST one DE gene at given distance (=> unique())
DEgenes.around.ctcf=unique(DEgenes.around.ctcf)[,.(fraction=.N/nb.ctcf),by=distance.binned]
##add missing distances
distance.binned.all=data.table(distance.binned=seq(-ceiling(plot.range/binsize)*binsize,plot.range,by=binsize))
DEgenes.around.ctcf=DEgenes.around.ctcf[distance.binned.all,on="distance.binned"]
DEgenes.around.ctcf[is.na(fraction),fraction:=0]


## non DE gene TSSs (same number as DE genes) around CTCF 
nonDEgenes.around.ctcf=rbindlist(lapply(1:nrow(chipseq.ctcf.wt.vs.ctcf0),function(i){
    chr1=chipseq.ctcf.wt.vs.ctcf0[i,chr]
    position1=chipseq.ctcf.wt.vs.ctcf0[i,best.pos]
    tmp=rnaseq.ctcf0.vs.wt.nonDE.sample[Chromosome==chr1&TSS>=position1-plot.range-binsize&TSS<=position1+plot.range+binsize] 
    tmp[,distance.binned:=Strand*round((TSS-position1)/binsize)*binsize]
    if(nrow(tmp)==0)
        return(cbind(i,data.table(distance.binned=NaN)))
    cbind(i,tmp[,.(distance.binned)])
}))
##fraction of CTCF with AT LEAST one non DE gene at given distance (=> unique())
nonDEgenes.around.ctcf=unique(nonDEgenes.around.ctcf)[,.(fraction=.N/nb.ctcf),by=distance.binned]
##add missing distances
distance.binned.all=data.table(distance.binned=seq(-ceiling(plot.range/binsize)*binsize,plot.range,by=binsize))
nonDEgenes.around.ctcf=nonDEgenes.around.ctcf[distance.binned.all,on="distance.binned"]
nonDEgenes.around.ctcf[is.na(fraction),fraction:=0]


##plot
pdf("output/figure_3g.pdf",10,5)
DEgenes.around.ctcf[,label:=paste0("DE in CTCF0\n(n=",nrow(rnaseq.ctcf0.vs.wt.DE)," genes)")]
nonDEgenes.around.ctcf[,label:=paste0("non-DE in CTCF0\n(n=",nrow(rnaseq.ctcf0.vs.wt.nonDE.sample)," genes)")]
alpha=c(0.3,1)
names(alpha)=c(paste0("non-DE in CTCF0\n(n=",nrow(rnaseq.ctcf0.vs.wt.nonDE.sample)," genes)"),paste0("DE in CTCF0\n(n=",nrow(rnaseq.ctcf0.vs.wt.DE)," genes)"))
fraction.DE.at.0=DEgenes.around.ctcf[distance.binned==0,fraction]
fraction.nonDE.average=nonDEgenes.around.ctcf[,mean(fraction)]
annotation=paste0(round(fraction.DE.at.0*100),"% (",round(fraction.DE.at.0/fraction.nonDE.average),"x enrichment)")
p=ggplot(aes(x=distance.binned/1000,y=fraction,alpha=label),data=rbind(DEgenes.around.ctcf,nonDEgenes.around.ctcf))
p=p+geom_line()
p=p+scale_alpha_manual("",values=alpha)
p=p+scale_x_continuous("Distance to CTCF peak (kb)",limits=c(-plot.range/1000,plot.range/1000))
p=p+scale_y_continuous(paste0("% of CTCF peaks with a gene TSS\nat a given distance (per ",binsize/1000,"kb)"),labels = scales::percent)
p=p+ggtitle("Enrichment of DE gene TSSs around CTCF peaks")
p=p+annotate("text",x=0,y=DEgenes.around.ctcf[,max(fraction,na.rm=TRUE)],label=annotation,hjust=-0.1,vjust=1)
p=p+theme_bw()
print(p)
dev.off()

##save data
fwrite(rbind(DEgenes.around.ctcf,nonDEgenes.around.ctcf)[,.(type=gsub("\n"," ",label),distance=format(distance.binned,scientific=FALSE,trim=TRUE),fraction)],file="output/figure_3g.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE,na="NA")




#############################################
## Sample same number of non-differentially
## expressed genes as differentially expressed (with logFC>0)
## genes with same distribution of mean expression
#############################################
plist=lapply(c("DOWN","UP"),function(direction){
    ##differentially expressed genes (in the specified direction)
    if(direction=="UP")
        rnaseq.ctcf0.vs.wt.DE=rnaseq.ctcf0.vs.wt[abs(log2FoldChange)>=rnaseq.threshold.log2fc&padj<=rnaseq.threshold.padj&log2FoldChange>0]
    if(direction=="DOWN")
        rnaseq.ctcf0.vs.wt.DE=rnaseq.ctcf0.vs.wt[abs(log2FoldChange)>=rnaseq.threshold.log2fc&padj<=rnaseq.threshold.padj&log2FoldChange<0]

    ##appoximate distribution of differentially expressed genes (in the specified direction) log average expression (baseMean)
    DE.genes.dens=density(rnaseq.ctcf0.vs.wt.DE[,log(baseMean)])
    DE.genes.dens_f=approxfun(DE.genes.dens$x,DE.genes.dens$y,yleft=0,yright=0)
    nonDE.genes.dens=density(rnaseq.ctcf0.vs.wt.nonDE[,log(baseMean)])
    nonDE.genes.dens_f=approxfun(nonDE.genes.dens$x,nonDE.genes.dens$y,yleft=0,yright=0)
    prob_f=function(x){ifelse(nonDE.genes.dens_f(x)<1e-8,0,DE.genes.dens_f(x)/nonDE.genes.dens_f(x))}               
    ##Sample same number of nonDE genes as DE genes, with same distribution of log baseMean as DE.genes
    rnaseq.ctcf0.vs.wt.nonDE.sample=rnaseq.ctcf0.vs.wt.nonDE[sample.int(nrow(rnaseq.ctcf0.vs.wt.nonDE),size=nrow(rnaseq.ctcf0.vs.wt.DE),replace=TRUE,prob=prob_f(rnaseq.ctcf0.vs.wt.nonDE[,log(baseMean)]))]
    ##Note: method used to generate random integers in sample.int() changed in R 3.6.0.

    
    #############################################
    ## CTCF around differentially expressed genes
    #############################################

    ## CTCF around DE gene TSSs 
    ctcf.around.DEgenes=rbindlist(lapply(1:nrow(rnaseq.ctcf0.vs.wt.DE),function(i){
        strand1=rnaseq.ctcf0.vs.wt.DE[i,Strand]
        chr1=rnaseq.ctcf0.vs.wt.DE[i,Chromosome]
        position1=rnaseq.ctcf0.vs.wt.DE[i,TSS]
        tmp=chipseq.ctcf.wt.vs.ctcf0[chr==chr1&best.pos>=position1-plot.range-binsize&best.pos<=position1+plot.range+binsize] 
        tmp[,distance.binned:=strand1*round((best.pos-position1)/binsize)*binsize] ##distance measured in the direction of gene transcription
        if(nrow(tmp)==0)
            return(cbind(i,data.table(distance.binned=NaN)))
        cbind(i,tmp[,.(distance.binned)])
    }))
    nb.DEgenes=nrow(rnaseq.ctcf0.vs.wt.DE)
    ##fraction of genes with AT LEAST one CTCF at given distance (=> unique())
    ctcf.around.DEgenes=unique(ctcf.around.DEgenes)[,.(fraction=.N/nb.DEgenes),by=distance.binned]
    ##add missing distances
    distance.binned.all=data.table(distance.binned=seq(-ceiling(plot.range/binsize)*binsize,plot.range,by=binsize))
    ctcf.around.DEgenes=ctcf.around.DEgenes[distance.binned.all,on="distance.binned"]
    ctcf.around.DEgenes[is.na(fraction),fraction:=0]

    ## CTCF around non DE gene TSSs (same number as DE genes)
    ctcf.around.nonDEgenes=rbindlist(lapply(1:nrow(rnaseq.ctcf0.vs.wt.nonDE.sample),function(i){
        strand1=rnaseq.ctcf0.vs.wt.nonDE.sample[i,Strand]
        chr1=rnaseq.ctcf0.vs.wt.nonDE.sample[i,Chromosome]
        position1=rnaseq.ctcf0.vs.wt.nonDE.sample[i,TSS]
        tmp=chipseq.ctcf.wt.vs.ctcf0[chr==chr1&best.pos>=position1-plot.range-binsize&best.pos<=position1+plot.range+binsize] 
        tmp[,distance.binned:=strand1*round((best.pos-position1)/binsize)*binsize] ##distance measured in the direction of gene transcription
        if(nrow(tmp)==0)
            return(cbind(i,data.table(distance.binned=NaN)))
        cbind(i,tmp[,.(distance.binned)])
    }))
    nb.nonDEgenes=nrow(rnaseq.ctcf0.vs.wt.nonDE.sample)
    ##fraction of genes with AT LEAST one CTCF at given distance (=> unique())
    ctcf.around.nonDEgenes=unique(ctcf.around.nonDEgenes)[,.(fraction=.N/nb.nonDEgenes),by=distance.binned]
    ##add missing distances
    distance.binned.all=data.table(distance.binned=seq(-ceiling(plot.range/binsize)*binsize,plot.range,by=binsize))
    ctcf.around.nonDEgenes=ctcf.around.nonDEgenes[distance.binned.all,on="distance.binned"]
    ctcf.around.nonDEgenes[is.na(fraction),fraction:=0]



    ##plot
    ctcf.around.DEgenes[,label:=paste0(direction," in CTCF0\n(n=",nrow(rnaseq.ctcf0.vs.wt.DE)," genes)")]
    ctcf.around.nonDEgenes[,label:=paste0("non-DE in CTCF0\n(n=",nrow(rnaseq.ctcf0.vs.wt.nonDE.sample)," genes)")]
    alpha=c(0.3,1)
    names(alpha)=c(paste0("non-DE in CTCF0\n(n=",nrow(rnaseq.ctcf0.vs.wt.nonDE.sample)," genes)"),paste0(direction," in CTCF0\n(n=",nrow(rnaseq.ctcf0.vs.wt.DE)," genes)"))
    fraction.DE.at.0=ctcf.around.DEgenes[distance.binned==0,fraction]
    fraction.nonDE.average=ctcf.around.nonDEgenes[,mean(fraction)]
    annotation=paste0(round(fraction.DE.at.0*100),"% (",round(fraction.DE.at.0/fraction.nonDE.average),"x enrichment)")
    p=ggplot(aes(x=distance.binned/1000,y=fraction,alpha=label),data=rbind(ctcf.around.DEgenes,ctcf.around.nonDEgenes))
    p=p+geom_line()
    p=p+scale_alpha_manual("",values=alpha)
    p=p+scale_x_continuous("Distance to gene TSS (kb)",limits=c(-plot.range/1000,plot.range/1000))
    p=p+scale_y_continuous(paste0("% of genes with CTCF peak\nat a given distance (per ",binsize/1000,"kb)"),labels = scales::percent)
    p=p+ggtitle(paste0("Enrichment of CTCF peaks around ",direction," gene TSSs"))
    p=p+annotate("text",x=0,y=ctcf.around.DEgenes[,max(fraction,na.rm=TRUE)],label=annotation,hjust=-0.1,vjust=1)
    p=p+theme_bw()

    ##save data
    fwrite(rbind(ctcf.around.DEgenes,ctcf.around.nonDEgenes)[,.(type=gsub("\n"," ",label),distance=format(distance.binned,scientific=FALSE,trim=TRUE),fraction)],file="output/figure_S3c.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=(direction=="DOWN"),na="NA",append=(direction=="UP"))

    p
})
pdf("output/figure_S3c.pdf",10,10)
print(plot_grid(plotlist=plist,ncol=1,align="v",axis="lr"))
dev.off()


#############################################
## Eigenvectors
#############################################
##ignore chr4
eigenvectors=eigenvectors[chr!="chr4"]
setkey(eigenvectors,chr,bin.start)
##DE genes: add eigenvector
rnaseq.ctcf0.vs.wt.DE.UP=rnaseq.ctcf0.vs.wt.DE[log2FoldChange>0,.(Chromosome,TSS.bin.start=as.integer(floor(TSS/binsize)*binsize))]
setkey(rnaseq.ctcf0.vs.wt.DE.UP,Chromosome,TSS.bin.start)
rnaseq.ctcf0.vs.wt.DE.UP=eigenvectors[,.(chr,bin.start,eigenvector.WT,eigenvector.CTCF0)][rnaseq.ctcf0.vs.wt.DE.UP]
rnaseq.ctcf0.vs.wt.DE.DOWN=rnaseq.ctcf0.vs.wt.DE[log2FoldChange<0,.(Chromosome,TSS.bin.start=as.integer(floor(TSS/binsize)*binsize))]
setkey(rnaseq.ctcf0.vs.wt.DE.DOWN,Chromosome,TSS.bin.start)
rnaseq.ctcf0.vs.wt.DE.DOWN=eigenvectors[,.(chr,bin.start,eigenvector.WT,eigenvector.CTCF0)][rnaseq.ctcf0.vs.wt.DE.DOWN]


pdf("output/figure_S3d.pdf",8,8)
smoothScatter(eigenvectors[,.(eigenvector.WT,eigenvector.CTCF0)],asp=1,bandwidth=0.0001,nbin=512,xlab="Eigenvector in WT",ylab="Eigenvector in CTCF0",nrpoints=0,main=paste0("CTCF0 versus WT eigenvector values for all ",binsize/1000,"kb bins in chr2L/R, chr3L/R and chrX (blue density plot)\nDE genes UP in CTCF0 (red) and DOWN in CTCF0 (green)." ),cex.main=0.8,useRaster=TRUE)
abline(v=0)
abline(h=0)
##add genes
points(rnaseq.ctcf0.vs.wt.DE.UP[,.(eigenvector.WT,eigenvector.CTCF0)],col="red",cex=0.5,pch=16)
points(rnaseq.ctcf0.vs.wt.DE.DOWN[,.(eigenvector.WT,eigenvector.CTCF0)],col="green",cex=0.5,pch=16)
dev.off()

##save data
fwrite(rbind(rnaseq.ctcf0.vs.wt.DE.UP[,.(type="DE gene (WT<CTCF0)",eigenvector.WT,eigenvector.CTCF0)],
      rnaseq.ctcf0.vs.wt.DE.DOWN[,.(type="DE gene (WT>CTCF0)",eigenvector.WT,eigenvector.CTCF0)],
      eigenvectors[,.(type=paste0("All ",binsize/1000,"kb bins"),eigenvector.WT,eigenvector.CTCF0)]),
      file="output/figure_S3d.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE,na="NA")
      
