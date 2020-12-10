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


library(reshape2)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(data.table)
library(R.utils)
library(readxl)
set.seed(83712)

#############################################
##paths
#############################################
path.CD.boundaries="data/SupplementaryData3_HiC_CD_boundaries.xlsx"
path.insulation.scores="data/SupplementaryData2_HiC_insulation_scores.xlsx"
path.chipseq.ctcf.wt.vs.ctcf0="data/SupplementaryData1_ChIPseq_CTCF_WT-CTCF0.xlsx"
path.chipseq.cp190.wt.vs.ctcf0="data/SupplementaryData7_ChIPseq_Cp190_WT-CTCF0.xlsx"
path.chipseq.cp190.ctcf0.vs.cp190KO="data/SupplementaryData8_ChIPseq_Cp190_CTCF0-Cp190KO.xlsx"
path.chipseq.cp190.wt.vs.cp190KO="data/SupplementaryData6_ChIPseq_Cp190_WT-Cp190KO.xlsx"
path.rnaseq.ctcf0.vs.wt="data/SupplementaryData4_RNAseq_CTCF0-WT.xlsx"
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
##chromosome sizes [dos Santos et al. (2015) Nucleic Acids Research]
chr.sizes=data.table(chr=c("chr2L","chr2R","chr3L","chr3R","chr4","chrX"),
                     length=c(23513712L,25286936L,28110227L,32079331L,1348131L,23542271L))


#############################################
##functions
#############################################
chr.sizes[,cumulative.start:=0L]
chr.sizes[,cumulative.end:=cumsum(length)]
chr.sizes[-1,cumulative.start:=chr.sizes[-nrow(chr.sizes),cumulative.end]]
##random positions with uniform distribution on the genome
get_random_positions=function(n){
    total.length=chr.sizes[,sum(length)]
    xi=as.integer(floor(runif(n,min=1,max=total.length+1)))
    rbindlist(lapply(xi,function(x){
        chr.sizes[x>cumulative.start&x<=cumulative.end,.(chr,position=x-cumulative.start)]
    }))   
}

blacklisted=function(ch,position){                
    apply(data.table(ch,position),1,function(x){
        blacklist[,any(chr==x[1]&as.integer(x[2])>=start&as.integer(x[2])<=end)]
    })
}            

#############################################
##read data
#############################################
CD.boundaries=as.data.table(read_xlsx(path.CD.boundaries))
insulation.scores=as.data.table(read_xlsx(path.insulation.scores))

##CTCF peaks in WT (WT-CTCF0)
chipseq.ctcf.wt.vs.ctcf0=as.data.table(read_xlsx(path.chipseq.ctcf.wt.vs.ctcf0))
##keep only "up" peaks with positive fold change (filter out false positive)
chipseq.ctcf.wt.vs.ctcf0=chipseq.ctcf.wt.vs.ctcf0[direction=="up"&best.logFC>0]
##keep only selected chrs
chipseq.ctcf.wt.vs.ctcf0=chipseq.ctcf.wt.vs.ctcf0[chr%in%chrs]

##Cp190 peaks in WT (WT-Cp190KO)
chipseq.cp190.wt.vs.cp190KO=as.data.table(read_xlsx(path.chipseq.cp190.wt.vs.cp190KO))
##keep only "up" peaks with positive fold change (filter out false positive)
chipseq.cp190.wt.vs.cp190KO=chipseq.cp190.wt.vs.cp190KO[direction=="up"&best.logFC>0]
##keep only selected chrs
chipseq.cp190.wt.vs.cp190KO=chipseq.cp190.wt.vs.cp190KO[chr%in%chrs]

##Cp190 peaks in CTCF0 (CTCF0-Cp190KO)
chipseq.cp190.ctcf0.vs.cp190KO=as.data.table(read_xlsx(path.chipseq.cp190.ctcf0.vs.cp190KO))
##keep only "up" peaks with positive fold change (filter out false positive)
chipseq.cp190.ctcf0.vs.cp190KO=chipseq.cp190.ctcf0.vs.cp190KO[direction=="up"&best.logFC>0]
##keep only selected chrs
chipseq.cp190.ctcf0.vs.cp190KO=chipseq.cp190.ctcf0.vs.cp190KO[chr%in%chrs]

##Cp190 differential binding (WT-CTCF0)
chipseq.cp190.wt.vs.ctcf0=as.data.table(read_xlsx(path.chipseq.cp190.wt.vs.ctcf0))
##keep only selected chrs
chipseq.cp190.wt.vs.ctcf0=chipseq.cp190.wt.vs.ctcf0[chr%in%chrs]


##blacklisted regions
blacklist=fread(path.blacklist)[,1:3,with=FALSE]
setnames(blacklist,c("chr","start","end"))
setkey(blacklist,chr,start,end)

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


#############################################
##ChIPseq andRNAseq around selected CD boundaries
#############################################

##add nearest CTCF peaks in WT (within +/-binsize) to CD.boundaries
CD.boundaries=rbindlist(lapply(1:nrow(CD.boundaries),function(i){
    chr1=CD.boundaries[i,chr]
    position1=CD.boundaries[i,position]
    ##CTCF peaks within +/-binsize
    tmp=chipseq.ctcf.wt.vs.ctcf0[chr==chr1&abs(best.pos-position1)<=binsize]
    if(nrow(tmp)==0)
        return(cbind(CD.boundaries[i,],data.table(CTCF.position=NaN,CTCF.start=NaN,CTCF.end=NaN,nb.neighboring.CTCF=0)))
    cbind(CD.boundaries[i,],tmp[,.(CTCF.position=best.pos[which.min(best.pos-position1)],CTCF.start=start[which.min(best.pos-position1)],CTCF.end=end[which.min(best.pos-position1)],nb.neighboring.CTCF=.N)])
}))

##add boundary.type
CD.boundaries[exist.wildtype==TRUE&exist.CTCF0==TRUE,boundary.type:="Common"]
CD.boundaries[exist.wildtype==TRUE&exist.CTCF0==FALSE,boundary.type:="WT-specific"]
CD.boundaries[exist.wildtype==FALSE&exist.CTCF0==TRUE,boundary.type:="CTCF0-specific"]

##order CD boundaries by delta physical insulation score
CD.boundaries[,score.diff:=score.CTCF0-score.wildtype]
setorder(CD.boundaries,score.diff)
CD.boundaries[,i:=1:.N]
setkey(CD.boundaries,i)

for(CD.boundaries.type in c("All","CTCF-bound"))
{
    if(CD.boundaries.type=="All")
        CD.boundaries.selected=CD.boundaries[,.(i,chr,position,boundary.type,score.diff)]
    if(CD.boundaries.type=="CTCF-bound")
        CD.boundaries.selected=CD.boundaries[nb.neighboring.CTCF>0,.(i,chr,position=CTCF.position,boundary.type,score.diff)]
        
    ## CTCF peaks in WT around selected CD boundaries 
    ctcf.wt.vs.ctcf0.around.boundaries=rbindlist(lapply(1:nrow(CD.boundaries.selected),function(i){
        chr1=CD.boundaries.selected[i,chr]
        position1=CD.boundaries.selected[i,position]
        tmp=chipseq.ctcf.wt.vs.ctcf0[chr==chr1&best.pos>=position1-plot.range-binsize&best.pos<=position1+plot.range+binsize] 
        tmp[,distance:=(best.pos-position1)]
        if(nrow(tmp)==0)
            return(cbind(i,data.table(distance=NaN,best.logFC=NaN)))
        cbind(i,tmp[,.(distance,best.logFC)])
    }))
    ## Cp190 peaks in WT around selected CD boundaries 
    cp190.wt.vs.cp190KO.around.boundaries=rbindlist(lapply(1:nrow(CD.boundaries.selected),function(i){
        chr1=CD.boundaries.selected[i,chr]
        position1=CD.boundaries.selected[i,position]
        tmp=chipseq.cp190.wt.vs.cp190KO[chr==chr1&best.pos>=position1-plot.range-binsize&best.pos<=position1+plot.range+binsize] 
        tmp[,distance:=(best.pos-position1)]
        if(nrow(tmp)==0)
            return(cbind(i,data.table(distance=NaN,best.logFC=NaN)))
        cbind(i,tmp[,.(distance,best.logFC)])
    }))
    ## Cp190 peaks in CTCF0 around selected CD boundaries 
    cp190.ctcf0.vs.cp190KO.around.boundaries=rbindlist(lapply(1:nrow(CD.boundaries.selected),function(i){
        chr1=CD.boundaries.selected[i,chr]
        position1=CD.boundaries.selected[i,position]
        tmp=chipseq.cp190.ctcf0.vs.cp190KO[chr==chr1&best.pos>=position1-plot.range-binsize&best.pos<=position1+plot.range+binsize] 
        tmp[,distance:=(best.pos-position1)]
        if(nrow(tmp)==0)
            return(cbind(i,data.table(distance=NaN,best.logFC=NaN)))
        cbind(i,tmp[,.(distance,best.logFC)])
    }))
    ## Cp190 differential binding regions (CTCF0-WT) around selected CD boundaries 
    cp190.wt.vs.ctcf0.around.boundaries=rbindlist(lapply(1:nrow(CD.boundaries.selected),function(i){
        chr1=CD.boundaries.selected[i,chr]
        position1=CD.boundaries.selected[i,position]
        tmp=chipseq.cp190.wt.vs.ctcf0[chr==chr1&best.pos>=position1-plot.range-binsize&best.pos<=position1+plot.range+binsize] 
        tmp[,distance:=(best.pos-position1)]
        if(nrow(tmp)==0)
            return(cbind(i,data.table(distance=NaN,best.logFC=NaN)))
        cbind(i,tmp[,.(distance,best.logFC)])
    }))
    ## RNAseq around selected CD boundaries 
    rnaseq.ctcf0.vs.wt.around.boundaries=rbindlist(lapply(1:nrow(CD.boundaries.selected),function(i){
        chr1=CD.boundaries.selected[i,chr]
        position1=CD.boundaries.selected[i,position]
        tmp=rnaseq.ctcf0.vs.wt[Chromosome==chr1&TSS>=position1-plot.range-binsize&TSS<=position1+plot.range+binsize] 
        tmp[,distance:=(TSS-position1)]
        if(nrow(tmp)==0)
            return(cbind(i,data.table(distance=NaN,padj=NaN,log2FoldChange=NaN)))
        cbind(i,tmp[,.(distance,padj,log2FoldChange)])
    }))
    rnaseq.ctcf0.vs.wt.around.boundaries[,gene.type:="non-DE"]
    rnaseq.ctcf0.vs.wt.around.boundaries[log2FoldChange>=rnaseq.threshold.log2fc&padj<=rnaseq.threshold.padj,gene.type:="DE (WT<CTCF0)"]
    rnaseq.ctcf0.vs.wt.around.boundaries[log2FoldChange<=rnaseq.threshold.log2fc&padj<=rnaseq.threshold.padj,gene.type:="DE (WT>CTCF0)"]
    ##insulation score around selected CD boundaries
    insulation.score.around.boundaries=rbindlist(lapply(1:nrow(CD.boundaries.selected),function(i){
        chr1=CD.boundaries.selected[i,chr]
        position1=CD.boundaries.selected[i,position]
        tmp=insulation.scores[chr==chr1&position>=position1-plot.range-binsize&position<=position1+plot.range+binsize] 
        tmp[,distance:=position-position1]
        if(nrow(tmp)==0)
            return(cbind(i,data.table(distance=NaN,score.diff=NaN)))
        cbind(i,tmp[,.(distance,score.diff=score.CTCF0-score.wildtype)])
    }))

    if(CD.boundaries.type=="CTCF-bound")
    {
        ##flag CTCF-bound CD boundaries as TSS proximal (+/-200bp) or TSS distal
        TSS.around.boundaries=rbindlist(lapply(1:nrow(CD.boundaries.selected),function(i){
            chr1=CD.boundaries.selected[i,chr]
            position1=CD.boundaries.selected[i,position]
            tmp=rnaseq.ctcf0.vs.wt[Chromosome==chr1&TSS>=position1-200&TSS<=position1+200] 
            cbind(i,CD.boundaries.selected[i,.(score.diff)],tmp[,.(TSS.proximal=(.N>0))])
        }))
    }


    ##prepare y axis labels (delta physical insulation score)
    y.breaks=as.integer(seq(1,nrow(CD.boundaries.selected),length.out=20))
    y.labels=as.character(signif(CD.boundaries.selected[y.breaks,score.diff],3))
    ##add breaks for specific values of score.diff
    for(value in c(-0.1,0,0.1))
    {
        b=CD.boundaries.selected[,.(score.diff,i=1:.N)][,.(score.diff.1=shift(sign(score.diff-value)),i.1=shift(i),score.diff.2=sign(score.diff-value),i.2=i)][is.finite(i.1)&score.diff.1*score.diff.2<=0,mean(c(i.1,i.2))]
        y.breaks=c(y.breaks,b)
        y.labels=c(y.labels,"")
    }

    chipseq.logFC.limit=quantile(c(ctcf.wt.vs.ctcf0.around.boundaries[,abs(best.logFC)],
                                   cp190.wt.vs.cp190KO.around.boundaries[,abs(best.logFC)],
                                   cp190.ctcf0.vs.cp190KO.around.boundaries[,abs(best.logFC)],
                                   cp190.wt.vs.ctcf0.around.boundaries[,abs(best.logFC)]),
                                 probs=0.99,na.rm=TRUE)



    p1=ggplot(aes(x=distance/1000,xend=distance/1000,y=i-0.5,yend=i+0.5,color=best.logFC),data=ctcf.wt.vs.ctcf0.around.boundaries)
    p1=p1+geom_segment(size=1)
    p1=p1+xlab("Distance to\nCD boundary (kb)")
    p1=p1+ylab(paste0(CD.boundaries.type," CD boundaries (n=",nrow(CD.boundaries.selected)," ranked by delta physical insulation score CTCF0-WT)"))
    p1=p1+scale_x_continuous(expand=c(0,0),limits=c(-plot.range/1000,plot.range/1000))
    p1=p1+scale_y_continuous(expand=c(0,0),breaks=y.breaks,labels=y.labels)
    p1=p1+scale_color_gradientn("CTCF\noccupancy\nin WT\n(log2FC)",colours=brewer.pal(9,"Greens"),limits=c(0,chipseq.logFC.limit),oob=scales::squish)
    p1=p1+theme_bw()
    p1=p1+theme(panel.grid=element_blank(),legend.position="top",legend.title=element_text(size=8),legend.key.width=unit(0.1,"inches"),legend.text=element_text(size=8,angle=90,vjust=0.5))

    p2=ggplot(aes(x=distance/1000,xend=distance/1000,y=i-0.5,yend=i+0.5,color=best.logFC),data=cp190.wt.vs.cp190KO.around.boundaries)
    p2=p2+geom_segment(size=1)
    p2=p2+xlab("Distance to\nCD boundary (kb)")
    p2=p2+ylab("")
    p2=p2+scale_x_continuous(expand=c(0,0),limits=c(-plot.range/1000,plot.range/1000))
    p2=p2+scale_y_continuous(expand=c(0,0))
    p2=p2+scale_color_gradientn("Cp190\noccupancy\nin WT\n(log2FC)",colours=brewer.pal(9,"Greens"),limits=c(0,chipseq.logFC.limit),oob=scales::squish)
    p2=p2+theme_bw()
    p2=p2+theme(panel.grid=element_blank(),legend.position="top",legend.title=element_text(size=8),axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.key.width=unit(0.1,"inches"),legend.text=element_text(size=8,angle=90,vjust=0.5))

    p3=ggplot(aes(x=distance/1000,xend=distance/1000,y=i-0.5,yend=i+0.5,color=best.logFC),data=cp190.ctcf0.vs.cp190KO.around.boundaries)
    p3=p3+geom_segment(size=1)
    p3=p3+xlab("Distance to\nCD boundary (kb)")
    p3=p3+ylab("")
    p3=p3+scale_x_continuous(expand=c(0,0),limits=c(-plot.range/1000,plot.range/1000))
    p3=p3+scale_y_continuous(expand=c(0,0))
    p3=p3+scale_color_gradientn("Cp190\noccupancy\nin CTCF0\n(log2FC)",colours=brewer.pal(9,"Greens"),limits=c(0,chipseq.logFC.limit),oob=scales::squish)
    p3=p3+theme_bw()
    p3=p3+theme(panel.grid=element_blank(),legend.position="top",legend.title=element_text(size=8),axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.key.width=unit(0.1,"inches"),legend.text=element_text(size=8,angle=90,vjust=0.5))

    p4=ggplot(aes(x=distance/1000,xend=distance/1000,y=i-0.5,yend=i+0.5,color=best.logFC),data=cp190.wt.vs.ctcf0.around.boundaries)
    p4=p4+geom_segment(size=1)
    p4=p4+xlab("Distance to\nCD boundary (kb)")
    p4=p4+ylab("")
    p4=p4+scale_x_continuous(expand=c(0,0),limits=c(-plot.range/1000,plot.range/1000))
    p4=p4+scale_y_continuous(expand=c(0,0))
    p4=p4+scale_color_gradientn("Cp190\ndifferential\nbinding\nlog2(WT/CTCF0)",colours=brewer.pal(11,"PRGn"),limits=c(-chipseq.logFC.limit,chipseq.logFC.limit),oob=scales::squish)
    p4=p4+theme_bw()
    p4=p4+theme(panel.grid=element_blank(),legend.position="top",legend.title=element_text(size=8),axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.key.width=unit(0.1,"inches"),legend.text=element_text(size=8,angle=90,vjust=0.5))

    p5=ggplot(aes(x=distance/1000,xend=distance/1000,y=i-0.5,yend=i+0.5,color=gene.type),data=rnaseq.ctcf0.vs.wt.around.boundaries)
    p5=p5+geom_segment(size=1,alpha=0.7)
    p5=p5+xlab("Distance to\nCD boundary (kb)")
    p5=p5+ylab("")
    p5=p5+scale_x_continuous(expand=c(0,0),limits=c(-plot.range/1000,plot.range/1000))
    p5=p5+scale_y_continuous(expand=c(0,0))
    p5=p5+scale_color_manual("RNA-seq",values=c("non-DE"="grey","DE (WT>CTCF0)"="blue","DE (WT<CTCF0)"="red"))
    p5=p5+theme_bw()
    p5=p5+theme(panel.grid=element_blank(),legend.position="top",legend.title=element_text(size=8),axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.key.width=unit(0.1,"inches"),legend.direction="vertical")#,legend.text=element_text(size=8,angle=90,vjust=0.5))

    p6=ggplot(aes(x=distance/1000,y=i,fill=score.diff),data=insulation.score.around.boundaries)
    if(CD.boundaries.type=="All")
        p6=p6+geom_raster()
    if(CD.boundaries.type=="CTCF-bound")
        p6=p6+geom_tile(width=binsize/1000)
    p6=p6+xlab("Distance to\nCD boundary (kb)")
    p6=p6+ylab("")
    p6=p6+scale_x_continuous(expand=c(0,0),limits=c(-plot.range/1000,plot.range/1000))
    p6=p6+scale_y_continuous(expand=c(0,0))
    m=insulation.scores[is.finite(score.CTCF0-score.wildtype),quantile(abs(score.CTCF0-score.wildtype),probs=0.99)]
    p6=p6+scale_fill_gradientn("Delta physical\ninsulation score\n(CTCF0-WT)",colours=rev(brewer.pal(11,"RdBu")),limits=c(-m,m),oob=scales::squish)
    p6=p6+theme_bw()
    p6=p6+theme(panel.grid=element_blank(),legend.position="top",legend.title=element_text(size=8),axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.key.width=unit(0.1,"inches"),legend.text=element_text(size=8,angle=90,vjust=0.5))

    p7=ggplot(aes(x=0,y=i,fill=boundary.type),data=CD.boundaries.selected[,.(i=1:.N,boundary.type)])
    p7=p7+geom_raster()
    p7=p7+scale_fill_manual("CD boundary type",values=c("none"="white","Common"="blue","WT-specific"="red","CTCF0-specific"="green"),breaks=c("none","Common","WT-specific","CTCF0-specific"))
    p7=p7+scale_x_continuous(expand=c(0,0))
    p7=p7+scale_y_continuous(expand=c(0,0))
    p7=p7+theme_bw()
    p7=p7+xlab("")
    p7=p7+ylab("")
    p7=p7+theme(axis.text=element_blank(),axis.ticks=element_blank(),legend.position="top",legend.direction="vertical",legend.title=element_text(size=8),legend.text=element_text(size=8))

    if(CD.boundaries.type=="CTCF-bound")
    {
        p8=ggplot(aes(x=0,y=i,fill=ifelse(TSS.proximal,"TSS proximal","TSS distal")),data=TSS.around.boundaries)
        p8=p8+geom_raster()
        p8=p8+scale_fill_manual("CTCF type",values=c("TSS distal"="grey90","TSS proximal"="black"),breaks=c("TSS distal","TSS proximal"))
        p8=p8+scale_x_continuous(expand=c(0,0))
        p8=p8+scale_y_continuous(expand=c(0,0))
        p8=p8+theme_bw()
        p8=p8+xlab("")
        p8=p8+ylab("")
        p8=p8+theme(axis.text=element_blank(),axis.ticks=element_blank(),legend.position="top",legend.direction="vertical",legend.title=element_text(size=8),legend.text=element_text(size=8))
    }

    if(CD.boundaries.type=="All")
    {
        pdf("output/figure_5c.pdf",15,12)
        print(plot_grid(p1,p2,p3,p4,p5,p6,p7,nrow=1,align="h",axis="tb",rel_widths=c(4,3,3,3,3,3,2)))
    }
    if(CD.boundaries.type=="CTCF-bound")
    {        
        pdf("output/figure_5d.pdf",16.4,12)
        print(plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,nrow=1,align="h",axis="tb",rel_widths=c(4,3,3,3,3,3,2,2)))
    }
    dev.off()

    ##save data
    if(CD.boundaries.type=="All")
    {
        filename="output/figure_5c.txt"
        tmp=rbind(
            ctcf.wt.vs.ctcf0.around.boundaries[,.(distance=format(distance,scientific=FALSE,trim=TRUE),CD.boundary.rank=i,CD.boundary.delta.physical.insulation.score=CD.boundaries.selected[i,score.diff],CD.boundary.type=CD.boundaries.selected[i,boundary.type],variable="CTCF occupancy in WT",value=best.logFC)],
            cp190.wt.vs.cp190KO.around.boundaries[,.(distance=format(distance,scientific=FALSE,trim=TRUE),CD.boundary.rank=i,CD.boundary.delta.physical.insulation.score=CD.boundaries.selected[i,score.diff],CD.boundary.type=CD.boundaries.selected[i,boundary.type],variable="Cp190 occupancy in WT",value=best.logFC)],
            cp190.ctcf0.vs.cp190KO.around.boundaries[,.(distance=format(distance,scientific=FALSE,trim=TRUE),CD.boundary.rank=i,CD.boundary.delta.physical.insulation.score=CD.boundaries.selected[i,score.diff],CD.boundary.type=CD.boundaries.selected[i,boundary.type],variable="Cp190 occupancy in CTCF0",value=best.logFC)],
            cp190.wt.vs.ctcf0.around.boundaries[,.(distance=format(distance,scientific=FALSE,trim=TRUE),CD.boundary.rank=i,CD.boundary.delta.physical.insulation.score=CD.boundaries.selected[i,score.diff],CD.boundary.type=CD.boundaries.selected[i,boundary.type],variable="Cp190 differential binding log2(WT/CTCF0)",value=best.logFC)],
            rnaseq.ctcf0.vs.wt.around.boundaries[,.(distance=format(distance,scientific=FALSE,trim=TRUE),CD.boundary.rank=i,CD.boundary.delta.physical.insulation.score=CD.boundaries.selected[i,score.diff],CD.boundary.type=CD.boundaries.selected[i,boundary.type],variable="RNAseq",value=gene.type)],
            insulation.score.around.boundaries[,.(distance=format(distance,scientific=FALSE,trim=TRUE),CD.boundary.rank=i,CD.boundary.delta.physical.insulation.score=CD.boundaries.selected[i,score.diff],CD.boundary.type=CD.boundaries.selected[i,boundary.type],variable="Delta physical insulation score",value=score.diff)])
    }
    if(CD.boundaries.type=="CTCF-bound")
    {
        filename="output/figure_5d.txt"
        tmp=rbind(
            ctcf.wt.vs.ctcf0.around.boundaries[,.(distance=format(distance,scientific=FALSE,trim=TRUE),CD.boundary.rank=i,CD.boundary.delta.physical.insulation.score=CD.boundaries.selected[i,score.diff],CD.boundary.type=CD.boundaries.selected[i,boundary.type],CTCF.type=TSS.around.boundaries[i,ifelse(TSS.proximal,"TSS proximal","TSS distal")],variable="CTCF occupancy in WT",value=best.logFC)],
            cp190.wt.vs.cp190KO.around.boundaries[,.(distance=format(distance,scientific=FALSE,trim=TRUE),CD.boundary.rank=i,CD.boundary.delta.physical.insulation.score=CD.boundaries.selected[i,score.diff],CD.boundary.type=CD.boundaries.selected[i,boundary.type],CTCF.type=TSS.around.boundaries[i,ifelse(TSS.proximal,"TSS proximal","TSS distal")],variable="Cp190 occupancy in WT",value=best.logFC)],
            cp190.ctcf0.vs.cp190KO.around.boundaries[,.(distance=format(distance,scientific=FALSE,trim=TRUE),CD.boundary.rank=i,CD.boundary.delta.physical.insulation.score=CD.boundaries.selected[i,score.diff],CD.boundary.type=CD.boundaries.selected[i,boundary.type],CTCF.type=TSS.around.boundaries[i,ifelse(TSS.proximal,"TSS proximal","TSS distal")],variable="Cp190 occupancy in CTCF0",value=best.logFC)],
            cp190.wt.vs.ctcf0.around.boundaries[,.(distance=format(distance,scientific=FALSE,trim=TRUE),CD.boundary.rank=i,CD.boundary.delta.physical.insulation.score=CD.boundaries.selected[i,score.diff],CD.boundary.type=CD.boundaries.selected[i,boundary.type],CTCF.type=TSS.around.boundaries[i,ifelse(TSS.proximal,"TSS proximal","TSS distal")],variable="Cp190 differential binding log2(WT/CTCF0)",value=best.logFC)],
            rnaseq.ctcf0.vs.wt.around.boundaries[,.(distance=format(distance,scientific=FALSE,trim=TRUE),CD.boundary.rank=i,CD.boundary.delta.physical.insulation.score=CD.boundaries.selected[i,score.diff],CD.boundary.type=CD.boundaries.selected[i,boundary.type],CTCF.type=TSS.around.boundaries[i,ifelse(TSS.proximal,"TSS proximal","TSS distal")],variable="RNAseq",value=gene.type)],
            insulation.score.around.boundaries[,.(distance=format(distance,scientific=FALSE,trim=TRUE),CD.boundary.rank=i,CD.boundary.delta.physical.insulation.score=CD.boundaries.selected[i,score.diff],CD.boundary.type=CD.boundaries.selected[i,boundary.type],CTCF.type=TSS.around.boundaries[i,ifelse(TSS.proximal,"TSS proximal","TSS distal")],variable="Delta physical insulation score",value=score.diff)])
    }
    fwrite(tmp,file=filename,sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE,na="NA")
    
}



#############################################
## wildtype CD boundaries around Cp190 peaks in WT 
#############################################

## wildtype CD boundaries around Cp190 peaks in WT
boundaries.around.cp190=rbindlist(lapply(1:nrow(chipseq.cp190.wt.vs.cp190KO),function(i){
    chr1=chipseq.cp190.wt.vs.cp190KO[i,chr]
    position1=chipseq.cp190.wt.vs.cp190KO[i,best.pos]
    ##consider only wildtype boundaries (exist.witdtype==TRUE)
    tmp=CD.boundaries[exist.wildtype==TRUE&chr==chr1&position>=position1-plot.range-binsize&position<=position1+plot.range+binsize] 
    tmp[,distance.binned:=round((position-position1)/binsize)*binsize]
    if(nrow(tmp)==0)
        return(cbind(i,data.table(distance.binned=NaN)))
    cbind(i,tmp[,.(distance.binned)])
}))
nb.cp190=nrow(chipseq.cp190.wt.vs.cp190KO)
##fraction of CP190 with AT LEAST one boundary at given distance (=> unique())
boundaries.around.cp190=unique(boundaries.around.cp190)[,.(fraction=.N/nb.cp190),by=distance.binned]
##add missing distances
distance.binned.all=data.table(distance.binned=seq(-ceiling(plot.range/binsize)*binsize,plot.range,by=binsize))
boundaries.around.cp190=boundaries.around.cp190[distance.binned.all,on="distance.binned"]
boundaries.around.cp190[is.na(fraction),fraction:=0]


## wildtype CD boundaries around random positions (same number as cp190 peaks)
chipseq.cp190.wt.vs.cp190KO.random=get_random_positions(nrow(chipseq.cp190.wt.vs.cp190KO))
boundaries.around.cp190.random=rbindlist(lapply(1:nrow(chipseq.cp190.wt.vs.cp190KO.random),function(i){
    chr1=chipseq.cp190.wt.vs.cp190KO.random[i,chr]
    position1=chipseq.cp190.wt.vs.cp190KO.random[i,position]
    ##consider only wildtype boundaries (exist.witdtype==TRUE)
    tmp=CD.boundaries[exist.wildtype==TRUE&chr==chr1&position>=position1-plot.range-binsize&position<=position1+plot.range+binsize] 
    tmp[,distance.binned:=round((position-position1)/binsize)*binsize]
    if(nrow(tmp)==0)
        return(cbind(i,data.table(distance.binned=NaN)))
    cbind(i,tmp[,.(distance.binned)])
}))
##fraction of random positions with AT LEAST one boundary at given distance (=> unique())
boundaries.around.cp190.random=unique(boundaries.around.cp190.random)[,.(fraction=.N/nb.cp190),by=distance.binned]
##add missing distances
distance.binned.all=data.table(distance.binned=seq(-ceiling(plot.range/binsize)*binsize,plot.range,by=binsize))
boundaries.around.cp190.random=boundaries.around.cp190.random[distance.binned.all,on="distance.binned"]
boundaries.around.cp190.random[is.na(fraction),fraction:=0]


##plot
pdf("output/figure_S5d.pdf",10,5)
boundaries.around.cp190[,label:="Cp190 peaks"]
boundaries.around.cp190.random[,label:="random regions"]
fraction.at.0=boundaries.around.cp190[distance.binned==0,fraction]
fraction.random.average=boundaries.around.cp190.random[,mean(fraction)]
annotation=paste0(round(fraction.at.0*100),"% (",round(fraction.at.0/fraction.random.average),"x enrichment)")
p=ggplot(aes(x=distance.binned/1000,y=fraction,alpha=label),data=rbind(boundaries.around.cp190,boundaries.around.cp190.random))
p=p+geom_line()
p=p+scale_alpha_manual("",values=c("random regions"=0.3,"Cp190 peaks"=1))
p=p+scale_x_continuous("Distance to Cp190 peak (kb)",limits=c(-plot.range/1000,plot.range/1000))
p=p+scale_y_continuous(paste0("% of Cp190 peaks with CD boundary\nat a given distance (per ",binsize/1000,"kb)"),labels = scales::percent)
p=p+ggtitle("Enrichment of CD boundaries around Cp190 peaks")
p=p+annotate("text",x=0,y=boundaries.around.cp190[,max(fraction,na.rm=TRUE)],label=annotation,hjust=-0.1,vjust=1)
p=p+theme_bw()
print(p)
dev.off()

##save data
fwrite(rbind(boundaries.around.cp190,boundaries.around.cp190.random)[,.(type=label,distance=format(distance.binned,scientific=FALSE,trim=TRUE),fraction)],file="output/figure_S5d.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE,na="NA")

#############################################
## Cp190 peaks in WT around wildtype CD boundaries
#############################################

## Cp190 peaks in WT around wildtype CD boundaries (exist.wildtype==TRUE)
cp190.around.boundaries=rbindlist(lapply(1:nrow(CD.boundaries[exist.wildtype==TRUE]),function(i){
    chr1=CD.boundaries[exist.wildtype==TRUE][i,chr]
    position1=CD.boundaries[exist.wildtype==TRUE][i,position]
    tmp=chipseq.cp190.wt.vs.cp190KO[chr==chr1&best.pos>=position1-plot.range-binsize&best.pos<=position1+plot.range+binsize] 
    tmp[,distance.binned:=round((best.pos-position1)/binsize)*binsize]
    if(nrow(tmp)==0)
        return(cbind(i,data.table(distance.binned=NaN)))
    cbind(i,tmp[,.(distance.binned)])
}))
nb.boundaries=nrow(CD.boundaries[exist.wildtype==TRUE])
##fraction of boundaries with AT LEAST one CP190 at given distance (=> unique())
cp190.around.boundaries=unique(cp190.around.boundaries)[,.(fraction=.N/nb.boundaries),by=distance.binned]
##add missing distances
distance.binned.all=data.table(distance.binned=seq(-ceiling(plot.range/binsize)*binsize,plot.range,by=binsize))
cp190.around.boundaries=cp190.around.boundaries[distance.binned.all,on="distance.binned"]
cp190.around.boundaries[is.na(fraction),fraction:=0]


## CP190 around random positions (same number as wildtype CD boundaires)
CD.boundaries.random=get_random_positions(nrow(CD.boundaries[exist.wildtype==TRUE]))
cp190.around.boundaries.random=rbindlist(lapply(1:nrow(CD.boundaries.random),function(i){
    chr1=CD.boundaries.random[i,chr]
    position1=CD.boundaries.random[i,position]
    tmp=chipseq.cp190.wt.vs.cp190KO[chr==chr1&best.pos>=position1-plot.range-binsize&best.pos<=position1+plot.range+binsize] 
    tmp[,distance.binned:=round((best.pos-position1)/binsize)*binsize]
    if(nrow(tmp)==0)
        return(cbind(i,data.table(distance.binned=NaN)))
    cbind(i,tmp[,.(distance.binned)])
}))
##fraction of boundaries with AT LEAST one random position at given distance (=> unique())
cp190.around.boundaries.random=unique(cp190.around.boundaries.random)[,.(fraction=.N/nb.boundaries),by=distance.binned]
##add missing distances
cp190.around.boundaries.random=cp190.around.boundaries.random[distance.binned.all,on="distance.binned"]
cp190.around.boundaries.random[is.na(fraction),fraction:=0]


##plot
pdf("output/figure_S5e.pdf",10,5)
cp190.around.boundaries[,label:="CD boundaries"]
cp190.around.boundaries.random[,label:="random regions"]
fraction.at.0=cp190.around.boundaries[distance.binned==0,fraction]
fraction.random.average=cp190.around.boundaries.random[,mean(fraction)]
annotation=paste0(round(fraction.at.0*100),"% (",round(fraction.at.0/fraction.random.average),"x enrichment)")
p=ggplot(aes(x=distance.binned/1000,y=fraction,alpha=label),data=rbind(cp190.around.boundaries,cp190.around.boundaries.random))
p=p+geom_line()
p=p+scale_alpha_manual("",values=c("random regions"=0.3,"CD boundaries"=1))
p=p+scale_x_continuous("Distance to CD boundary (kb)",limits=c(-plot.range/1000,plot.range/1000))
p=p+scale_y_continuous(paste0("% of CD boundaries with Cp190 peak\nat a given distance (per ",binsize/1000,"kb)"),labels = scales::percent)
p=p+ggtitle("Enrichment of Cp190 peaks around CD boundaries")
p=p+annotate("text",x=0,y=cp190.around.boundaries[,max(fraction,na.rm=TRUE)],label=annotation,hjust=-0.1,vjust=1)
p=p+theme_bw()
print(p)
dev.off()

##save data
fwrite(rbind(cp190.around.boundaries,cp190.around.boundaries.random)[,.(type=label,distance=format(distance.binned,scientific=FALSE,trim=TRUE),fraction)],file="output/figure_S5e.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE,na="NA")



#############################################
## Residual CD boundaries in CTCF0
## Residual Cp190 peaks in CTCF0
## and TSS proximity to CTCF bound CD boundaries
#############################################

##consider only CTCF bound CD boundaries
CTCF.bound.CD.boundaries=copy(CD.boundaries[nb.neighboring.CTCF>0,.(chr,CTCF.start,CTCF.end,boundary.type,score.diff)])

##check if CTCF overlap  cp190 peaks in CTCF0
CTCF.bound.CD.boundaries=rbindlist(lapply(1:nrow(CTCF.bound.CD.boundaries),function(i){
    chr1=CTCF.bound.CD.boundaries[i,chr]
    start1=CTCF.bound.CD.boundaries[i,CTCF.start]
    end1=CTCF.bound.CD.boundaries[i,CTCF.end]
    ##Cp190 peaks in CTCF0 overlapping CTCF
    tmp=chipseq.cp190.ctcf0.vs.cp190KO[chr==chr1&start<=end1&end>=start1]
    cbind(i,CTCF.bound.CD.boundaries[i],tmp[,.(overlap.cp190.ctcf0.vs.cp190KO=(.N>0))])
}))

##Note:
## Residual CD boundaries in CTCF0 = Non WT-specific CD boundaries
## Residual Cp190 peaks in CTCF0 = Cp190 peaks in CTCF0 (overlapping CTCF)
contingency.table=CTCF.bound.CD.boundaries[,.(.N),by=.(residual.cp190.peak=overlap.cp190.ctcf0.vs.cp190KO,residual.CD.boundary=(boundary.type!="WT-specific"))]
fisher=fisher.test(reshape2::acast(contingency.table,residual.CD.boundary~residual.cp190.peak,value.var="N"))

pdf("output/figure_5e.pdf",5,5)
p=ggplot(aes(x=ifelse(residual.cp190.peak,"yes","no"),y=ifelse(residual.CD.boundary,"yes","no"),fill=N),data=contingency.table)
p=p+geom_tile(color="black")
p=p+geom_text(aes(label=N))
p=p+xlab("Residual Cp190 peak in CTCF0")
p=p+ylab("Residual CD boundary in CTCF0")
p=p+scale_x_discrete(expand=c(0,0))
p=p+scale_y_discrete(expand=c(0,0))
p=p+scale_fill_gradientn(colors=c("white","grey"))
p=p+ggtitle(paste0("CTCF-bound CD boundaries (n=",nrow(CTCF.bound.CD.boundaries),")\n","p-value=",signif(fisher$p.value,2)," (",fisher$method,", ",fisher$alternative,")"))
p=p+theme_bw()+theme(plot.title=element_text(size=9),legend.position="none")
print(p)
dev.off()

##save data
fwrite(contingency.table,file="output/figure_5e.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE,na="NA")



##Delta physical insulation score distribution for CTCF-bound CD boundary with/without residual Cp190 peak in CTCF0
n_fun=function(x){
    return(data.frame(y=median(x),label=paste0("n=",length(x))))
}
wtest1=wilcox.test(CTCF.bound.CD.boundaries[overlap.cp190.ctcf0.vs.cp190KO==FALSE,score.diff],
                   CTCF.bound.CD.boundaries[overlap.cp190.ctcf0.vs.cp190KO==TRUE,score.diff])
p1=ggplot(aes(x=factor(ifelse(residual.cp190.peak,"residual peak","no peak"),levels=c("residual peak","no peak")),y=score.diff),data=CTCF.bound.CD.boundaries[,.(residual.cp190.peak=overlap.cp190.ctcf0.vs.cp190KO,score.diff)])
p1=p1+geom_boxplot()
p1=p1+stat_summary(,fun.data=n_fun,geom="text",size=3,vjust = -0.5,color="grey50")
p1=p1+xlab("Cp190 peak in CTCF0")
p1=p1+ylab("Delta physical insulation score (CTCF0-wildtype)")
p1=p1+ggtitle(paste0("CTCF-bound CD boundaries (n=",nrow(CTCF.bound.CD.boundaries),")\n","p-value=",signif(wtest1$p.value,2),"\n(",wtest1$method,", ",wtest1$alternative,", W-statistic=",signif(wtest1$statistic,5),")"))
p1=p1+theme_bw()+theme(plot.title=element_text(size=8),legend.position="none")

##Delta physical insulation score distribution for CTCF-bound CD boundary proximal/distal to TSS
wtest2=wilcox.test(TSS.around.boundaries[TSS.proximal==TRUE,score.diff],
             TSS.around.boundaries[TSS.proximal==FALSE,score.diff])
p2=ggplot(aes(x=factor(ifelse(TSS.proximal,"proximal","distal"),levels=c("proximal","distal")),y=score.diff),data=TSS.around.boundaries)
p2=p2+geom_boxplot()
p2=p2+stat_summary(,fun.data=n_fun,geom="text",size=3,vjust = -0.5,color="grey50")
p2=p2+xlab("Proximity to TSS")
p2=p2+ylab("Delta physical insulation score (CTCF0-wildtype)")
p2=p2+ggtitle(paste0("CTCF-bound CD boundaries (n=",nrow(TSS.around.boundaries),")\n","p-value=",signif(wtest2$p.value,2),"\n(",wtest2$method,", ",wtest2$alternative,", W-statistic=",signif(wtest2$statistic,5),")"))
p2=p2+theme_bw()+theme(plot.title=element_text(size=8),legend.position="none")


pdf("output/figure_5f.pdf",10,7)
print(plot_grid(p1,p2,nrow=1,align="h",axis="tb"))
dev.off()

##save data
tmp=merge(TSS.around.boundaries[,.(i,proximity.to.TSS=ifelse(TSS.proximal,"proximal","distal"))],CTCF.bound.CD.boundaries[,.(i,residual.cp190.peak=overlap.cp190.ctcf0.vs.cp190KO,delta.physical.insulation.score=score.diff)],by="i")[,.(proximity.to.TSS,residual.cp190.peak,delta.physical.insulation.score)]
fwrite(tmp[order(delta.physical.insulation.score)],file="output/figure_5f.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE,na="NA")
