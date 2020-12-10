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
library(RColorBrewer)
library(data.table)
library(readxl)
set.seed(83712)

#############################################
##paths
#############################################
path.CD.boundaries="data/SupplementaryData3_HiC_CD_boundaries.xlsx"
path.insulation.scores="data/SupplementaryData2_HiC_insulation_scores.xlsx"
path.chipseq.ctcf.wt.vs.ctcf0="data/SupplementaryData1_ChIPseq_CTCF_WT-CTCF0.xlsx"
dir.create("output/", showWarnings = FALSE,recursive=TRUE)

#############################################
##variables
#############################################
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

#############################################
##read data
#############################################
CD.boundaries=as.data.table(read_xlsx(path.CD.boundaries))
insulation.scores=as.data.table(read_xlsx(path.insulation.scores))

##CTCF peaks in WT
chipseq.ctcf.wt.vs.ctcf0=as.data.table(read_xlsx(path.chipseq.ctcf.wt.vs.ctcf0))
##keep only "up" peaks with positive fold change (filter out false positive)
chipseq.ctcf.wt.vs.ctcf0=chipseq.ctcf.wt.vs.ctcf0[direction=="up"&best.logFC>0]
##keep only selected chrs
chipseq.ctcf.wt.vs.ctcf0=chipseq.ctcf.wt.vs.ctcf0[chr%in%chrs]

#############################################
## wildtype CD boundaries around CTCF
#############################################

## wildtype CD boundaries around CTCF 
boundaries.around.ctcf=rbindlist(lapply(1:nrow(chipseq.ctcf.wt.vs.ctcf0),function(i){
    chr1=chipseq.ctcf.wt.vs.ctcf0[i,chr]
    position1=chipseq.ctcf.wt.vs.ctcf0[i,best.pos]
    ##consider only wildtype boundaries (exist.witdtype==TRUE)
    tmp=CD.boundaries[exist.wildtype==TRUE&chr==chr1&position>=position1-plot.range-binsize&position<=position1+plot.range+binsize] 
    tmp[,distance.binned:=round((position-position1)/binsize)*binsize]
    if(nrow(tmp)==0)
        return(cbind(i,data.table(distance.binned=NaN)))
    cbind(i,tmp[,.(distance.binned)])
}))
nb.ctcf=nrow(chipseq.ctcf.wt.vs.ctcf0)
##fraction of CTCF with AT LEAST one boundary at given distance (=> unique())
boundaries.around.ctcf=unique(boundaries.around.ctcf)[,.(fraction=.N/nb.ctcf),by=distance.binned]
##add missing distances
distance.binned.all=data.table(distance.binned=seq(-ceiling(plot.range/binsize)*binsize,plot.range,by=binsize))
boundaries.around.ctcf=boundaries.around.ctcf[distance.binned.all,on="distance.binned"]
boundaries.around.ctcf[is.na(fraction),fraction:=0]


## wildtype CD boundaries around random positions (same number as ctcf peaks)
chipseq.ctcf.wt.vs.ctcf0.random=get_random_positions(nrow(chipseq.ctcf.wt.vs.ctcf0))
boundaries.around.ctcf.random=rbindlist(lapply(1:nrow(chipseq.ctcf.wt.vs.ctcf0.random),function(i){
    chr1=chipseq.ctcf.wt.vs.ctcf0.random[i,chr]
    position1=chipseq.ctcf.wt.vs.ctcf0.random[i,position]
    ##consider only wildtype boundaries (exist.witdtype==TRUE)
    tmp=CD.boundaries[exist.wildtype==TRUE&chr==chr1&position>=position1-plot.range-binsize&position<=position1+plot.range+binsize] 
    tmp[,distance.binned:=round((position-position1)/binsize)*binsize]
    if(nrow(tmp)==0)
        return(cbind(i,data.table(distance.binned=NaN)))
    cbind(i,tmp[,.(distance.binned)])
}))
##fraction of random positions with AT LEAST one boundary at given distance (=> unique())
boundaries.around.ctcf.random=unique(boundaries.around.ctcf.random)[,.(fraction=.N/nb.ctcf),by=distance.binned]
##add missing distances
distance.binned.all=data.table(distance.binned=seq(-ceiling(plot.range/binsize)*binsize,plot.range,by=binsize))
boundaries.around.ctcf.random=boundaries.around.ctcf.random[distance.binned.all,on="distance.binned"]
boundaries.around.ctcf.random[is.na(fraction),fraction:=0]


##plot
pdf("output/figure_2a.pdf",10,5)
boundaries.around.ctcf[,label:="CTCF peaks"]
boundaries.around.ctcf.random[,label:="random regions"]
fraction.at.0=boundaries.around.ctcf[distance.binned==0,fraction]
fraction.random.average=boundaries.around.ctcf.random[,mean(fraction)]
annotation=paste0(round(fraction.at.0*100),"% (",round(fraction.at.0/fraction.random.average),"x enrichment)")
p=ggplot(aes(x=distance.binned/1000,y=fraction,alpha=label),data=rbind(boundaries.around.ctcf,boundaries.around.ctcf.random))
p=p+geom_line()
p=p+scale_alpha_manual("",values=c("random regions"=0.3,"CTCF peaks"=1))
p=p+scale_x_continuous("Distance to CTCF peak (kb)",limits=c(-plot.range/1000,plot.range/1000))
p=p+scale_y_continuous(paste0("% of CTCF peaks with CD boundary\nat a given distance (per ",binsize/1000,"kb)"),labels = scales::percent)
p=p+ggtitle("Enrichment of CD boundaries around CTCF peaks")
p=p+annotate("text",x=0,y=boundaries.around.ctcf[,max(fraction,na.rm=TRUE)],label=annotation,hjust=-0.1,vjust=1)
p=p+theme_bw()
print(p)
dev.off()
##save data
fwrite(rbind(boundaries.around.ctcf,boundaries.around.ctcf.random)[,.(type=label,distance=format(distance.binned,scientific=FALSE,trim=TRUE),fraction)],file="output/figure_2a.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE,na="NA")


#############################################
## CTCF around wildtype CD boundaries
#############################################

## CTCF around wildtype CD boundaries (exist.wildtype==TRUE)
ctcf.around.boundaries=rbindlist(lapply(1:nrow(CD.boundaries[exist.wildtype==TRUE]),function(i){
    chr1=CD.boundaries[exist.wildtype==TRUE][i,chr]
    position1=CD.boundaries[exist.wildtype==TRUE][i,position]
    tmp=chipseq.ctcf.wt.vs.ctcf0[chr==chr1&best.pos>=position1-plot.range-binsize&best.pos<=position1+plot.range+binsize] 
    tmp[,distance.binned:=round((best.pos-position1)/binsize)*binsize]
    if(nrow(tmp)==0)
        return(cbind(i,data.table(distance.binned=NaN)))
    cbind(i,tmp[,.(distance.binned)])
}))
nb.boundaries=nrow(CD.boundaries[exist.wildtype==TRUE])
##fraction of boundaries with AT LEAST one CTCF at given distance (=> unique())
ctcf.around.boundaries=unique(ctcf.around.boundaries)[,.(fraction=.N/nb.boundaries),by=distance.binned]
##add missing distances
distance.binned.all=data.table(distance.binned=seq(-ceiling(plot.range/binsize)*binsize,plot.range,by=binsize))
ctcf.around.boundaries=ctcf.around.boundaries[distance.binned.all,on="distance.binned"]
ctcf.around.boundaries[is.na(fraction),fraction:=0]


## CTCF around random positions (same number as wildtype CD boundaires)
CD.boundaries.random=get_random_positions(nrow(CD.boundaries[exist.wildtype==TRUE]))
ctcf.around.boundaries.random=rbindlist(lapply(1:nrow(CD.boundaries.random),function(i){
    chr1=CD.boundaries.random[i,chr]
    position1=CD.boundaries.random[i,position]
    tmp=chipseq.ctcf.wt.vs.ctcf0[chr==chr1&best.pos>=position1-plot.range-binsize&best.pos<=position1+plot.range+binsize] 
    tmp[,distance.binned:=round((best.pos-position1)/binsize)*binsize]
    if(nrow(tmp)==0)
        return(cbind(i,data.table(distance.binned=NaN)))
    cbind(i,tmp[,.(distance.binned)])
}))
##fraction of boundaries with AT LEAST one random position at given distance (=> unique())
ctcf.around.boundaries.random=unique(ctcf.around.boundaries.random)[,.(fraction=.N/nb.boundaries),by=distance.binned]
##add missing distances
ctcf.around.boundaries.random=ctcf.around.boundaries.random[distance.binned.all,on="distance.binned"]
ctcf.around.boundaries.random[is.na(fraction),fraction:=0]


##plot
pdf("output/figure_2b.pdf",10,5)
ctcf.around.boundaries[,label:="CD boundaries"]
ctcf.around.boundaries.random[,label:="random regions"]
fraction.at.0=ctcf.around.boundaries[distance.binned==0,fraction]
fraction.random.average=ctcf.around.boundaries.random[,mean(fraction)]
annotation=paste0(round(fraction.at.0*100),"% (",round(fraction.at.0/fraction.random.average),"x enrichment)")
p=ggplot(aes(x=distance.binned/1000,y=fraction,alpha=label),data=rbind(ctcf.around.boundaries,ctcf.around.boundaries.random))
p=p+geom_line()
p=p+scale_alpha_manual("",values=c("random regions"=0.3,"CD boundaries"=1))
p=p+scale_x_continuous("Distance to CD boundary (kb)",limits=c(-plot.range/1000,plot.range/1000))
p=p+scale_y_continuous(paste0("% of CD boundaries with CTCF peak\nat a given distance (per ",binsize/1000,"kb)"),labels = scales::percent)
p=p+ggtitle("Enrichment of CTCF peaks around CD boundaries")
p=p+annotate("text",x=0,y=ctcf.around.boundaries[,max(fraction,na.rm=TRUE)],label=annotation,hjust=-0.1,vjust=1)
p=p+theme_bw()
print(p)
dev.off()
##save data
fwrite(rbind(ctcf.around.boundaries,ctcf.around.boundaries.random)[,.(type=label,distance=format(distance.binned,scientific=FALSE,trim=TRUE),fraction)],file="output/figure_2b.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE,na="NA")


#############################################
## CTCF around all CD boundaries
#############################################

##order CD boundaries by delta physical insulation score
CD.boundaries[,score.diff:=score.CTCF0-score.wildtype]
setorder(CD.boundaries,score.diff)
CD.boundaries[,i:=1:.N]
setkey(CD.boundaries,i)

##add boundary.type
CD.boundaries[exist.wildtype==TRUE&exist.CTCF0==TRUE,boundary.type:="Common"]
CD.boundaries[exist.wildtype==TRUE&exist.CTCF0==FALSE,boundary.type:="WT-specific"]
CD.boundaries[exist.wildtype==FALSE&exist.CTCF0==TRUE,boundary.type:="CTCF0-specific"]

## CTCF around all CD boundaries 
ctcf.around.boundaries=rbindlist(lapply(1:nrow(CD.boundaries),function(i){
    chr1=CD.boundaries[i,chr]
    position1=CD.boundaries[i,position]
    tmp=chipseq.ctcf.wt.vs.ctcf0[chr==chr1&best.pos>=position1-plot.range-binsize&best.pos<=position1+plot.range+binsize] 
    tmp[,distance:=(best.pos-position1)]
    if(nrow(tmp)==0)
        return(cbind(i,data.table(distance=NaN)))
    cbind(i,tmp[,.(distance)])
}))

##prepare y axis labels (delta physical insulation score)
y.breaks=as.integer(seq(1,nrow(CD.boundaries),length.out=20))
y.labels=as.character(signif(CD.boundaries[y.breaks,score.diff],3))
##add breaks for specific values of score.diff
for(value in c(-0.1,0,0.1))
{
    b=CD.boundaries[,.(score.diff,i=1:.N)][,.(score.diff.1=shift(sign(score.diff-value)),i.1=shift(i),score.diff.2=sign(score.diff-value),i.2=i)][is.finite(i.1)&score.diff.1*score.diff.2<=0,mean(c(i.1,i.2))]
    y.breaks=c(y.breaks,b)
    y.labels=c(y.labels,"")
}
    
pdf("output/figure_2d.pdf",10,12)
p1=ggplot(aes(x=distance/1000,xend=distance/1000,y=i-0.5,yend=i+0.5),data=ctcf.around.boundaries)
p1=p1+geom_segment(size=1)
p1=p1+xlab("Distance to CD boundary (kb)")
p1=p1+ylab("CD boundaries (ranked by delta physical insulation score CTCF0-WT)")
p1=p1+scale_x_continuous(expand=c(0,0),limits=c(-plot.range/1000,plot.range/1000))
p1=p1+scale_y_continuous(expand=c(0,0),breaks=y.breaks,labels=y.labels)
p1=p1+theme_bw()
p1=p1+ggtitle("CTCF peaks around CD boundaries")
p1=p1+theme(panel.grid=element_blank(),legend.position="top")

p2=ggplot(aes(x=0,y=i,fill=boundary.type),data=CD.boundaries[,.(i=1:.N,boundary.type)])
p2=p2+geom_raster()
p2=p2+scale_fill_manual("CD boundary type",values=c("none"="white","Common"="blue","WT-specific"="red","CTCF0-specific"="green"),breaks=c("none","Common","WT-specific","CTCF0-specific"))
p2=p2+scale_x_continuous(expand=c(0,0))
p2=p2+scale_y_continuous(expand=c(0,0))
p2=p2+theme_bw()
p2=p2+xlab("")
p2=p2+ylab("")
p2=p2+theme(axis.text=element_blank(),axis.ticks=element_blank())

print(plot_grid(p1,p2,nrow=1,align="h",axis="tb",rel_widths=c(3,1)))
dev.off()

##save data
fwrite(ctcf.around.boundaries[,.(distance=format(distance,scientific=FALSE,trim=TRUE),CD.boundary.rank=i,CD.boundary.delta.physical.insulation.score=CD.boundaries[i,score.diff],CD.boundary.type=CD.boundaries[i,boundary.type])],file="output/figure_2d.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE,na="NA")

#############################################
## Physical insulation score around CTCF
#############################################
##order CTCF by best.logFC
setorder(chipseq.ctcf.wt.vs.ctcf0,best.logFC)
chipseq.ctcf.wt.vs.ctcf0[,i:=1:.N]
setkey(chipseq.ctcf.wt.vs.ctcf0,i)

##add delta physical insulation score
insulation.scores[,score.diff:=score.CTCF0-score.wildtype]

##insulation score around CTCF
insulation.score.around.ctcf=rbindlist(lapply(1:nrow(chipseq.ctcf.wt.vs.ctcf0),function(i){
    chr1=chipseq.ctcf.wt.vs.ctcf0[i,chr]
    position1=chipseq.ctcf.wt.vs.ctcf0[i,floor(best.pos/binsize)*binsize+binsize/2] ##center of the bin containing best.pos
    ##consider only wildtype boundaries (exist.witdtype==TRUE)
    tmp=insulation.scores[chr==chr1&position>=position1-plot.range-binsize&position<=position1+plot.range+binsize] 
    tmp[,distance:=position-position1]
    if(nrow(tmp)==0)
        return(cbind(i,data.table(distance=NaN,score.wildtype=NaN,score.CTCF0=NaN,score.diff=NaN)))
    cbind(i,tmp[,.(distance,score.wildtype,score.CTCF0,score.diff)])
}))


##prepare y axis labels (CTCF best.logFC)
y.breaks=as.integer(seq(1,nrow(chipseq.ctcf.wt.vs.ctcf0),length.out=20))
y.labels=as.character(signif(chipseq.ctcf.wt.vs.ctcf0[y.breaks,best.logFC],3))


##Delta physical insulation score around CTCF
pdf("output/figure_2e.pdf",8,10)
p=ggplot(aes(x=distance/1000,y=i,fill=score.diff),data=insulation.score.around.ctcf)
p=p+geom_raster()
m=insulation.scores[is.finite(score.diff),quantile(abs(score.diff),probs=0.99)]
p=p+scale_fill_gradientn("Delta physical\ninsulation score",colours=rev(brewer.pal(11,"RdBu")),limits=c(-m,m),oob=scales::squish)
p=p+xlab("Distance to CTCF peak (kbp)")
p=p+ylab("CTCF peak (ranked by logFC)")
p=p+scale_x_continuous(expand=c(0,0))
p=p+scale_y_continuous(expand=c(0,0),breaks=y.breaks,labels=y.labels)
p=p+theme_bw()
p=p+ggtitle("Delta physical insulation score around CTCF peaks")
print(p)
dev.off()

##save data
fwrite(insulation.score.around.ctcf[,.(distance=format(distance,scientific=FALSE,trim=TRUE),CTCF.rank=i,CTCF.logFC=chipseq.ctcf.wt.vs.ctcf0[i,best.logFC],delta.physical.insulation.score=score.diff)],file="output/figure_2e.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE,na="NA")

##Physical insulation score around CTCF
average.insulation.score.around.ctcf=rbind(
    insulation.score.around.ctcf[is.finite(score.wildtype),.(genotype="WT",score=mean(score.wildtype)),by=distance],
    insulation.score.around.ctcf[is.finite(score.CTCF0),.(genotype="CTCF0",score=mean(score.CTCF0)),by=distance]
)
pdf("output/figure_2f.pdf",10,5)
p=ggplot(aes(x=distance/1000,y=score,color=genotype),data=average.insulation.score.around.ctcf)
p=p+geom_line()
p=p+xlab("Distance to CTCF peak (kbp)")
p=p+ylab("Average physical insulation score")
p=p+scale_x_continuous(limits=c(-plot.range/1000,plot.range/1000))
p=p+scale_color_manual("",values=c("WT"="black","CTCF0"="red"))
p=p+theme_bw()
p=p+ggtitle("Physical insulation score around CTCF peaks")
print(p)
dev.off()

##save data
fwrite(average.insulation.score.around.ctcf[,.(genotype,distance=format(distance,scientific=FALSE,trim=TRUE),average.physical.insulation.score=score)],file="output/figure_2f.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE,na="NA")



##Delta physical insulation score at CTCF peaks (nearest bin boundary position)
insulation.score.at.ctcf=merge(chipseq.ctcf.wt.vs.ctcf0[,.(chr,position=round(best.pos/binsize)*binsize,best.logFC)],insulation.scores[,.(chr,position,score.diff)],by=c("chr","position"),all.x=TRUE,all.y=FALSE)
##group best.logFC for boxplots
n.boxes=4
group.size=insulation.score.at.ctcf[,(max(best.logFC)-min(best.logFC))/n.boxes]+1e-8#2
insulation.score.at.ctcf[,group.best.logFC:=floor((best.logFC-min(best.logFC))/group.size)*group.size+group.size/2+min(best.logFC)]


pdf("output/figure_S2d.pdf",6,10)
p=ggplot(aes(x=best.logFC,y=score.diff),data=insulation.score.at.ctcf)
p=p+geom_point(alpha=0.2)
p=p+geom_boxplot(aes(x=group.best.logFC,group=group.best.logFC),outlier.shape=NA,fill=NA,color="red",width=0.99*group.size) #no need to show outliers since all points are plotted
p=p+geom_text(aes(x=group.best.logFC,label=label),data=insulation.score.at.ctcf[,.(label=paste0("N=",.N),score.diff=quantile(score.diff,0.62)),by=group.best.logFC],size=3,vjust =0.5,color="red")
p=p+xlab("CTCF occupancy\nlog2(WT/CTCF0)")
p=p+ylab("Delta physical insulation score (CTCF0-WT)")
p=p+ggtitle("Insulation score defects in CTCF0 at former CTCF peaks")
p=p+theme_bw()
print(p)
dev.off()

##save data
fwrite(insulation.score.at.ctcf[,.(CTCF.logFC=best.logFC,delta.physical.insulation.score=score.diff)],file="output/figure_S2d.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE,na="NA")
