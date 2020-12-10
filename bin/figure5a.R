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


library(data.table)
library(GenomicRanges)
library(readxl)
library(eulerr)

#############################################
##paths
#############################################

path.chipseq.ctcf.wt.vs.ctcf0="data/SupplementaryData1_ChIPseq_CTCF_WT-CTCF0.xlsx"
path.chipseq.cp190.wt.vs.cp190KO="data/SupplementaryData6_ChIPseq_Cp190_WT-Cp190KO.xlsx"
path.chipseq.cp190.wt.vs.ctcf0="data/SupplementaryData7_ChIPseq_Cp190_WT-CTCF0.xlsx"

dir.create("output/", showWarnings = FALSE,recursive=TRUE)

#############################################
##variables
#############################################
chrs=c("chr2L","chr2R","chr3L","chr3R","chr4","chrX")


#############################################
##read data
#############################################

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

##Cp190 differential binding (WT-CTCF0)
chipseq.cp190.wt.vs.ctcf0=as.data.table(read_xlsx(path.chipseq.cp190.wt.vs.ctcf0))
##keep only "up" peaks with positive fold change (WT>CTCF0)
chipseq.cp190.wt.vs.ctcf0=chipseq.cp190.wt.vs.ctcf0[direction=="up"&best.logFC>0]
##keep only selected chrs
chipseq.cp190.wt.vs.ctcf0=chipseq.cp190.wt.vs.ctcf0[chr%in%chrs]

#############################################
##Overlaps between ChIP-seq peaks/DB regions
#############################################
##convert to GRange
data1=GRanges(seqnames=chipseq.ctcf.wt.vs.ctcf0$chr,ranges=IRanges(start=chipseq.ctcf.wt.vs.ctcf0$start,end=chipseq.ctcf.wt.vs.ctcf0$end))

data2=GRanges(seqnames=chipseq.cp190.wt.vs.cp190KO$chr,ranges=IRanges(start=chipseq.cp190.wt.vs.cp190KO$start,end=chipseq.cp190.wt.vs.cp190KO$end))

data3=GRanges(seqnames=chipseq.cp190.wt.vs.ctcf0$chr,ranges=IRanges(start=chipseq.cp190.wt.vs.ctcf0$start,end=chipseq.cp190.wt.vs.ctcf0$end))


hits12=findOverlaps(data1,data2)
hits23=findOverlaps(data2,data3)
hits31=findOverlaps(data3,data1)

data12=as.data.table(hits12)[,.(id1=queryHits,id2=subjectHits)]
data23=as.data.table(hits23)[,.(id2=queryHits,id3=subjectHits)]
data31=as.data.table(hits31)[,.(id3=queryHits,id1=subjectHits)]

##merge overlaps
data123=merge(data12,data23,by=c("id2"),all=TRUE)
data123=merge(data123,data31,by=c("id1","id3"),all=TRUE)
setcolorder(data123,c("id1","id2","id3"))
setorder(data123,id1,id2,id3)
##add missing
data123=merge(data123,data.table(id1=1:length(data1)),by="id1",all=TRUE)
data123=merge(data123,data.table(id2=1:length(data2)),by="id2",all=TRUE)
data123=merge(data123,data.table(id3=1:length(data3)),by="id3",all=TRUE)

##count splitted peaks/DB regions
tt=data123[,.N,by=id1][is.finite(id1),table(N)]
cat("CTCF peaks (WT):",length(data1)," peaks\n")
cat(" ",tt["1"]," peaks unchanged\n")
for(n in names(tt)[names(tt)!="1"])
{
    cat(" ",tt[n]," peaks have been split in ",n," sub-peaks\n")
}
cat(" using ",data123[is.finite(id1),.N]," sub-peaks\n")

tt=data123[,.N,by=id2][is.finite(id2),table(N)]
cat("Cp190 peaks (WT):",length(data2)," peaks\n")
cat(" ",tt["1"]," peaks unchanged\n")
for(n in names(tt)[names(tt)!="1"])
{
    cat(" ",tt[n]," peaks have been split in ",n," sub-peaks\n")
}
cat(" using ",data123[is.finite(id2),.N]," sub-peaks\n")

tt=data123[,.N,by=id3][is.finite(id3),table(N)]
cat("Cp190 differential binding (WT>CTCF0):",length(data3)," differential binding regions\n")
cat(" ",tt["1"]," differential binding regions unchanged\n")
for(n in names(tt)[names(tt)!="1"])
{
    cat(" ",tt[n]," differential binding regions have been split in ",n," sub-regions\n")
}
cat(" using ",data123[is.finite(id3),.N]," sub-regions\n")


data123=data123[,.(in1=is.finite(id1),in2=is.finite(id2),in3=is.finite(id3))]
setorder(data123,in1,in2,in3)
setnames(data123,c("CTCF peak (WT)","Cp190 peak (WT)","Cp190 differential binding (WT>CTCF0)"))

pdf("output/figure_5a.pdf",10,5)
vd= euler(data123)
plot(vd,quantities=list(fontsize=8,col="black"),legend = TRUE)
dev.off()

##save data
fwrite(data123[,.N,by=c("CTCF peak (WT)","Cp190 peak (WT)","Cp190 differential binding (WT>CTCF0)")],file="output/figure_5a.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE,na="NA")


