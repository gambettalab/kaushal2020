#!/usr/bin/env Rscript


## Generic script to analyse Cell analyzer results of transected cells with a reporter construct          
## Expected input : raw data saved from Novocyte machine (.csv) placed in /data folder 
## Or in this example an Excel spreadsheet containig each sample in separte sheets
## Output files produced in /outut folder:                                                      
## - Fluorescence scatter plot with gating indication for each sample as shown in Fig. S4a  
## - Violin plot of fluorescence signal and ratio distribution as shown in Fig. 4b and S4b       
## - Mean and median table of ratio distributions used as value of insulator strength in Fig. 4c 

##################################################################################################

## Copyright (C) 2020 by Pascal Cousin, Centre for Integrative Genomics, University of Lausanne, Switzerland.
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






# Load some libraries
library(SDMTools)     # install from https://CRAN.R-project.org/package=SDMTools if missing
library(ggplot2)      # use install.packages("ggplot2") if missing
require(gridExtra)    # use install.packages("gridExtra") if missing
library(readxl)
# Load raw data exported from Novocyte

path_to_data <- "./data" 
path_to_result <- "./output"

# file list only for this exemple
# files <- c("C11_Spacer(rep1).csv","C12_Spacer(rep2).csv","D1_Gypsy(rep1).csv","D2_Gypsy(rep2).csv")
# Uncomment next line for a more flexible file list
# files <- list.files(path = path_to_data, pattern = ".csv") # Raw Novocyte files expected in /data folder in csv format
files <- excel_sheets(paste(path_to_data,"Fig4b_Cell_analyzer_data.xlsx",sep="/"))
# Initialise some lists
AllDataTemp <- list()
AllData <- list()
S2_Gated_data <- list()
Final_Gated_data <- list()
ratio <- list()
final_ratio <- list()
sample_name  <- list()
unique_file_name <- list()
# Define gates
S2_gate = cbind(x=c(log10(1607914),log10(1154420),log10(2876799),log10(6397704),log10(5455329)),y=c(log10(17412),log10(83087),log10(345638),log10(424628),log10(62667)))
Final_gate = cbind(x=c(0,3,5,4.3,3,0),y=c(0.0,0.0,2,4.8,7,5))
Final_gate.df <- data.frame(Final_gate)
# Extract sample name from file name 
count <- 0
for (sample in files) {
  count = count + 1
  temp <- strsplit(files[count],"_")
  sample_name[count] <- strsplit(temp[[1]][2],".csv")
}

count <-0
for (name in unique(sample_name)){
  count <- count +1
  unique_file_name[[count]] <- name
}
# Extract data and apply gate
count <- 0
for (sample in files) {
  count = count + 1
  # read raw data from .csv
  # temp_file <- read.csv(paste(path_to_data,sample, sep="/")) 
  # read raw data from .xlsx
  # temp_file <- read_excel(paste(path_to_data,"Fig4b_Cell_analyzer_data.xlsx",sep="/"), sheet = sample, .name_repair = "universal" )
  temp_file <- data.frame(read_excel(paste(path_to_data,"Fig4b_Cell_analyzer_data.xlsx",sep="/"), sheet = sample), check.names = "TRUE")
  # select columns containing usfull detection chanels
  AllDataTemp[[count]] <- temp_file[1:4] 
  # remove aberent negative values
  AllData[[count]] <- subset(AllDataTemp[[count]],AllDataTemp[[count]]$PE.Texas.Red.H>0 & AllDataTemp[[count]]$FITC.H>0) 
  # Select S2 cells and keep
  S2_gated = pnt.in.poly(cbind(log10(AllData[[count]]$FSC.H),log10(AllData[[count]]$SSC.H)),S2_gate) 
  AllData[[count]]$S2gated <- S2_gated$pip
  S2_Gated_data[[count]] <- subset(AllData[[count]],AllData[[count]]$S2gated == 1 )
  # select untransfected cells and remove
  untransfeced = pnt.in.poly(cbind(log10(S2_Gated_data[[count]]$FITC.H),log10(S2_Gated_data[[count]]$PE.Texas.Red.H)),Final_gate) 
  S2_Gated_data[[count]]$masked <- untransfeced$pip
  Final_Gated_data[[count]] <- subset(S2_Gated_data[[count]], S2_Gated_data[[count]]$masked ==0)
  
  
  
  # Scatter plot ------------------------------------------------------------
   df <- data.frame(x = log10(S2_Gated_data[[count]]$FITC.H), y = log10(S2_Gated_data[[count]]$PE.Texas.Red.H),d = densCols(log10(S2_Gated_data[[count]]$FITC.H), log10(S2_Gated_data[[count]]$PE.Texas.Red.H),nbin = 2000, colramp = colorRampPalette(rev(rainbow(10, end = 4/6)))))
  plot_name <- paste(path_to_result, paste(paste("Scatter_plot",strsplit(sample, "_")[[1]][2],strsplit(sample, "_")[[1]][3],sep="_"),"pdf",sep="."),sep="/")
  pdf(plot_name, height = 7, width = 7)
  d <- ggplot(df) +
    geom_point(aes(log10(S2_Gated_data[[count]]$FITC.H), log10(S2_Gated_data[[count]]$PE.Texas.Red.H), col = d), size = 0.01) +
    scale_color_identity() +
    geom_polygon(data = Final_gate.df, aes(x=x, y=y), colour="black", fill="NA")+
    ggtitle(sample_name[count]) + theme_classic() + scale_x_continuous(name= "log10(eGFP)",expand = c(0, 0), limits = c(0,8)) + 
    scale_y_continuous(name="log10(mCherry)",expand = c(0, 0), limits = c(0,8)) +
    theme(axis.ticks.length = unit(7, "pt"))
  print(d)
  dev.off()
  
  # S2 gating plot ----------------------------------------------------------
  # optional, only for quality check
  
  # plot(log10(AllData[[count]]$FSC.H),log10(AllData[[count]]$SSC.H),col=AllData[[count]]$S2gated +2, pch=".",main = sample_name[count],xlab = "log10 FSC",ylab =  "log10 SSC",xlim = c(5.8,7),ylim=c(3,6))
   
  # polygon(S2_gate)
  
  
  # Ratio log2(mCherry/eGFP)
  Final_Gated_data[[count]]$ratio  <- log2((Final_Gated_data[[count]]$PE.Texas.Red.H)/(Final_Gated_data[[count]]$FITC.H))
  Final_Gated_data[[count]]$name  <- as.factor(sample_name[count])
  Final_Gated_data[[count]]$sample_num <- count
  
  
}


# Violin plot extention
# Adapted from https://rdrr.io/cran/UCSCXenaShiny/src/R/GeomSplitViolin.R

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, draw_group = function(self, data, ..., draw_quantiles = NULL){
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1,'group']
  newdata <- plyr::arrange(transform(data, x = if(grp%%2==1) xminv else xmaxv), if(grp%%2==1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 
                                              1))
    quantiles <- create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})

geom_split_violin <- function (mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, show.legend = show.legend, inherit.aes = inherit.aes, params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}
median.quartile <- function(x){
  out <- quantile(x, probs = c(0.25,0.5,0.75))
  names(out) <- c("ymin","y","ymax")
  return(out) 
}

# Prepare data for the violin plot
all_sample <- do.call(rbind,Final_Gated_data)
m_Cherry <- all_sample[,c(4,8)]
m_Cherry$dye <- as.factor("mCherry")
names(m_Cherry)[1] <- "fluo"
GFP <- all_sample[,c(3,8)]
GFP$dye <- as.factor("eGFP")
names(GFP)[1] <- "fluo"
split_plot.df <- rbind(m_Cherry,GFP)

# Violin plots
pdf(paste(path_to_result,"Fig. 4b.pdf",sep ="/"), height = 7, width = 7)
pp  <- ggplot(split_plot.df,aes(name,log10(fluo),  fill = dye)) +  
  stat_summary(fun.data=median.quartile, geom="crossbar", width=0.15 ,position = "dodge") +
  ylab("log10(fluo)") + theme_classic() + geom_split_violin( trim=FALSE, alpha=0.3) +
  scale_fill_manual(values = c("red", "green")) + theme(legend.position='none')


p2 <- ggplot(all_sample,aes(x=name,y=ratio)) +geom_violin(fill = grey(0.8), trim = TRUE) + 
  stat_summary(fun.data=median.quartile, geom="crossbar", width=0.05 ) +
  scale_fill_grey(start = 0, end = 1, aesthetics = "fill") +
  geom_boxplot(width = 0.05, outlier.shape = NA) + ylab("log2(mCherry/GFP)") +theme_classic()

grid.arrange(pp, p2, nrow = 2)
dev.off()

# Save mean and median values of log2(mCherry/GFP)

count <- 0
median_table = list()
for (sample in files) {
  count <- count+1
  median_table[[count]] <- c(sample,summary(Final_Gated_data[[count]]$ratio)[3],summary(Final_Gated_data[[count]]$ratio)[4])
}
final_table = do.call(rbind,median_table)
write.csv(final_table,paste(path_to_result,"Mean_median_table.csv", sep="/"))

