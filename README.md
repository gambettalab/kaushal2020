Collection of scripts used to prepare figures for the article:

Kaushal et al. In preparation.

## Dependencies

To use these scripts, [R](https://www.R-project.org) (v3.5.1) must be installed along with the following packages:

* [reshape2](https://CRAN.R-project.org/package=reshape2) (v1.4.3).
* [ggplot2](http://ggplot2.org) (v2.2.0).
* [RColorBrewer](https://CRAN.R-project.org/package=RColorBrewer) (v1.1-2).
* [cowplot](https://CRAN.R-project.org/package=cowplot) (v0.9.4).
* [data.table](https://CRAN.R-project.org/package=data.table) (v1.12.2).
* [R.utils](https://CRAN.R-project.org/package=R.utils) (v2.7.0).
* [readxl](https://CRAN.R-project.org/package=readxl) (v1.1.0).
* [eulerr](https://CRAN.R-project.org/package=eulerr) (v6.0.0).
* [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html) (v1.34.0).
* [gridExtra](https://CRAN.R-project.org/package=gridExtra) (v2.3).
* [SDMTools](https://CRAN.R-project.org/package=SDMTools) (v1.1-221).

Version number used to generate the figures in the article are specified above (parenthesis). 

## Installation

Assuming that all dependencies are installed, scripts can be run without installation (apart from downloading and decompressing).

These scripts assume that all input data are located in a directory `data/` in the current working directory. Input data can be downloaded from supplementary information section of the article and should have the following names:

```
data/SupplementaryData1_ChIPseq_CTCF_WT-CTCF0.xlsx          
data/SupplementaryData2_HiC_insulation_scores.xlsx          
data/SupplementaryData3_HiC_CD_boundaries.xlsx      
data/SupplementaryData4_RNAseq_CTCF0-WT.xlsx        
data/SupplementaryData5_HiC_eigenvectors.xlsx       
data/SupplementaryData6_ChIPseq_Cp190_WT-Cp190KO.xlsx   
data/SupplementaryData7_ChIPseq_Cp190_WT-CTCF0.xlsx     
data/SupplementaryData8_ChIPseq_Cp190_CTCF0-Cp190KO.xlsx
data/SupplementaryData9_RNAseq_Cp190KO-WT.xlsx          
```

In addition, internet access is required since the scripts will download the list of blacklisted regions from <https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/dm6-blacklist.v2.bed.gz> ([Amemiya  et al. Scientific Reports (2019)](https://doi.org/10.1038/s41598-019-45839-z) and [ENCODE Project Consortium, Nature (2012)](https://doi.org/10.1038/nature11247) )


## Usage

Scripts can be run from the shell (e.g. bash):

```
bin/figure2.R
```

```
bin/figure3.R
```

```
bin/figure4.R
```

```
bin/figure5a.R
```

```
bin/figure5.R
```

All output will be saved in the directory `output/`. Total runtime should be approximately 10 minutes on a desktop computer.


