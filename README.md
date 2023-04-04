# ASAShortCourse-MultiOmics



![Picture1](https://user-images.githubusercontent.com/6655031/229678834-9fec0e0b-042d-40bc-85ce-5e42d72ea864.jpg)

Figure from Tarazona et al. (2021) Nature Computational Science

## Directors:
* George Tseng (University of Pittsburgh)
* Katerina Kechris (University of Colorado Anschutz Medical Campus)

## Instructors: 
* Wenjia Wang (University of Pittsburgh)
* Sierra Niemiec (University of Colorado Anschutz Medical Campus)
* Jack Pattee (University of Colorado Anschutz Medical Campus)
* Rick Chang (University of Pittsburgh) 

## All sessions are 3:00 – 4:30 PM (EASTERN)

## Session 1 (April 11) - Horizontal data integration 

## Session 2 (April 13) - Unsupervised clustering of multi-omics data 

## Session 3 (April 18) - Dimension reduction for multi-omics data 

RMarkdown code and the corresponding .html report are located in the 'Lecture 3' subdirectory. Data is loaded through the 'r.jive' package, and so no data download is required.

Instructions for obtaining the packages required for the lab portion of the lecture are below.

The following code can be used to install these packages to your default package directory. If you do not have BiocManager already installed, you will need to do so in order to download the TCGAbiolinks package. The other packages are available on CRAN.

install.packages("r.jive")
install.packages("cluster")
install.packages("mclust")
BiocManager::install("TCGAbiolinks")
install.packages("survival")
install.packages("ggplot2")
install.packages("ggfortify")


## Session 4 (April 20) - Multi-omics causal mediation analysis and single cell multi-omics analysis 
