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

## Lecture 1 (April 11) - Horizontal data integration 

All the materials of the first lecture including slides, leukemia data and the tutorial for **MetaOmics** tool are located in the **Lecture 1** subdirectory.

The lab session about **MetaOmics** pipeline will use **docker**. You are suggested to download and install it before class following the instructions below or more details in the tutorial file 'MetaOmics_Tutorial.pdf':

**Install docker**: 

1. **For Macbook user**, you can directly download and install **docker** from [https://www.docker.com/](https://www.docker.com/);
2. **For Windows user**, make sure you have installed **WSL 2**. See [prerequisites and instructions](https://docs.docker.com/desktop/install/windows-install/) of installing **docker** on Windows. 
   
*Potential trouble shooting on Windows: refer to the video tutorials for [installing **WSL 2**](https://www.youtube.com/watch?v=_fntjriRe48) and [installing **docker**](https://www.youtube.com/watch?v=5RQbdMn04Oc) if you have any trouble. Some PCs may need to further enable virtualization in BIOS, and you can follow [these steps](https://www.simplilearn.com/enable-virtualization-windows-10-article) or [this video](https://www.youtube.com/watch?v=X2fKuPS3yIM).
 
**Run command in terminal to create a container and run MetaOmics tool:**

* `docker pull metaomics/app`
* `docker run --rm --name metaOmics -p 3838:3838 metaomics/app`

Then you can access the MetaOmics pipeline by opening [http://127.0.0.1:3838/metaOmics/](http://127.0.0.1:3838/metaOmics/) on your web browser.

*(Go to [https://docs.docker.com/engine/reference/commandline/run/](https://docs.docker.com/engine/reference/commandline/run/) or use `docker run –help` in terminal for more instructions about **docker**.)*

## Lecture 2 (April 13) - Unsupervised clustering of multi-omics data 

A juptyer lab notebook and corresponding Rscript and pre-processed data for this lab are located in the **Lecture 2** subdirectory. Please have the data downloaded and ready to analyze prior to the lab. The data was originally sourced from the r.jive (omics data) and TCGAbiolinks (clinical data) packages.

For this lab, we’ll be using the R package **MOVICS** which is a wrapper for multiple other packages and thus has many dependencies. In the most current version of R (4.2.3 “Shortstop Beagle”), you can install this package with the following command in R:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
if (!require("devtools")) 
    install.packages("devtools")
    
devtools::install_github("xlucpu/MOVICS")
```

Unfortunately, since this package has many dependencies, it will fail if any of them are unable to load. If you are presented errors that a package is missing, please install that package and try again. It may take several times. **We highly recommend ensuring this package is install prior to the lab and give ample time to ensure all dependencies are loaded.**

Furthermore, **heatmap.plus** has been removed from CRAN and you may need to install it manually from source (version >= 2.5.5). To do so, go to the link below and download the files, then use the following commands to install in R:

heatmap.plus: https://cran.r-project.org/src/contrib/Archive/heatmap.plus/

```
install.packages(path_to_file, repos = NULL, type="source") # path to file is for the package source code you just downloaded
```

If you continue to have trouble, **MOVICS** provides some troubleshooting guidance here: https://github.com/xlucpu/MOVICS

Finally, if you'd like to download our example data directly from Github for the analysis, we'll have code to do so using the package *Rfssa*. Alternatively, you can directly download the "brca_dat.Rdata" file under Lecture 2.

**Update**: Another attendee noted that there may be issues in installing the CIMLR package (a dependency of the MOVICS package). The user received an error regarding their gfortran library, which was installed through the homebrew gcc library on Mac. The user found that adding the following lines to their ~/.R/Makecars file resolved their issue:
 
```
FC = /opt/homebrew/Cellar/gcc/12.2.0/bin/gfortran

F77 = /opt/homebrew/Cellar/gcc/12.2.0/bin/gfortran

FLIBS = -L/opt/homebrew/Cellar/gcc/12.2.0/lib/gcc/12
```

The user referenced the following Stack Overflow exchange that may be able to provide further details: https://stackoverflow.com/questions/29992066/rcpp-warning-directory-not-found-for-option-l-usr-local-cellar-gfortran-4-8/29993906#29993906


## Lecture 3 (April 18) - Dimension reduction for multi-omics data 

RMarkdown code and the corresponding .html report are located in the **Lecture 3** subdirectory. Data is loaded through the 'r.jive' package, and so no data download is required.

If there are issues downloading these files individually, one of the following will likely work: right-clicking the download link and selecting 'Save Link As', or downloading the entire repository as a .zip file and navigating to the relevant documents in the unzipped directory.

The .html file is the preferred lab document; however, if you have issues downloading and viewing the .html document, the .pdf document will also work.

Instructions for obtaining the packages required for the lab portion of the lecture are below.

The following code can be used to install these packages to your default package directory. If you do not have BiocManager already installed, you will need to do so in order to download the TCGAbiolinks package. The other packages are available on CRAN.

install.packages("r.jive")

install.packages("cluster")

install.packages("mclust")

BiocManager::install("TCGAbiolinks")

install.packages("survival")

install.packages("ggplot2")

install.packages("ggfortify")


## Lecture 4 (April 20) - Multi-omics causal mediation analysis and single cell multi-omics analysis

Lecture notes and the materials for the lab session are located in the **Lecture 4** subdirectory.
Similar to **Lecture 2 & 3** the following code can be used to install these packages to your default package directory.

install.packages("HIMA")

install.packages("knitr")

install.packages("ggplot2")

install.packages("ggrepel")

install.packages("Seurat")

install.packages("devtools")

devtools::install_github('satijalab/seurat-data')

install.packages("cowplot")

install.packages("dplyr")
