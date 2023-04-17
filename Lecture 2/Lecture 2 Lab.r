## to install MOVICS
### note: there are many dependencies; you may get an error
### for a missing package; download and try again
### it may take several times
### you can refer to their IMPORTS file from their DESCRIPTION file on Github
### for a list of dependencies:
### https://github.com/xlucpu/MOVICS/blob/master/DESCRIPTION

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# if (!require("devtools")) 
#     install.packages("devtools")
# devtools::install_github("xlucpu/MOVICS")

install.packages("~/Downloads/SNFtool_2.3.0.tar.gz", repos = NULL, type="source") 



library(MOVICS)
set.seed(4444)

library("jpeg")
jj <- readJPEG("MOVICS_pipeline.jpeg",native=TRUE)
plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
rasterImage(jj,0,0,1,1)


jj <- readJPEG("methods_comparison.jpeg",native=TRUE)
plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
rasterImage(jj,0,0,1,1)


# install.packages("Rfssa") if you want to download through github
library(Rfssa)
url <- "https://github.com/KechrisLab/ASAShortCourse-MultiOmics/blob/main/Lecture%202/brca_dat.Rdata"
load_github_data(url)

# load("brca_dat.Rdata")

# let's get a quick look at our data
names(brca_dat)
paste("dim of clinical data:", dim(brca_dat[["clinical"]]))
head(brca_dat[["clinical"]])

# check sample names all match
identical(brca_dat[["clinical"]]$bcr_patient_barcode, colnames(brca_dat[["MO"]][["Expression"]]))
identical(brca_dat[["clinical"]]$bcr_patient_barcode, colnames(brca_dat[["MO"]][["Methylation"]]))
identical(brca_dat[["clinical"]]$bcr_patient_barcode, colnames(brca_dat[["MO"]][["miRNA"]]))
                                                                                
identical(colnames(brca_dat[["MO"]][["Expression"]]), colnames(brca_dat[["MO"]][["Methylation"]]))
identical(colnames(brca_dat[["MO"]][["Expression"]]), colnames(brca_dat[["MO"]][["miRNA"]]))
identical(colnames(brca_dat[["MO"]][["Methylation"]]), colnames(brca_dat[["MO"]][["miRNA"]]))
# should all be TRUE (6)

paste("names of MO data:", names(brca_dat[["MO"]]))
paste("dim of mRNA data:", dim(brca_dat[["MO"]][["Expression"]]))
paste("dim of methylation data:", dim(brca_dat[["MO"]][["Methylation"]]))
paste("dim of miRNA data:", dim(brca_dat[["MO"]][["miRNA"]]))

# data checking -- are there any missing values?
sum(is.na(brca_dat[["MO"]][["Expression"]]))
sum(is.na(brca_dat[["MO"]][["Methylation"]]))
sum(is.na(brca_dat[["MO"]][["miRNA"]]))

# exp
range(brca_dat[["MO"]][["Expression"]])
plot(density(brca_dat[["MO"]][["Expression"]]), main = "Expression")

# methyl
range(brca_dat[["MO"]][["Methylation"]])
plot(density(brca_dat[["MO"]][["Methylation"]]), main = "Methylation")

# miRNA
range(brca_dat[["MO"]][["miRNA"]])
plot(density(brca_dat[["MO"]][["miRNA"]]), main = "miRNA")


jj <- readJPEG("ex_breastcancer_pic.jpeg")
plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
rasterImage(jj,0,0,1,1)
# https://doi.org/10.1016/B978-0-12-800886-7.00021-2

optk.brca <- getClustNum(data        = brca_dat[["MO"]],
                         is.binary   = c(F,F,F), # all omics data is continuous (not binary)
                         try.N.clust = 2:8, # try cluster number from 2 to 8
                         fig.name    = "CLUSTER NUMBER OF TCGA-BRCA")



optk.brca

# what if we use the suggested k=3 clusters?
# you don't need to run this during lab; I'm just presenting it as an example
mo_rslts_3 <- getMOIC(data = brca_dat[["MO"]],
                         methodslist = list("SNF", "PINSPlus", "NEMO", "LRAcluster", "IntNMF"),
                         N.clust     = 3,
                         type        = c("gaussian", "gaussian", "gaussian"))

cmoic.brca_3 <- getConsensusMOIC(moic.res.list = mo_rslts_3,
                               fig.name      = "CONSENSUS HEATMAP - 3 Clusters",
                               distance      = "euclidean",
                               linkage       = "average")

getSilhouette(sil      = cmoic.brca_3$sil, # a sil object returned by getConsensusMOIC()
              fig.path = getwd(),
              fig.name = "SILHOUETTE",
              height   = 5.5,
              width    = 5)

mo_rslts <- getMOIC(data = brca_dat[["MO"]],
                         methodslist = list("SNF", "PINSPlus", "NEMO", "LRAcluster", "IntNMF"),
                         N.clust     = 4, # set number of clusters
                         type        = c("gaussian", "gaussian", "gaussian")) # what is the distribution of the datasets in MO list (same order)


cmoic.brca <- getConsensusMOIC(moic.res.list = mo_rslts,
                               fig.name      = "CONSENSUS HEATMAP - 4 Clusters",
                               distance      = "euclidean",
                               linkage       = "ward.D")

getSilhouette(sil      = cmoic.brca$sil, # a sil object returned by getConsensusMOIC()
              fig.path = getwd(),
              fig.name = "SILHOUETTE",
              height   = 5.5,
              width    = 5)

# convert beta value to M value for stronger signal
std_dat <- brca_dat[["MO"]]
std_dat[["Methylation"]] <- log2(std_dat[["Methylation"]] / (1 - std_dat[["Methylation"]]))

# data normalization for heatmap
plotdata <- getStdiz(data       = std_dat,
                     halfwidth  = c(2,2,2), # no truncation for mutation
                     centerFlag = c(T,T,T), # no center for mutation
                     scaleFlag  = c(T,T,T)) # no scale for mutation

mRNA.col   <- c("#00FF00", "#008000", "#000000", "#800000", "#FF0000")
meth.col   <- c("#0074FE", "#96EBF9", "#FEE900", "#F00003")
miRNA.col <- c("#6699CC", "white"  , "#FF3C38")
col.list   <- list(mRNA.col, meth.col, miRNA.col)

# extract PAM50, pathologic stage for sample annotation
annCol    <- brca_dat[["clinical"]][,c("BRCA_Subtype_PAM50"), drop = FALSE]

# generate corresponding colors for sample annotation
annColors <- list(
                  BRCA_Subtype_PAM50  = c("Basal" = "blue",
                            "Her2"   = "red",
                            "LumA"   = "yellow",
                            "LumB"   = "green",
                            "Normal" = "black")
                )



# comprehensive heatmap
getMoHeatmap(data          = plotdata,
             row.title     = names(std_dat),
             is.binary     = c(F,F,F), # we don't have any binary omics data (ex mutation)
             legend.name   = c("mRNA","M value","miRNA"),
             clust.res     = mo_rslts$SNF$clust.res, # cluster results for SNF
             color         = col.list,
             # annCol        = annCol, # annotation for samples (if you want to show PAM50 classes too)
             # annColors     = annColors, # annotation color
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP OF SNF")

getMoHeatmap(data          = plotdata,
             row.title     = names(std_dat),
             is.binary     = c(F,F,F), # all data is continuous
             legend.name   = c("mRNA","M value","miRNA"),
             clust.res     = mo_rslts$PINSPlus$clust.res, # cluster results for PINSPlus
             color         = col.list,
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP OF PINSPlus")

# comprehensive heatmap (may take a while)
getMoHeatmap(data          = plotdata,
             row.title     = names(std_dat),
             is.binary     = c(F,F,F),
             legend.name   = c("mRNA","M value","miRNA"),
             clust.res     = mo_rslts$NEMO$clust.res, # cluster results for NEMO
             color         = col.list,
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP OF PINSPlus")

# comprehensive heatmap (may take a while)
getMoHeatmap(data          = plotdata,
             row.title     = names(std_dat),
             is.binary     = c(F,F,F),
             legend.name   = c("mRNA","M value","miRNA"),
             clust.res     = mo_rslts$LRAcluster$clust.res, # cluster results for LRAcluster
             color         = col.list,
             width         = 10, 
             height        = 5, 
             fig.name      = "COMPREHENSIVE HEATMAP OF PINSPlus")

# comprehensive heatmap (may take a while)
getMoHeatmap(data          = plotdata,
             row.title     = names(std_dat),
             is.binary     = c(F,F,F),
             legend.name   = c("mRNA","M value","miRNA"),
             clust.res     = mo_rslts$IntNMF$clust.res, # cluster results for intNMF
             color         = col.list,
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP OF IntNMF")

getMoHeatmap(data          = plotdata,
             row.title     = names(plotdata),
             is.binary     = c(F,F,F), # no binary omics data
             legend.name   = c("mRNA","M value","miRNA"),
             clust.res     = cmoic.brca$clust.res, # consensusMOIC results
             clust.dend    = NULL, # show no dendrogram for samples
             show.colnames = FALSE, # show no sample names
             show.row.dend = c(T,T,T), # show dendrogram for features
             annRow        = NULL, # no selected features
             color         = col.list,
             annCol        = annCol, # annotation for samples
             annColors     = annColors, # annotation color
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP OF CONSENSUSMOIC")

clust_rslts_SNF_df <- data.frame(mo_rslts$SNF$clust.res)
colnames(clust_rslts_SNF_df) <- c("samID", "SNF")
clust_rslts_CM_df <- data.frame(cmoic.brca$clust.res)
colnames(clust_rslts_CM_df) <- c("samID", "Consensus")

clust_rslts_df <- merge(clust_rslts_SNF_df, clust_rslts_CM_df, by="samID")
# head(clust_rslts_df)

table(clust_rslts_df$SNF, clust_rslts_df$Consensus)

# survival comparison
brca_dat[["clinical"]]$futime = as.numeric(brca_dat[["clinical"]]$futime)
head(brca_dat[["clinical"]])
surv.brca <- compSurv(moic.res         = cmoic.brca,
                      surv.info        = brca_dat[["clinical"]],
                      convt.time       = "m", # convert day unit to month
                      surv.median.line = "h", # draw horizontal line at median survival
                      xyrs.est         = c(5,10), # estimate 5 and 10-year survival
                      fig.name         = "KAPLAN-MEIER CURVE OF CONSENSUSMOIC")


print(surv.brca)

# clinVars_df <- brca_dat[["clinical"]][,c("ajcc_pathologic_stage", "age_at_diagnosis","ajcc_pathologic_t", "ajcc_pathologic_n","ajcc_pathologic_m")]
clinVars_df <- brca_dat[["clinical"]]
rownames(clinVars_df) <- clinVars_df$bcr_patient_barcode
clinVars_df <- clinVars_df[,c( "ajcc_pathologic_stage", "ajcc_pathologic_t", "ajcc_pathologic_n", "ajcc_pathologic_m", "age_at_diagnosis")]
head(clinVars_df)


clin.brca <- compClinvar(moic.res      = cmoic.brca,
                         var2comp      = clinVars_df, # data.frame needs to summarize (must has row names of samples)
                         strata        = "Subtype", # stratifying variable (e.g., Subtype in this example)
                         # factorVars    = c("ajcc_pathologic_stage"), # features that are considered categorical variables
                         factorVars    = c("ajcc_pathologic_stage", "ajcc_pathologic_t", "ajcc_pathologic_n", "ajcc_pathologic_m"), # features that are considered categorical variables
                         nonnormalVars = "age_at_diagnosis", # feature(s) that are considered using nonparametric test
                         exactVars     = NULL, # feature(s) that are considered using exact test
                         doWord        = FALSE, # generate .docx file in local path
                         tab.name      = "SUMMARIZATION OF CLINICAL FEATURES")
clin.brca

# compare agreement with other subtypes
sub_df = data.frame(
                    BRCA_Subtype_PAM50 = brca_dat[["clinical"]][,c("BRCA_Subtype_PAM50")])
rownames(sub_df) = brca_dat[["clinical"]][,c("bcr_patient_barcode")]
head(sub_df)

# agreement comparison (support up to 6 classifications include current subtype)
agree.brca <- compAgree(moic.res  = cmoic.brca,
                        subt2comp = sub_df,
                        doPlot    = TRUE,
                        box.width = 0.2,
                        fig.name  = "AGREEMENT OF CONSENSUSMOIC WITH PAM50 Subtype")

# run DEA with limma
runDEA(dea.method = "limma",
       expr       = brca_dat[["MO"]][["Expression"]], # normalized expression data
       moic.res   = cmoic.brca,
       overwt = T,
       res.path   = getwd(), # path to save marker files
       prefix     = "de_TCGA-BRCA")


# choose limma result to identify subtype-specific DOWN-regulated biomarkers
marker.dn <- runMarker(moic.res      = cmoic.brca,
                       dea.method    = "limma",
                       prefix        = "de_TCGA-BRCA",
                       dirct         = "down",
                       dat.path = getwd(),
                       res.path = getwd(),
                       p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                       p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs                       
                       n.marker      = 200, # number of biomarkers for each subtype
                       doplot        = T,
                       annCol        = annCol,
                       annColors     = annColors,
                       norm.expr     = brca_dat[["MO"]][["Expression"]],
                       fig.name      = "UPREGULATED BIOMARKER HEATMAP")


# subtype-specific UP-regulated biomarkers
marker.up <- runMarker(moic.res      = cmoic.brca,
                       dea.method    = "limma",
                       prefix        = "de_TCGA-BRCA",
                       dirct         = "up",
                       dat.path = getwd(),
                       res.path = getwd(),
                       p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                       p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs                       
                       n.marker      = 200, # number of biomarkers for each subtype
                       doplot        = T,
                       annCol        = annCol,
                       annColors     = annColors,
                       norm.expr     = brca_dat[["MO"]][["Expression"]],
                       fig.name      = "UPREGULATED BIOMARKER HEATMAP")


# MUST locate ABSOLUTE path of msigdb file
MSIGDB.FILE <- system.file("extdata", "c5.bp.v7.1.symbols.xls", package = "MOVICS", mustWork = TRUE)

# # run GSEA to identify DOWN-regulated GO pathways using results from edgeR
gsea.down <- runGSEA(moic.res     = cmoic.brca,
                   dea.method   = "limma", # name of DEA method
                   prefix       = "de_TCGA-BRCA", # MUST be the same of argument in runDEA()
                   dat.path     = getwd(), # path of DEA files
                   res.path     = getwd(), # path to save GSEA files
                   msigdb.path  = MSIGDB.FILE, # MUST be the ABSOLUTE path of msigdb file
                   norm.expr    = brca_dat[["MO"]][["Expression"]], # use normalized expression to calculate enrichment score
                   dirct        = "down", # direction of dysregulation in pathway
                   p.cutoff     = 0.05, # p cutoff to identify significant pathways
                   p.adj.cutoff = 0.25, # padj cutoff to identify significant pathways
                   gsva.method  = "gsva", # method to calculate single sample enrichment score
                   norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                   fig.name     = "DOWNREGULATED PATHWAY HEATMAP")


data.frame(gsea.down$gsea.list$CS2[1:6,3:6])

head(round(gsea.down$grouped.es,3))

# # run GSEA to identify up-regulated GO pathways using results from limma
gsea.up <- runGSEA(moic.res     = cmoic.brca,
                   dea.method   = "limma", # name of DEA method
                   prefix       = "detesting_TCGA-BRCA", # MUST be the same of argument in runDEA()
                   dat.path     = getwd(), # path of DEA files
                   res.path     = getwd(), # path to save GSEA files
                   msigdb.path  = MSIGDB.FILE, # MUST be the ABSOLUTE path of msigdb file
                   norm.expr    = brca_dat[["MO"]][["Expression"]], # use normalized expression to calculate enrichment score
                   dirct        = "up", # direction of dysregulation in pathway
                   p.cutoff     = 0.05, # p cutoff to identify significant pathways
                   p.adj.cutoff = 0.25, # padj cutoff to identify significant pathways
                   gsva.method  = "gsva", # method to calculate single sample enrichment score
                   norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                   fig.name     = "UPREGULATED PATHWAY HEATMAP")


# MUST locate ABSOLUTE path of gene set file
GSET.FILE <- 
  system.file("extdata", "gene sets of interest.gmt", package = "MOVICS", mustWork = TRUE)

# run GSVA to estimate single sample enrichment score based on given gene set of interest
gsva.res <- 
  runGSVA(moic.res      = cmoic.brca,
          norm.expr     = brca_dat[["MO"]][["Expression"]],
          gset.gmt.path = GSET.FILE, # ABSOLUTE path of gene set file
          gsva.method   = "gsva", # method to calculate single sample enrichment score
          annCol        = annCol,
          annColors     = annColors,
          fig.path      = getwd(),
          fig.name      = "GENE SETS OF INTEREST HEATMAP",
          height        = 5,
          width         = 10)


message("check raw enrichment score")
gsva.res$raw.es[1:3,1:3]

message("check z-scored and truncated enrichment score")
gsva.res$scaled.es[1:3,1:3]

