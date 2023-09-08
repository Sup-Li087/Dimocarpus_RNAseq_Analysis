################################################################################
##### This work used several tutorials
### 1. TutorialRNA differential expression analysis
### https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html
### 2. https://bioinformatics-core-shared-training.github.io/RNAseq-R/rna-
### seq-preprocessing.nb.html
### 3. very useful paper about design and contrasts
### https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/
### inst/doc/designmatrices.html#design-and-contrast-matrices
### 4. Workflow article
### https://f1000research.com/articles/5-1438
### number 4 is also used for generating heatmap plots 
### 5. EdgeR user guide
### https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc
### /edgeRUsersGuide.pdf
################################################################################

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")

BiocManager::install("edgeR")
BiocManager::install("baySeq")

library(edgeR)
library(baySeq)

###############################################################################
### Steps
## 1. Make appropriate dataframe - I still leave the replicates as it is
## 2. Subset to make each file for comparing each pairs
################################################################################
##### import read counts
count <- read.csv("All_annot for R.csv", header = T)
str(count)

matrix <- as.matrix(count[,-1])
head(matrix)
groups <- c("D1","D1","D6_N", "D6_N", "D6_F", "D6_F", "D8_N",
            "D8_N", "D8_M", "D8_M", "D8_F", "D8_F")

levels(groups)
groups <- as.factor(groups)

genes <- (count$gene_id)

library(edgeR)
d <- DGEList(counts = matrix, group = groups, genes = genes)
d
d.full$samples

################################################################################
### filtering data
################################################################################
dim(d)

d.full <- d ## keep the old one in case we mess up

head(d$counts)

head(cpm(d))

apply(d.full$counts, 2, sum) # total gene counts per sample

#There are a few ways to filter out lowly expressed genes. 
#When there are biological replicates in each group, 
#in this case we have a sample size of 2 in each group, 
#we favour filtering on a minimum counts per million threshold present 
#in at least 2 samples. Two represents the smallest sample size 
#for each group in our experiment. 
#In this dataset, we choose to retain genes if they are expressed at 
#a counts-per-million (CPM) above 0.5 in at least two samples.
#
#We’ll use the cpm function from the edgeR library (M D Robinson, 
#McCarthy, and Smyth 2010) to generate the CPM values and then filter. 
#Note that by converting to CPMs we are normalising for the different 
#sequencing depths for each sample.

#A CPM of 0.5 is used as it corresponds to a count of 10-15 
#for the library sizes in this data set. If the count is any smaller, 
#it is considered to be very low, indicating that the associated gene 
#is not expressed in that sample. A requirement for expression in two 
#or more libraries is used as each group contains two replicates. 
#This ensures that a gene will be retained if it is only expressed 
#in one group. Smaller CPM thresholds are usually appropriate for 
#larger libraries. As a general rule, a good threshold can be chosen 
#by identifying the CPM that corresponds to a count of 10, 
#which in this case is about 0.5. You should filter with CPMs rather than 
#filtering on the counts directly, as the latter does not account 
#for differences in library sizes between samples.

## in our case cpm 0.65

keep <- rowSums(cpm(d.full) > 0.65) >= 2
d <- d.full[keep,]
dim(d)

### reset the library sizes
d$samples$lib.size <- colSums(d$counts)
d$samples

################################################################################
####Normalizing the data
################################################################################
d <- calcNormFactors(d)
d
################################################################################
#### Data Exploration
## Before proceeding with the computations for differential expression, 
## it is possible to produce a plot showing the sample relations 
## based on multidimensional scaling.
################################################################################
plotMDS(d, col = as.numeric(d$samples$group))

plotMD(d, column = 1)
abline(h=0, col="red", lty=2, lwd=2)

### Estimating the Dispersion
##
### The first major step in the analysis of DGE data using the NB model 
### is to estimate the dispersion parameter for each tag, 
### a measure of the degree of inter-library variation for that tag. 
### Estimating the common dispersion gives an idea of overall variability 
### across the genome for this dataset.
################################################################################
d1 <- estimateCommonDisp(d, verbose = T, robust = T)
names(d1)

d1 <- estimateTagwiseDisp(d1)
names(d1)

plotBCV(d1)

################################################################################
#### GLM estimares of dispersion
#### here we used the model without the intercept term. 
#### In our case the variable is factor so the models 
#### with or without intercept is the same
################################################################################

design.mat <- model.matrix(~ 0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)
levels(groups) ## D1 = control = 1st entry OK

d2 <- estimateGLMCommonDisp(d, design.mat)
d2 <- estimateGLMTrendedDisp(d2,design.mat, method="power")

# You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" 
# otherwise.

d2 <- estimateGLMTagwiseDisp(d2, design.mat)
plotBCV(d2)
################################################################################
#### Differential expression
#### Classic approach
################################################################################
et12 <- exactTest(d1, pair=c(1,2)) # compare groups D1 control and D6F
et13 <- exactTest(d1, pair=c(1,3)) # compare groups D1 control and D6N
et23 <- exactTest(d1, pair=c(2,3)) # compare groups D6F and D6N
et14 <- exactTest(d1, pair=c(1,4)) # compare groups D1 control and D8F
et15 <- exactTest(d1, pair=c(1,5)) # compare groups D1 control and D8M
et16 <- exactTest(d1, pair=c(1,6)) # compare groups D1 control and D8N
et46 <- exactTest(d1, pair=c(4,6)) # compare groups D8N and D8F
et56 <- exactTest(d1, pair=c(5,6)) # compare groups D8N and D8M
et25 <- exactTest(d1, pair=c(2,5)) # compare groups D6F and D8M
et24 <- exactTest(d1, pair=c(2,4)) # compare groups D6F and D8F


topD1.D6F <- topTags(et12)
topTags(et13, n=10)
topTags(et23, n=10)
topTags(et14, n=10)

################################################################################
### Testing for differential expression
### First, we fit genewise glms:
### GLM approach
################################################################################
fit <- glmQLFit(d1, design.mat)

my.contrasts <- makeContrasts(D1vsD6_F = D1-D6_F, 
                              D1vsD6_N = D1-D6_N, 
                              D1vsD8_F = D1-D8_F, 
                              D1vsD8_M = D1-D8_M,
                              D1vsD8_N = D1-D8_N,
                              D6_NvsD6_F = D6_N-D6_F,
                              D8_NvsD8_M = D8_N-D8_M,
                              D8_NvsD8_F = D8_N-D8_F,
                              D6_FvsD8_M = D6_F-D8_M,
                              D6_FvsD8_F = D6_F-D8_F,
                              levels = design.mat)


qlf.d1d6f <- glmQLFTest(fit, contrast = my.contrasts[,"D1vsD6_F"])
topTags(qlf.d1d6f)
all.d1d6f <- topTags(qlf.d1d6f, n = Inf)
up.d1d6f <- head(all.d1d6f[all.d1d6f$table$logFC > 0,], n= 20L)
write.csv(up.d1d6f , "d1 vs d6f.csv")
is.de <- decideTestsDGE(qlf.d1d6f)
summary(is.de)

plotMD(qlf.d1d6f, status = is.de, values=c(1,-1), col=c("red","blue"),
       legend="topright")

### try to reduce the number of expressed gene
tr <- glmTreat(fit, contrast = my.contrasts[,"D1vsD6_F"], 
               lfc=log2(1.5))
topTags(tr)
is.de <- decideTestsDGE(tr)
summary(is.de)

plotMD(tr, status = is.de, values=c(1,-1), col=c("red","blue"),
legend="topright")
table(keep)

################################################################################
### alternative approach
### We use this for further plotting
################################################################################
qlf.d1d6f.tr <- glmQLFTest(tr, contrast = my.contrasts[,"D1vsD6_F"])
topTags(qlf.d1d6f)


qlf.d1d6n <- glmQLFTest(fit, contrast = my.contrasts[,"D1vsD6_N"])
all.d1d6n <- topTags(qlf.d1d6n, n = Inf)
up.d1d6n <- head(all.d1d6n[all.d1d6n$table$logFC > 0,], n= 20L)
write.csv(up.d1d6n, "d1 vs d6n.csv")

qlf.d1d8f <- glmQLFTest(fit, contrast = my.contrasts[,"D1vsD8_F"])
all.d1d8f <- topTags(qlf.d1d8f, n = Inf)
up.d1d8f <- head(all.d1d8f[all.d1d8f$table$logFC > 0,], n= 20L)
write.csv(up.d1d8f, "d1 vs d8f.csv")

qlf.d1d8m <- glmQLFTest(fit, contrast = my.contrasts[,"D1vsD8_M"])
all.d1d8m <- topTags(qlf.d1d8m, n = Inf)
up.d1d8m <- head(all.d1d8m[all.d1d8m$table$logFC > 0,], n= 20L)
write.csv(up.d1d8m , "d1 vs d8m.csv")

qlf.d1d8n <- glmQLFTest(fit, contrast = my.contrasts[,"D1vsD8_N"])
all.d1d8n <- topTags(qlf.d1d8n, n = Inf)
up.d1d8n<- head(all.d1d8n[all.d1d8n$table$logFC > 0,], n= 20L)
write.csv(up.d1d8n, "d1 vs d8n.csv")

qlf.d6fd6n <- glmQLFTest(fit, contrast = my.contrasts[,"D6_NvsD6_F"])
all.d6fd6n <- topTags(qlf.d6fd6n, n = Inf)
up.d6fd6n <- head(all.d6fd6n[all.d6fd6n$table$logFC > 0,], n= 20L)
write.csv(up.d6fd6n, "d6n vs d6f.csv")

qlf.d8nd8m <- glmQLFTest(fit, contrast = my.contrasts[,"D8_NvsD8_M"])
all.d8nd8m <- topTags(qlf.d8nd8m, n = Inf)
up.d8nd8m<- head(all.d8nd8m[all.d8nd8m$table$logFC > 0,], n= 20L)
write.csv(up.d8nd8m, "d8n vs d8m.csv")

qlf.d8nd8f <- glmQLFTest(fit, contrast = my.contrasts[,"D8_NvsD8_F"])
all.d8nd8f  <- topTags(qlf.d8nd8f , n = Inf)
up.d8nd8f <- head(all.d8nd8f [all.d8nd8f $table$logFC > 0,], n= 20L)
write.csv(up.d8nd8f, "d8n vs d8f.csv")

qlf.d6fd8m <- glmQLFTest(fit, contrast = my.contrasts[,"D6_FvsD8_M"])
all.d6fd8m  <- topTags(qlf.d6fd8m , n = Inf)
up.d6fd8m <- head(all.d6fd8m [all.d6fd8m $table$logFC > 0,], n= 20L)
write.csv(up.d6fd8m , "d6f vs d8n.csv")

qlf.d6fd8f <- glmQLFTest(fit, contrast = my.contrasts[,"D6_FvsD8_F"])
all.d6fd8f <- topTags(qlf.d6fd8f, n = Inf)
up.d6fd8f<- head(all.d6fd8f[all.d6fd8f$table$logFC > 0,], n= 20L)
write.csv(up.d6fd8f, "d6f vs d8f.csv")

################################################################################
## Make a list of all results
################################################################################
DEresult.list <- list(all.d1d6f, all.d1d6n, all.d1d8f,
                      all.d1d8m, all.d1d8n, all.d6fd6n,
                      all.d6fd8f, all.d6fd8m, all.d8nd8f,
                      all.d8nd8m)
names(DEresult.list) <- c('all.d1d6f', 'all.d1d6n', 'all.d1d8f',
'all.d1d8m', 'all.d1d8n', 'all.d6fd6n',
'all.d6fd8f', 'all.d6fd8m', 'all.d8nd8f',
'all.d8nd8m')
file.name <- c('all.d1d6f', 'all.d1d6n', 'all.d1d8f',
               'all.d1d8m', 'all.d1d8n', 'all.d6fd6n',
               'all.d6fd8f', 'all.d6fd8m', 'all.d8nd8f',
               'all.d8nd8m')

################################################################################
### Save into files in a new folder
################################################################################
dir.create(paste0("all results")) ### create a new folder
for (i in 1:length(DEresult.list)) {
  
  write.csv(DEresult.list[[i]], 
            paste0("C:/Users/pimonrat.t/OneDrive - Chiang Mai University/Research project/Suparat/RNA/all results/",
                   file.name[i],
                   ".csv"))
} 

################################################################################
########## volcano plot
### https://biostatsquid.com/volcano-plots-r-tutorial/
### cutoff is logFC > 1 and Pvalue < 0.05
### cutoff can be changed according to what you want.
################################################################################
head(all.d1d6f)

all.d1d6f.df <- as.data.frame(all.d1d6f)
all.d1d6f.df$diffexpressed <- "No"
all.d1d6f.df$diffexpressed[all.d1d6f.df$logFC > 1 & 
                             all.d1d6f.df$PValue < 0.05] <- "Up"
all.d1d6f.df$diffexpressed[all.d1d6f.df$logFC < -1 & 
                             all.d1d6f.df$PValue < 0.05] <- "Down"
all.d1d6f.df$diffexpressed <- as.factor(all.d1d6f.df$diffexpressed)

################################################################################
#### create a new list of df for volcano plots
### and put some more columns for label purposes
################################################################################
volcano.dflist <- list()
for (i in 1:10) {
  volcano.dflist[[i]] <- as.data.frame(DEresult.list[[i]])
  volcano.dflist[[i]]$diffexpressed <- "No"
  volcano.dflist[[i]]$diffexpressed[volcano.dflist[[i]]$logFC > 1 & 
                                      volcano.dflist[[i]]$PValue < 0.05] <- "Up"
  volcano.dflist[[i]]$diffexpressed[volcano.dflist[[i]]$logFC < -1 & 
                                      volcano.dflist[[i]]$PValue < 0.05] <- "Down"
  volcano.dflist[[i]]$diffexpressed <- as.factor(volcano.dflist[[i]]$diffexpressed)
  volcano.dflist[[i]]$delabel <- ifelse(volcano.dflist[[i]]$genes 
                                 %in% head(volcano.dflist[[i]][order(volcano.dflist[[i]]$PValue), 
                                                        "genes"], 30), volcano.dflist[[i]]$genes, NA)
  # last line create label of top 30
}

names(volcano.dflist) <- c('all.d1d6f', 'all.d1d6n', 'all.d1d8f',
                          'all.d1d8m', 'all.d1d8n', 'all.d6fd6n',
                          'all.d6fd8f', 'all.d6fd8m', 'all.d8nd8f',
                          'all.d8nd8m')

################################################################################
#### draw plots in to a list and save in to files.
################################################################################
install.packages("ggplot2")
library(ggplot2)
library(ggrepel)
#BiocManager::install('EnhancedVolcano') not in use for now
#library(EnhancedVolcano)

volcano.plotlist <- list()
for (i in 1:10) {
volcano.plotlist[[i]] <- 
  ggplot(data = volcano.dflist[[i]], 
         aes(x = logFC, y = -log10(PValue), col = diffexpressed, 
             label = delabel)) + 
    labs(color = 'Regulation', #legend_title, 
         x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
    geom_point() + theme_minimal() +
    ggtitle(names(volcano.dflist[i])) +
    geom_text_repel(max.overlaps = Inf, colour = "black") + # To show all labels
  theme_classic() +
  theme(axis.line.x = element_line(colour = 'black', linewidth =0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', linewidth =0.5, linetype='solid')) +
  scale_colour_manual(values = c("darkgreen", "grey", "red")) # to set the colours of our variable<br /><br /><br />
                      #labels = c("Downregulated", "Not significant", "Upregulated"))
  
}

names(volcano.plotlist) <- c('all.d1d6f', 'all.d1d6n', 'all.d1d8f',
                           'all.d1d8m', 'all.d1d8n', 'all.d6fd6n',
                           'all.d6fd8f', 'all.d6fd8m', 'all.d8nd8f',
                           'all.d8nd8m')
################################################################################
#### write plots to a file
#### change working directory accordingly
################################################################################
setwd("C:/Users/pimonrat.t/OneDrive - Chiang Mai University/Research project/Suparat/RNA/Volcano plots")
for(i in seq_along(volcano.plotlist)){
  ggsave(filename = paste0(names(volcano.dflist[i]), ".png"),
  plot = volcano.plotlist[[i]],
  width = 11,
  height = 7,
  units = "in")
}
################################################################################
##### Making heatmaps
##### data are in an object called d (DGEList already)
### Heatmaps are a popular way to display differential 
### expression results for publication purposes. To create a heatmap, 
### we first convert the read counts into log2-counts-per-million (logCPM) 
### values. This can be done with the cpm function:
################################################################################
logCPM1 <- cpm(d, prior.count = 2, log = TRUE) #prior count = 2 to avoid log 0
rownames(logCPM1) <- d$genes$genes 
colnames(logCPM1) <- paste(d$samples$group, 1:2, sep="-")
head(logCPM1)
## group = Levels: D1 D6_F D6_N D8_F D8_M D8_N
## 
#tr2 <- glmTreat(fit, contrast = my.contrasts, 
#               lfc=log2(1.5))

o <- order(tr$table$PValue)
logCPM1 <- logCPM1[o[1:30],]

logCPM1 <- t(scale(t(logCPM1)))

install.packages("gplots")
library(gplots)
col.pan <- colorpanel(100, "blue", "white", "red")
heatmap.2(logCPM1, col = col.pan, Rowv = TRUE, scale = "none",
          trace = "none", dendrogram = "both", cexRow = 1, cexCol = 1.4,
          density.info="none",
          margin = c(7,10), lhei = c(2,10), lwid = c(2,6),
          main = "Heat map across all samples using the top 30 
          most DE genes between the D1 and the D6_F group")
par(cex.main = 1) 

################################################################################

