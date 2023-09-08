#############################################
### ALL KCLO3 analysis
### No replicate
#############################################
library(edgeR)

### read count data
countkclo3 <- read.csv("ALL_KCLO3_RC_02.csv", header = T)
str(countkclo3)
head(countkclo3)

matrix.k <- as.matrix(countkclo3[,-1])
head(matrix.k)
colnames(matrix.k)
groups.k <- colnames(matrix.k)

groups.k <- as.factor(groups.k)
levels(groups.k)

genes.k <- (countkclo3$gene_id)
length(genes.k)

library(edgeR)
dk <- DGEList(counts = matrix.k, group = groups.k, 
              genes = genes.k)
dk

## keep the old one in case we mess up
dk.full <- dk
dk.full$samples
################################################################################
### filtering data
################################################################################
dim(dk)

head(dk$counts)

head(cpm(dk))

apply(dk.full$counts, 2, sum)

keep <- rowSums(cpm(dk.full) > 0.5) >= 1
table(keep)

dk <- dk.full[keep, , keep.lib.sizes=FALSE]
dim(dk)

### reset the library sizes
dk$samples$lib.size <- colSums(dk$counts)
dk$samples

################################################################################
####Normalizing the data
################################################################################
dk <- calcNormFactors(dk)
dk

################################################################################
#### Data Exploration
## Before proceeding with the computations for differential expression, 
## it is possible to produce a plot showing the sample relations 
## based on multidimensional scaling.
################################################################################
colors <- rep(seq(1:15),2)
plotMDS(dk, col = colors[groups.k])

plotMD(dk, column = 2)
abline(h=0, col="red", lty=2, lwd=2)

##############################################################
### No replication 
### Options 2 in topic  2.12 EdgeR manual
### The results will be sensitive to bcv 
############################################################
bcv <- 0.1 

et1 <- exactTest(dk, dispersion = bcv^2, pair=c(1,16))
et2 <- exactTest(dk, dispersion = bcv^2, pair=c(2,17))
et3 <- exactTest(dk, dispersion = bcv^2, pair=c(3,18))
et4 <- exactTest(dk, dispersion = bcv^2, pair=c(4,19))
et5 <- exactTest(dk, dispersion = bcv^2, pair=c(5,20))
et6 <- exactTest(dk, dispersion = bcv^2, pair=c(6,21))
et7 <- exactTest(dk, dispersion = bcv^2, pair=c(7,22))
et8 <- exactTest(dk, dispersion = bcv^2, pair=c(8,23))
et9 <- exactTest(dk, dispersion = bcv^2, pair=c(9,24))
et10 <- exactTest(dk, dispersion = bcv^2, pair=c(10,25))
et11 <- exactTest(dk, dispersion = bcv^2, pair=c(11,26))
et12 <- exactTest(dk, dispersion = bcv^2, pair=c(12,27))
et13 <- exactTest(dk, dispersion = bcv^2, pair=c(13,28))
et14 <- exactTest(dk, dispersion = bcv^2, pair=c(14,29))
et15 <- exactTest(dk, dispersion = bcv^2, pair=c(15,30))

################################################################################
### make data frame of all comparisons
### and write to csv files
################################################################################
all.et1 <- topTags(et1, n = Inf)
all.et2 <- topTags(et2, n = Inf)
all.et3 <- topTags(et3, n = Inf)
all.et4 <- topTags(et4, n = Inf)
all.et5 <- topTags(et5, n = Inf)
all.et6 <- topTags(et6, n = Inf)
all.et7 <- topTags(et7, n = Inf)
all.et8 <- topTags(et8, n = Inf)
all.et9 <- topTags(et9, n = Inf)
all.et10 <- topTags(et10, n = Inf)
all.et11 <- topTags(et11, n = Inf)
all.et12 <- topTags(et12, n = Inf)
all.et13 <- topTags(et13, n = Inf)
all.et14 <- topTags(et14, n = Inf)
all.et15 <- topTags(et15, n = Inf)

DEresult.kclo3 <- list(all.et1, all.et2, all.et3,
                       all.et4, all.et5, all.et6,
                       all.et7, all.et8, all.et9,
                       all.et10, all.et11, all.et12,
                       all.et13, all.et14, all.et15)
names(DEresult.kclo3) <- c("all.et1", "all.et2", "all.et3",
                           "all.et4", "all.et5", "all.et6",
                           "all.et7", "all.et8", "all.et9",
                           "all.et10", "all.et11", "all.et12",
                           "all.et13", "all.et14", "all.et15")
file.name <- c("all.et1", "all.et2", "all.et3",
               "all.et4", "all.et5", "all.et6",
               "all.et7", "all.et8", "all.et9",
               "all.et10", "all.et11", "all.et12",
               "all.et13", "all.et14", "all.et15")

################################################################################
### Save into files in a new folder
################################################################################
dir.create(paste0("all kclo3 results June21")) ### create a new folder
for (i in 1:length(DEresult.kclo3)) {
  
  write.csv(DEresult.kclo3[[i]], 
            paste0("C:/Users/pimonrat.t/OneDrive - Chiang Mai University/Research project/Suparat/RNA/all kclo3 results June21/",
                   file.name[i],
                   ".csv"))
} 

################################################################################
#### create a new list of df for volcano plots
### and put some more columns for label purposes
################################################################################
volcano.dflist <- list()
for (i in 1:15) {
  volcano.dflist[[i]] <- as.data.frame(DEresult.kclo3[[i]])
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

names(volcano.dflist) <- c("all.et1", "all.et2", "all.et3",
                           "all.et4", "all.et5", "all.et6",
                           "all.et7", "all.et8", "all.et9",
                           "all.et10", "all.et11", "all.et12",
                           "all.et13", "all.et14", "all.et15")

################################################################################
#### draw plots in to a list and save in to files.
################################################################################
install.packages("ggplot2")
library(ggplot2)
library(ggrepel)
#BiocManager::install('EnhancedVolcano') not in use for now
#library(EnhancedVolcano)

volcano.plotlist <- list()
for (i in 1:15) {
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

names(volcano.plotlist) <- c("all.et1", "all.et2", "all.et3",
                             "all.et4", "all.et5", "all.et6",
                             "all.et7", "all.et8", "all.et9",
                             "all.et10", "all.et11", "all.et12",
                             "all.et13", "all.et14", "all.et15")

################################################################################
#### write plots to a file
#### change working directory accordingly
################################################################################
dir.create(paste0("all kclo3 volcano plots June21")) 
setwd("C:/Users/pimonrat.t/OneDrive - Chiang Mai University/Research project/Suparat/RNA/all kclo3 volcano plots June21")
for(i in seq_along(volcano.plotlist)){
  ggsave(filename = paste0(names(volcano.dflist[i]), ".png"),
         plot = volcano.plotlist[[i]],
         width = 11,
         height = 7,
         units = "in")
}

