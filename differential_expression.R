#************************************************************************************
# Data groups:
#             1) CT (28_CT, 29_CT, 33_CT, 38_CT, 39_CT, 42_CT, 46_CT, 47_CT, 48_CT)
#             2) AO (5_AO, 6_AO, 12_AO, 16_AO, 18_AO, 19_AO, 21_AO, 22_AO, 23_AO)
#             3) EO (1_EO, 2_EO, 3_EO, 10_EO, 11_EO, 20_EO, 24_EO, 26_EO, 27_EO)
#
# Contrasts:
#             1) AO vs. CT
#             2) EO vs. CT
#             3) EO vs. AO
#
# Total 27 microarrays
#************************************************************************************

library(tidyverse)
library(oligo)

#-----------------#
#### Read data ####
#-----------------#

# Get the raw data
setwd("~/Documents/INMEGEN/Ataxia/")
# setwd("~/Projects/Ataxia/")
celfiles <- list.celfiles('./celfiles/', full.names = T)
raw.data <- read.celfiles(celfiles)

## Label the samples accordingly
# Reorder the samples
ct.idx <- c(16, 17, 19:25)
ao.idx <- c(4:7, 10:12, 26, 27)
eo.idx <- c(1:3, 8, 9, 13:15, 18)
raw.data <- raw.data[, c(ct.idx, ao.idx, eo.idx)]

# Labels
contr.lab <- c("CT", "AO", "EO")
sampleNames(raw.data) <- sampleNames(raw.data) %>%
  str_replace(., "\\.CEL", "")
palette <- c(rep("steelblue", 9), rep("violet", 9), rep("orangered", 9))

# Directories
plot.dir <- "~/Documents/INMEGEN/Ataxia/resultados/Plots/expression/"
html.dir <- "~/Documents/INMEGEN/Ataxia/resultados/HTML/expression/"
csv.dir <- "~/Documents/INMEGEN/Ataxia/resultados/CSV/expression/"

# plot.dir <- "~/Projects/Ataxia/resultados/Plots/expression/"
# html.dir <- "~/Projects/Ataxia/resultados/HTML/expression/"
# csv.dir <- "~/Projects/Ataxia/resultados/CSV/expression/"

#----------------#
#### QC plots ####
#----------------#

## Data before the normalization
# target="core" is used to visualize data at transcript level
pdf(paste0(plot.dir, "raw_data_boxplot.pdf"), width = 9, height = 16)
oligo::boxplot(raw.data, target = "core", col = palette,
               main = "Raw data distribution", cex = 0.8, las = 2)
dev.off()

pdf(paste0(plot.dir, "raw_data_density.pdf"), width = 9, height = 16)
oligo::hist(raw.data, target = "core", col = palette,
            main = "Raw data density curves", cex = 0.8, lwd = 2)
dev.off()

#---------------------#
#### Normalization ####
#---------------------#

## Normalization by quantiles
groups_list <- list(CT = 1:9, AO = 10:18, EO = 19:27)
norm.m <- "quantile"
data.norm <- raw.data

q.norm <- normalize(raw.data[, groups_list$CT], method = norm.m)
pm(data.norm)[, groups_list$CT] <- pm(q.norm)

q.norm <- normalize(raw.data[, groups_list$AO], method = norm.m)
pm(data.norm)[, groups_list$AO] <- pm(q.norm)

q.norm <- normalize(raw.data[, groups_list$EO], method = norm.m)
pm(data.norm)[, groups_list$EO] <- pm(q.norm)

# Plots after quantiles normalization
pdf(paste0(plot.dir, "qnorm_boxplot.pdf"), width = 9, height = 16)
oligo::boxplot(data.norm, target = "core", col = palette,
               main = "Quantile-normalized by group", cex = 0.8, las = 2)
dev.off()

pdf(paste0(plot.dir, "qnorm_density.pdf"), width = 9, height = 16)
oligo::hist(data.norm, target = "core", col = palette,
            main = "Quantile-normalized by group", cex = 0.8, lwd = 2)
dev.off()

## Normalization by RMA
data.RMA <- rma(data.norm, target = "core")
expr.mat <- exprs(data.RMA)

# Plots after RMA normalization
pdf(paste0(plot.dir, "RMA_boxplot.pdf"), width = 9, height = 16)
oligo::boxplot(data.RMA, col = palette, main = "RMA-normalized data (all)",
               las = 2, ylab = "Expression levels")
dev.off()

pdf(paste0(plot.dir, "RMA_density.pdf"), width = 9, height = 16)
oligo::hist(data.RMA, col = palette, main = "RMA-normalized data (all)",
            cex = 0.8, lwd = 2)
dev.off()

#-------------------------------#
#### Differential expression ####
#-------------------------------#
library(limma)
library(annotate)
library(gplots)
library(clariomdhumantranscriptcluster.db)

## Design matrix
n_groups <- length(groups_list)
n_array <- lengths(groups_list) %>% sum()
design <- matrix(rep(0, n_groups*n_array), nrow = n_array)
colnames(design) <- contr.lab
design[groups_list$CT, 1] <- 1
design[groups_list$AO, 2] <- 1
design[groups_list$EO, 3] <- 1

## Contrasts
# Contrast description
cd.1 <- "AO vs. CT"
cd.2 <- "EO vs. CT"
cd.3 <- "EO vs. AO"
# Contrast matrix
contrasts <- makeContrasts("AO - CT", "EO - CT", "EO - AO", levels = design)

## Model fitting
annDB <- "clariomdhumantranscriptcluster.db"
# b.thresh <- -3.3 # B statistic threshold (genes)
# mirna.b <- -3.6 # B statistic threshold (miRNA)
# lncrna.b <- -4 # B statistic threshold (lncRNA)
p.thresh <- 0.01 # p-value threshold (genes)
mirna.p <- 0.01 # p-value threshold (miRNA)
lncrna.p <- 0.05 # p-value threshold (lncRNA)
lfc.thresh <- 0.5 # Log fold change threshold (genes)
mirna.lfc <- 0.45 # Log fold change threshold (miRNA)
lncrna.lfc <- 0.3 # Log fold change threshold (lncRNA)
lncrna.lfc.2 <- 0.37 # Log fold change threshold (lncRNA) - first contrast only
n_transcripts <- 542500 # According to the data sheet from Thermofisher
fit <- lmFit(expr.mat, design)
fitC <- contrasts.fit(fit, contrasts)
fitCB <- eBayes(fitC)

## Differential expressed genes
diff.contrast <- function(fit.mod, cont.ind, f.name, html.title, max.n, p.thresh, lfc.thresh, 
                          mirna.p, mirna.lfc, lncrna.p, lncrna.lfc) {
  # Sort by ranking
  TT <- topTable(fit = fit.mod, coef = cont.ind, adjust="fdr", sort.by="logFC", number = max.n)
  
  # Get the thresholds
  p.val <- max(p.thresh, mirna.p, lncrna.p) # We want the maximum to avoid discarding entries (< p.val)
  lfc.val <- min(lfc.thresh, mirna.lfc, lncrna.lfc) # We don't want to discard entries (> lfc.val)
  
  # Filter transcripts according to the lowest thresholds
  selected <- TT[TT$P.Value < p.val & abs(TT$logFC) > lfc.val, ]
  
  # Annotate the transcripts
  probe.labs <- rownames(selected)
  sym <- getSYMBOL(probe.labs, data = annDB)
  diff.expr <- data.frame(sym, selected)
  indx <- !is.na(diff.expr$sym)
  diff.expr <- diff.expr[indx, ] # Only get thoese with a gene symbol attached
  
  ## Get genes, miRNA and lncRNA with their corresponding thresholds
  # Genes in general
  diff.genes <- diff.expr[diff.expr$P.Value < p.thresh & abs(diff.expr$logFC) > lfc.thresh, ]
  # Using miRNA thresholds (DOES NOT RETURN miRNA ONLY ENTIRES)
  diff.mirna <- diff.expr[diff.expr$P.Value < mirna.p & abs(diff.expr$logFC) > mirna.lfc, ]
  # Using lncRNA thresholds (DOES NOT RETURN lncRNA ONLY ENTIRES)
  diff.lncrna <- diff.expr[diff.expr$P.Value < lncrna.p & abs(diff.expr$logFC) > lncrna.lfc, ]
  
  # Change p values and adjusted p values to scientific notation
  diff.genes$P.Value <- format(diff.genes$P.Value, scientific = TRUE)
  diff.mirna$P.Value <- format(diff.mirna$P.Value, scientific = TRUE)
  diff.lncrna$P.Value <- format(diff.lncrna$P.Value, scientific = TRUE)
  
  diff.genes$adj.P.Val <- format(diff.genes$adj.P.Val, scientific = TRUE)
  diff.mirna$adj.P.Val <- format(diff.mirna$adj.P.Val, scientific = TRUE)
  diff.lncrna$adj.P.Val <- format(diff.lncrna$adj.P.Val, scientific = TRUE)
  
  # Write results (only genes in general)
  write.csv(diff.genes, file = paste0(csv.dir, f.name, ".csv"))
  genes.entID <- rownames(diff.genes) %>% getEG(., data = annDB)
  htmlpage(genelist = list(genes.entID), filename = paste0(html.dir, f.name, ".html"), 
           title = html.title, 
           othernames = diff.genes, table.head = c("Entrez ID", colnames(diff.genes)), 
           table.center = TRUE)
  
  # Return differentially expressed matrices for later use
  list(genes = diff.genes, mirna = diff.mirna, lncrna = diff.lncrna) %>% return()
}

# First contrast (AO - CT)
diff.list1 <- diff.contrast(fit.mod = fitCB, cont.ind = 1, f.name = "AO-CT", max.n = n_transcripts, 
                            html.title = paste0("Differential expression: ", cd.1), 
                            p.thresh = p.thresh, lfc.thresh = lfc.thresh, 
                            mirna.p = mirna.p, mirna.lfc = mirna.lfc, 
                            lncrna.p = lncrna.p, lncrna.lfc = lncrna.lfc)
diff.expr1 <- diff.list1$genes

# Second contrast (EO - CT)
diff.list2 <- diff.contrast(fit.mod = fitCB, cont.ind = 2, f.name = "EO-CT", max.n = n_transcripts, 
                            html.title = paste0("Differential expression: ", cd.2), 
                            p.thresh = p.thresh, lfc.thresh = lfc.thresh, 
                            mirna.p = mirna.p, mirna.lfc = mirna.lfc, 
                            lncrna.p = lncrna.p, lncrna.lfc = lncrna.lfc)
diff.expr2 <- diff.list2$genes

# Third contrast (EO - AO)
diff.list3 <- diff.contrast(fit.mod = fitCB, cont.ind = 3, f.name = "EO-AO", max.n = n_transcripts, 
                            html.title = paste0("Differential expression: ", cd.3), 
                            p.thresh = p.thresh, lfc.thresh = lfc.thresh, 
                            mirna.p = mirna.p, mirna.lfc = mirna.lfc, 
                            lncrna.p = lncrna.p, lncrna.lfc = lncrna.lfc)
diff.expr3 <- diff.list3$genes

#---------------------#
#### Volcano plots ####
#---------------------#

contrast.volcano <- function(fit.mod, title, cont.ind, pch.st = 20, p.thresh, lfc.thresh) {
  # Transform p-value threshold to -log(p-value)
  logP.thresh <- -log10(p.thresh) # The smaller a p-value gets, the bigger the -log(p-value)
  
  # Matrices of coefficients and logP
  coef.mat <- as.matrix(fit.mod$coef) # logFC
  logP.mat <- -log10(fit.mod$p.value) %>% as.matrix() # -log(p-value)
  
  # Set the x and y-axis
  xl <- min(coef.mat[, cont.ind]) - 0.5
  xu <- max(coef.mat[, cont.ind]) + 0.5
  yl <- 0
  yu <- max(logP.mat[, cont.ind]) + 0.5
  
  # Volcano base plot
  volcanoplot(fit.mod, col="blue", ylim=c(yl, yu), xlim=c(xl, xu), coef = cont.ind,
              main = title, cex.lab = 1.3, ylab = "-log(p-value)", style = "p-value", pch = pch.st)
  
  # Section the plot
  par(new=T) # So the next plotted thing overlaps
  abline(v=-lfc.thresh, col="brown", ylab="", xlab="")
  par(new=T)
  abline(v=lfc.thresh, col="brown", ylab="", xlab="")
  par(new=T)
  abline(h=logP.thresh, col="black", ylab="", xlab="")
  
  # Obtain indices of interest
  ind1 = abs(coef.mat[, cont.ind]) > lfc.thresh # Values above log fold change threshold
  ind2 = logP.mat[, cont.ind] > logP.thresh # Values above p-value threshold
  ind3 = (coef.mat[, cont.ind] > lfc.thresh & logP.mat[, cont.ind] > logP.thresh) # Upper right quadrant 
  ind4 = (coef.mat[, cont.ind] < -lfc.thresh & logP.mat[, cont.ind] > logP.thresh) # Upper left quadrant
  
  # Values above log fold change threshold
  x = coef.mat[ind1, cont.ind]
  y = logP.mat[ind1, cont.ind]
  par(new=T)
  plot(x, y, col="magenta",ylim=c(yl,yu), xlim=c(xl,xu),main="", pch = 20, 
       xlab="", ylab="",cex.lab=1.3, cex = 0.35)
  # Values above B threshold
  x = coef.mat[ind2, cont.ind]
  y = logP.mat[ind2, cont.ind]
  par(new=T)
  plot(x, y, col="orange",  ylim=c(yl,yu), xlim=c(xl,xu), main="", pch = 20, 
       xlab="", ylab="",cex.lab=1.3, cex = 0.35)
  # Upper right quadrant
  x = coef.mat[ind3, cont.ind]
  y = logP.mat[ind3, cont.ind]
  par(new=T)
  plot(x, y, col="red",  ylim=c(yl,yu), xlim=c(xl,xu), main="", pch = 20, 
       xlab="", ylab="",cex.lab=1.3, cex = 0.35)
  # Upper left quadrant
  x = coef.mat[ind4, cont.ind]
  y = logP.mat[ind4, cont.ind]
  par(new=T)
  plot(x, y, col="darkgreen", ylim=c(yl,yu), xlim=c(xl,xu), main="", pch = 20, 
       xlab="", ylab="",cex.lab=1.3, cex = 0.35)
}

## Remove probes with no gene symbol
# Annotate the transcripts
probe.labs <- rownames(fitCB)
sym <- getSYMBOL(probe.labs, data = annDB)
indx <- !is.na(sym)
fitCBS <- fitCB[indx, ]

# First contrast (AO - CT)
pdf(paste0(plot.dir, "volcano_AO_v_CT.pdf"), width = 14, height = 8)
contrast.volcano(fit.mod = fitCBS, cont.ind = 1, 
                 title = paste0("Differential expression: ", cd.1),
                 p.thresh = p.thresh, lfc.thresh = lfc.thresh)
dev.off()

# Second contrast (EO - CT)
pdf(paste0(plot.dir, "volcano_EO_v_CT.pdf"), width = 14, height = 8)
contrast.volcano(fit.mod = fitCBS, cont.ind = 2, 
                 title = paste0("Differential expression: ", cd.2),
                 p.thresh = p.thresh, lfc.thresh = lfc.thresh)
dev.off()

# Third contrast (EO - AO)
pdf(paste0(plot.dir, "volcano_EO_v_AO.pdf"), width = 14, height = 8)
contrast.volcano(fit.mod = fitCBS, cont.ind = 3, 
                 title = paste0("Differential expression: ", cd.3),
                 p.thresh = p.thresh, lfc.thresh = lfc.thresh)
dev.off()

#----------------#
#### Heatmaps ####
#----------------#

contrast.heatmap <- function(..., diff.expr, title, f.name) {
  # Expand the first arguments
  # substitute just create a language object, it need to be evaluated
  arrays <- substitute(list(...)) %>% eval() # This will have strings to select the columns from expr.mat
  array.regex <- paste0("(", paste0(arrays, collapse="|"), ")") # Turn it into a single regex string
  # Select corresponding columns
  cont.ind <- colnames(expr.mat) %>% grep(array.regex, .)
  # Get expression values from the differentially expressed genes
  data.clus <- expr.mat[rownames(diff.expr), cont.ind] # These will be the ones for the cluster
  # Obtain their gene symbols
  gene.syms <- diff.expr$sym
  # Set outer margin (oma) and margin (mar) arguments for the plot
  par(oma = c(3,1,3,4), mar = c(12,5,2,2)+0.1)
  # Heatmap plot 
  ind.hmap <- heatmap.2(data.clus, col=greenred(75), scale = "row", key = TRUE, symkey = FALSE, 
                        density.info = "none", trace = "none", cexRow = 0.55, cexCol = 0.7, 
                        main = title, labRow = gene.syms, ColSideColors = palette[cont.ind])
  
  # Genes 
  idx <- rev(ind.hmap$rowInd) # Reverse order of indices heatmap matrix
  genes.hmap <- cbind(data.clus[idx, ], diff.expr[idx, -1], diff.expr[idx, "sym"])
  
  # Write results
  write.csv(genes.hmap, file = paste0(csv.dir, f.name, ".csv"))
}

# First contrast (AO - CT)
pdf(paste0(plot.dir, "heatmap_AO_v_CT.pdf"), width = 9, height = 16)
contrast.heatmap("AO", "CT", diff.expr = diff.expr1, 
                 title = paste0("Expression heatmap: ", cd.1), 
                 f.name = "expression-heatmap_AO-CT")
dev.off()

# Second contrast (EO - CT)
pdf(paste0(plot.dir, "heatmap_EO_v_CT.pdf"), width = 9, height = 16)
contrast.heatmap("EO", "CT", diff.expr = diff.expr2, 
                 title = paste0("Expression heatmap: ", cd.2), 
                 f.name = "expression-heatmap_EO-CT")
dev.off()

# Third contrast (EO - AO)
pdf(paste0(plot.dir, "heatmap_EO_v_AO.pdf"), width = 9, height = 16)
contrast.heatmap("EO", "AO", diff.expr = diff.expr3, 
                 title = paste0("Expression heatmap: ", cd.3), 
                 f.name = "expression-heatmap_EO-AO")
dev.off()

### miRNA ###

# Get probes belonging to miRNAs (for the volcano plots)
mirna.sym <- base::grepl("^MIR", sym)
fitCBM <- fitCB[mirna.sym, ]

# First contrast with miRNA (AO - CT)
diff.mirna1 <- diff.list1$mirna[base::grepl("^MIR", diff.list1$mirna$sym), ]
write.csv(diff.mirna1, file = paste0(csv.dir, "miRNA-genes_AO-CT.csv"))
mirna1.entID <- rownames(diff.mirna1) %>% getEG(., data = annDB)
htmlpage(genelist = list(mirna1.entID), filename = paste0(html.dir, "mirna_AO-CT.html"), 
         title = paste0("miRNA differential expression: ", cd.1), 
         othernames = diff.mirna1, table.head = c("Entrez ID", colnames(diff.mirna1)), 
         table.center = TRUE)

pdf(paste0(plot.dir, "mirna_heatmap_AO_v_CT.pdf"), width = 9, height = 16)
contrast.heatmap("AO", "CT", diff.expr = diff.mirna1,
                 title = paste0("miRNA differential expression: ", cd.1),
                 f.name = "miRNA-heatmap_AO-CT")
dev.off()

pdf(paste0(plot.dir, "mirna_volcano_AO_v_CT.pdf"), width = 14, height = 8)
contrast.volcano(fitCBM, cont.ind = 1,
                 title = paste0("miRNA differential expression: ", cd.1),
                 p.thresh = mirna.p, lfc.thresh = mirna.lfc)
dev.off()

# Second contrast with miRNA (EO - CT)
diff.mirna2 <- diff.list2$mirna[base::grepl("^MIR", diff.list2$mirna$sym), ]
write.csv(diff.mirna2, file = paste0(csv.dir, "miRNA-genes_EO-CT.csv"))
mirna2.entID <- rownames(diff.mirna2) %>% getEG(., data = annDB)
htmlpage(genelist = list(mirna2.entID), filename = paste0(html.dir, "mirna_EO-CT.html"), 
         title = paste0("miRNA differential expression: ", cd.2), 
         othernames = diff.mirna2, table.head = c("Entrez ID", colnames(diff.mirna2)), 
         table.center = TRUE)

pdf(paste0(plot.dir, "mirna_heatmap_EO_v_CT.pdf"), width = 9, height = 16)
contrast.heatmap("EO", "CT", diff.expr = diff.mirna2,
                 title = paste0("miRNA differential expression: ", cd.2), 
                 f.name = "miRNA-heatmap_EO-CT")
dev.off()

pdf(paste0(plot.dir, "mirna_volcano_EO_v_CT.pdf"), width = 14, height = 8)
contrast.volcano(fitCBM, cont.ind = 2,
                 title = paste0("miRNA differential expression: ", cd.2),
                 p.thresh = mirna.p, lfc.thresh = mirna.lfc)
dev.off()

# Third contrast with miRNA (EO - AO)
diff.mirna3 <- diff.list3$mirna[base::grepl("^MIR", diff.list3$mirna$sym), ]
write.csv(diff.mirna3, file = paste0(csv.dir, "miRNA-genes_EO-AO.csv"))
mirna3.entID <- rownames(diff.mirna3) %>% getEG(., data = annDB)
htmlpage(genelist = list(mirna3.entID), filename = paste0(html.dir, "mirna_EO-AO.html"), 
         title = paste0("miRNA differential expression: ", cd.3), 
         othernames = diff.mirna3, table.head = c("Entrez ID", colnames(diff.mirna3)), 
         table.center = TRUE)

pdf(paste0(plot.dir, "mirna_heatmap_EO_v_AO.pdf"), width = 9, height = 16)
contrast.heatmap("EO", "AO", diff.expr = diff.mirna3,
                 title = paste0("miRNA differential expression: ", cd.3), 
                 f.name = "miRNA-heatmap_EO-AO")
dev.off()

pdf(paste0(plot.dir, "mirna_volcano_EO_v_AO.pdf"), width = 14, height = 8)
contrast.volcano(fitCBM, cont.ind = 3,
                 title = paste0("miRNA differential expression: ", cd.3),
                 p.thresh = mirna.p, lfc.thresh = mirna.lfc)
dev.off()

### lncRNA ###

# Get probes belonging to lncRNAs (for the volcano plots)
lncrna.sym <- base::grepl("^LINC", sym)
fitCBL <- fitCB[lncrna.sym, ]

# First contrast with lncRNA (AO - CT)
diff.lncrna1 <- diff.list1$lncrna[base::grepl("^LINC", diff.list1$lncrna$sym), ]
write.csv(diff.lncrna1, file = paste0(csv.dir, "lncRNA-genes_AO-CT.csv"))
lncrna1.entID <- rownames(diff.lncrna1) %>% getEG(., data = annDB)
htmlpage(genelist = list(lncrna1.entID), filename = paste0(html.dir, "lncrna_AO-CT.html"), 
         title = paste0("lncRNA differential expression: ", cd.1), 
         othernames = diff.lncrna1, table.head = c("Entrez ID", colnames(diff.lncrna1)), 
         table.center = TRUE)

pdf(paste0(plot.dir, "lncrna_heatmap_AO_v_CT.pdf"), width = 9, height = 16)
contrast.heatmap("AO", "CT", diff.expr = diff.lncrna1,
                 title = paste0("lncRNA differential expression: ", cd.1),
                 f.name = "lncRNA-heatmap_AO-CT")
dev.off()

pdf(paste0(plot.dir, "lncrna_volcano_AO_v_CT.pdf"), width = 14, height = 8)
contrast.volcano(fitCBL, cont.ind = 1,
                 title = paste0("lncRNA differential expression: ", cd.1),
                 p.thresh = lncrna.p, lfc.thresh = lncrna.lfc)
dev.off()

# Second contrast with lncRNA (EO - CT)
diff.lncrna2 <- diff.list2$lncrna[base::grepl("^LINC", diff.list2$lncrna$sym), ]
write.csv(diff.lncrna2, file = paste0(csv.dir, "lncRNA-genes_EO-CT.csv"))
lncrna2.entID <- rownames(diff.lncrna2) %>% getEG(., data = annDB)
htmlpage(genelist = list(lncrna2.entID), filename = paste0(html.dir, "lncrna_EO-CT.html"), 
         title = paste0("lncRNA differential expression: ", cd.2), 
         othernames = diff.lncrna2, table.head = c("Entrez ID", colnames(diff.lncrna2)), 
         table.center = TRUE)

pdf(paste0(plot.dir, "lncrna_heatmap_EO_v_CT.pdf"), width = 9, height = 16)
contrast.heatmap("EO", "CT", diff.expr = diff.lncrna2,
                 title = paste0("lncRNA differential expression: ", cd.2), 
                 f.name = "lncRNA-heatmap_EO-CT")
dev.off()

pdf(paste0(plot.dir, "lncrna_volcano_EO_v_CT.pdf"), width = 14, height = 8)
contrast.volcano(fitCBL, cont.ind = 2,
                 title = paste0("lncRNA differential expression: ", cd.2),
                 p.thresh = lncrna.p, lfc.thresh = lncrna.lfc)
dev.off()

# Third contrast with lncRNA (EO - AO)
diff.lncrna3 <- diff.list3$lncrna[base::grepl("^LINC", diff.list3$lncrna$sym), ]
write.csv(diff.lncrna3, file = paste0(csv.dir, "lncRNA-genes_EO-AO.csv"))
lncrna3.entID <- rownames(diff.lncrna3) %>% getEG(., data = annDB)
htmlpage(genelist = list(lncrna3.entID), filename = paste0(html.dir, "lncrna_EO-AO.html"), 
         title = paste0("lncRNA differential expression: ", cd.3), 
         othernames = diff.lncrna3, table.head = c("Entrez ID", colnames(diff.lncrna3)), 
         table.center = TRUE)

pdf(paste0(plot.dir, "lncrna_heatmap_EO_v_AO.pdf"), width = 9, height = 16)
contrast.heatmap("EO", "AO", diff.expr = diff.lncrna3,
                 title = paste0("lncRNA differential expression: ", cd.3), 
                 f.name = "lncRNA-heatmap_EO-AO")
dev.off()

pdf(paste0(plot.dir, "lncrna_volcano_EO_v_AO.pdf"), width = 14, height = 8)
contrast.volcano(fitCBL, cont.ind = 3,
                 title = paste0("lncRNA differential expression: ", cd.3),
                 p.thresh = lncrna.p, lfc.thresh = lncrna.lfc)
dev.off()























