library(tidyverse)
library(oligo)
library(RColorBrewer)
library(limma)
library(annotate)
library(gplots)
# library(ragene10stprobeset.db)
library(pd.ragene.1.0.st.v1)
library(affycoretools)

#-----------------#
#### Read data ####
#-----------------#

setwd("/home/svillicaña/INMEGEN_LB")

# Directories
# rds.dir <- "/home/svillicaña/INMEGEN_LB/resultados/RDS/expression/"
# plot.dir <- "/home/svillicaña/INMEGEN_LB/resultados/Plots/expression/"
# csv.dir <- "/home/svillicaña/INMEGEN_LB/resultados/CSV/expression/"

rds.dir <- "/home/svillicaña/INMEGEN_LB/resultados/RDS/expression/"
plot.dir <- "/home/jpeña/INMEGEN_LB/resultados/Plots/expression/"
csv.dir <- "/home/jpeña/INMEGEN_LB/resultados/CSV/expression/"

# Sample sheet
celfiles_raw <- read.table('data/sindrome_metabolico_group_file.txt',
                       header = T, as.is = T)

# Group technical replicates
celfiles <- celfiles_raw %>%
  mutate(group = ifelse(str_detect(sample, "_RT$"), "RT", group) )

# List of groups
contr.lab <- unique(celfiles$group)
groups_list <- lapply(as.list(contr.lab),
                      function(x) which(celfiles$group == x))
names(groups_list) <- contr.lab
palette <- factor(celfiles$group, 
                  labels = brewer.pal(n = n_distinct(celfiles$group), name = "Dark2")) %>%
  as.character()

# Load normalized data
data.RMA <- readRDS(paste0(rds.dir, "data.RMA.rds"))
expr.mat <- exprs(data.RMA)

#-------------------------------#
#### Differential expression ####
#-------------------------------#

## Design matrix
n_groups <- length(groups_list)
n_array <- lengths(groups_list) %>% sum()
design <- matrix(rep(0, n_groups*n_array), nrow = n_array)
colnames(design) <- contr.lab
design[groups_list$hipotalamo, 1] <- 1
design[groups_list$pulmon, 2] <- 1
design[groups_list$stem_cells, 3] <- 1
design[groups_list$tejido_adiposo, 4] <- 1
design[groups_list$RT, 5] <- 1

## Contrasts
# Contrast description
cd.1 <- "hipotalamo vs. pulmon"
cd.2 <- "hipotalamo vs. stem_cells"
cd.3 <- "hipotalamo vs. tejido_adiposo"
cd.4 <- "pulmon vs. stem_cells"
cd.5 <- "pulmon vs. tejido_adiposo"
cd.6 <- "stem_cells vs. tejido_adiposo"
# Contrast matrix
contrasts <- makeContrasts("hipotalamo - pulmon",
                           "hipotalamo - stem_cells",
                           "hipotalamo - tejido_adiposo",
                           "pulmon - stem_cells",
                           "pulmon - tejido_adiposo",
                           "stem_cells - tejido_adiposo",
                           levels = design)

## Model fitting
# annDB <- "ragene10stprobeset.db"
p.thresh <- 0.05 # p-value threshold
lfc.thresh <- 0.3 # Log fold change threshold
n_transcripts <- 29214 # Maximum given by topTable
fit <- lmFit(expr.mat, design)
fitC <- contrasts.fit(fit, contrasts)
fitCB <- eBayes(fitC)

# Load annotation information
all.eset <- annotateEset(data.RMA, pd.ragene.1.0.st.v1)
f.data <- fData(all.eset)
f.data <- f.data[, -1]  # Remove PROBESETID column
colnames(f.data) <- c("ID", "Symbol", "Gene.Name")

## Differentially expressed genes
diff.contrast <- function(fit.mod, cont.ind, f.name, max.n, p.thresh, lfc.thresh) {
  
  # Sort by ranking
  TT <- topTable(fit = fit.mod, coef = cont.ind, adjust="fdr", sort.by="logFC", number = max.n)
  
  # Filter transcripts according to the lowest thresholds
  selected <- TT[TT$P.Value <= p.thresh & abs(TT$logFC) >= lfc.thresh, ]
  
  # Annotate the probesets
  probe.labs <- rownames(selected)
  anno.info <- f.data[probe.labs, ]
  diff.expr <- data.frame(anno.info, selected)
  indx <- !is.na(diff.expr$Symbol)
  diff.expr <- diff.expr[indx, ] # Only get those with a gene symbol attached
  
  ## Get genes below specified p-value and above specified lfc threshold
  diff.genes <- diff.expr[diff.expr$P.Value < p.thresh & abs(diff.expr$logFC) > lfc.thresh, ]
  
  # Change p values and adjusted p values to scientific notation
  diff.genes$P.Value <- format(diff.genes$P.Value, scientific = TRUE)
  diff.genes$adj.P.Val <- format(diff.genes$adj.P.Val, scientific = TRUE)
  
  # Write results (only genes in general)
  write.csv(diff.genes, file = paste0(csv.dir, f.name, ".csv"))

  # Return differentially expressed genes table
  return(diff.genes)
}

# Get the differential expressed genes for each contrast
diff.list <- c()  # Vector to store differentially expressed data frames

contr.names <- colnames(contrasts) %>% str_replace_all(" ", "")
for (contr.idx in 1:length(contr.names)){
  
  diff.genes <- diff.contrast(fit.mod = fitCB, cont.ind = contr.idx, 
                              f.name = contr.names[contr.idx], max.n = n_transcripts,
                              p.thresh = p.thresh, lfc.thresh = lfc.thresh)
  
  diff.list <- c(diff.list, diff.genes)
}

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

# Remove probes with no gene symbol
probe.labs <- rownames(fitCB)
syms <- f.data[probe.labs, "Symbol"]
indx <- !is.na(syms)
fitCBS <- fitCB[indx, ]

# Volcano plots for all contrasts
cd.list <- c(cd.1, cd.2, cd.3, cd.4, cd.5, cd.6)
for (contr.idx in 1:length(contr.names)){
  
  # Save the plot as PDF
  pdf(paste0(plot.dir, "volcano_", contr.names[contr.idx], ".pdf"), width = 14, height = 8)
  contrast.volcano(fit.mod = fitCBS, cont.ind = contr.idx,
                   title = paste0("Differential expression: ", cd.list[contr.idx]),
                   p.thresh = p.thresh, lfc.thresh = lfc.thresh)
  dev.off()
}

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
  gene.syms <- diff.expr$Symbol
  
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

# Heatmaps for all contrasts
for (contr.idx in 1:length(contr.names)){
  
  # Split contrast groups
  groups <- str_split(contr.names[contr.idx], "-")
  
  # Save the plot as PDF
  pdf(paste0(plot.dir, "heatmap_", contr.names[contr.idx], ".pdf"), width = 9, height = 16)
  contrast.heatmap(groups[1], groups[2], diff.expr = diff.list[contr.idx], 
                   title = paste0("Expression heatmap: ", cd.list[contr.idx]), 
                   f.name = paste0("expression-heatmap_", contr.names[contr.idx]))
  dev.off()
}







