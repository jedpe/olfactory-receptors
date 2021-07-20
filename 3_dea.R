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
rds.dir <- "/home/svillicaña/INMEGEN_LB/resultados/RDS/expression/"
plot.dir <- "/home/svillicaña/INMEGEN_LB/resultados/Plots/expression/"
csv.dir <- "/home/svillicaña/INMEGEN_LB/resultados/CSV/expression/"

# rds.dir <- "/home/svillicaña/INMEGEN_LB/resultados/RDS/expression/"
# plot.dir <- "/home/jpeña/INMEGEN_LB/resultados/Plots/expression/"
# csv.dir <- "/home/jpeña/INMEGEN_LB/resultados/CSV/expression/"

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
  
  # Annotate the probesets
  probe.labs <- rownames(TT)
  anno.info <- f.data[probe.labs, ]
  diff.expr <- data.frame(anno.info, TT)
  indx <- !is.na(diff.expr$Symbol)
  diff.genes <- diff.expr[indx, ] # Only get those with a gene symbol attached
  
  # Write results
  write.csv(diff.genes, file = paste0(csv.dir, f.name, ".csv"))

  # Return differentially expressed genes table
  return(diff.genes)
}

# Get the differential expressed genes for each contrast
diff.list <- list()  # Vector to store differentially expressed data frames

contr.names <- colnames(contrasts) %>% str_replace_all(" ", "")
for (contr.idx in 1:length(contr.names)){
  
  diff.genes <- diff.contrast(fit.mod = fitCB, cont.ind = contr.idx, 
                              f.name = contr.names[contr.idx], max.n = n_transcripts,
                              p.thresh = p.thresh, lfc.thresh = lfc.thresh)
  
  diff.list[[contr.idx]] <- diff.genes
}

# Filter Orl
diff.list.orl <- mapply(function(diff.genes, f.name){
  diff.genes.orl <- diff.genes %>%
    dplyr::filter(str_detect(Gene.Name, "olfactory receptor")) %>%
    mutate(adj.P.Val = p.adjust(P.Value, method = "fdr"))
  
  # Write results
  write.csv(diff.genes.orl, file = paste0(csv.dir, f.name, "_olr_genes.csv"))
}, diff.list, contr.names)

#---------------------#
#### Volcano plots ####
#---------------------#
  
contrast.volcano <- function(diff.genes, title, p.thresh, lfc.thresh) {
  log.p.thresh <- -log10(p.thresh)
  
  # Prepare data
  diff.genes <- diff.genes %>%
    arrange(adj.P.Val) %>%
    transmute(Symbol, logFC, log.adj.P.Val = -log10(adj.P.Val))
  diff.genes$Group <- "NA"
  diff.genes$Group[diff.genes$logFC >= lfc.thresh &
                     diff.genes$log.adj.P.Val >= log.p.thresh] <- "UP"
  diff.genes$Group[diff.genes$logFC <= -lfc.thresh &
                     diff.genes$log.adj.P.Val >= log.p.thresh] <- "DOWN"
  
  # Plot
  p <- ggplot(diff.genes, aes(x = logFC, y = log.adj.P.Val)) +
    geom_point(aes(colour = Group)) +
    ggtitle(title) +
    xlab("Log2 Fold Change") + 
    ylab("-log(p-value)") +
    scale_colour_manual(values = c("DOWN" = "firebrick", "NA" = "darkgrey", "UP" = "seagreen")) +
    geom_hline(yintercept = log.p.thresh, linetype="dashed") +
    geom_vline(xintercept = lfc.thresh) + 
    geom_vline(xintercept = -lfc.thresh) + 
    theme(legend.position = "none")
  
  return(p)
}

# Volcano plots for all contrasts
cd.list <- c(cd.1, cd.2, cd.3, cd.4, cd.5, cd.6)
for (contr.idx in 1:length(contr.names)){
  # Save the plot as PDF
  volcano.idx <- contrast.volcano(diff.list[[contr.idx]],
                                  title = paste0("Differential expression: ", cd.list[contr.idx]),
                                  p.thresh = p.thresh, lfc.thresh = lfc.thresh)
  ggsave(paste0(plot.dir, "volcano_", contr.names[contr.idx], ".pdf"),
         volcano.idx, width = 14, height = 8, device = "pdf")
  
  # Orl
  volcano.idx <- contrast.volcano(diff.list.orl[[contr.idx]],
                                  title = paste0("Differential expression of Orl genes: ", cd.list[contr.idx]),
                                  p.thresh = p.thresh, lfc.thresh = lfc.thresh)
  ggsave(paste0(plot.dir, "volcano_", contr.names[contr.idx], "_olr_genes.pdf"),
         volcano.idx, width = 14, height = 8, device = "pdf")
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







