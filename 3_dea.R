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
html.dir <- "/home/svillicaña/INMEGEN_LB/resultados/HTML/expression/"
csv.dir <- "/home/svillicaña/INMEGEN_LB/resultados/CSV/expression/"

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
#paste0(html.dir, f.name, ".html")
cd.1 <- "hipotalamo vs. pulmon"
cd.2 <- "hipotalamo vs. stem_cells"
cd.3 <- "hipotalamo vs. tejido_adiposo"

# First contrast (hipotalamo vs. pulmon)
diff.list1 <- diff.contrast(fit.mod = fitCB, cont.ind = 1, f.name = "hipotalamo-pulmon",
                            max.n = n_transcripts, 
                            html.title = paste0("Differetial expression: ", cd.1),
                            p.thresh = p.thresh, lfc.thresh = lfc.thresh)

# Second contrast (hipotalamo vs. stem cells)

# Third contrast (hipotalamo vs. tejido adiposo)
















