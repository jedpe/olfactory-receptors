library(tidyverse)
library(oligo)
library(RColorBrewer)
library(limma)
library(annotate)
library(gplots)
library(ragene10stprobeset.db)

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

# Normalized data
data.RMA <- readRDS(paste0(rds.dir, "data.RMA.rds"))

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
annDB <- "ragene10stprobeset.db"
p.thresh <- 0.05 # p-value threshold (genes)
lfc.thresh <- 0.3 # Log fold change threshold (genes)
n_transcripts <- 542500 # According to the data sheet from Thermofisher
fit <- lmFit(expr.mat, design)
fitC <- contrasts.fit(fit, contrasts)
fitCB <- eBayes(fitC)

