library(tidyverse)
library(oligo)
library(RColorBrewer)

#-----------------#
#### Read data ####
#-----------------#

# Get the raw data
setwd("/home/svillicaña/INMEGEN_LB")
celfiles <- read.table('data/sindrome_metabolico_group_file.txt',
                       header = T, as.is = T)
raw.data <- read.celfiles(celfiles$file)

# Labels
contr.lab <- unique(celfiles$group)
sampleNames(raw.data) <- celfiles$sample
palette <- factor(celfiles$group, 
                  labels = brewer.pal(n = n_distinct(celfiles$group), name = "Dark2")) %>%
  as.character()

# Directories
rds.dir <- "/home/svillicaña/INMEGEN_LB/resultados/RDS/expression/"
plot.dir <- "/home/svillicaña/INMEGEN_LB/resultados/Plots/expression/"
html.dir <- "/home/svillicaña/INMEGEN_LB/resultados/HTML/expression/"
csv.dir <- "/home/svillicaña/INMEGEN_LB/resultados/CSV/expression/"

#----------------#
#### QC plots ####
#----------------#

## Data before the normalization
# target="core" is used to visualize data at transcript level
# Boxplot
pdf(paste0(plot.dir, "raw_data_boxplot.pdf"), width = 16, height = 9)
oligo::boxplot(raw.data, target = "core", col = palette,
               main = "Raw data distribution", cex = 0.8, las = 2)
dev.off()

# Density plot
pdf(paste0(plot.dir, "raw_data_density.pdf"), width = 16, height = 9)
oligo::hist(raw.data, target = "core", col = palette,
            main = "Raw data density curves", cex = 0.8, lwd = 2)
dev.off()

# PCA
exp.raw <- log2(Biobase::exprs(raw.data))
PCA.raw <- prcomp(t(exp.raw), scale. = FALSE)

percentVar <- round(100*PCA.raw$sdev^2/sum(PCA.raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA.raw$x[,1], PC2 = PCA.raw$x[,2],
                     Group = celfiles$group)

pdf(paste0(plot.dir, "raw_data_pca.pdf"), width = 16, height = 9)
ggplot(dataGG, aes(PC1, PC2)) +
       geom_point(aes(colour = Group)) +
  ggtitle("PCA plot of the log-transformed raw data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio) +
  scale_color_brewer(palette="Dark2")
dev.off()

