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

#----------------------------------#
#### Normalization by quantiles ####
#----------------------------------#

# Create list of index by group
groups_list <- lapply(as.list(contr.lab),
                      function(x) which(celfiles$group == x))
names(groups_list) <- contr.lab
norm.m <- "quantile"
data.norm <- raw.data

# Apply normalization by group
for(x in groups_list){
  q.norm <- normalize(raw.data[, x], method = norm.m)
  pm(data.norm)[, x] <- pm(q.norm)
}

## Plots after quantiles normalization
# Boxplot
pdf(paste0(plot.dir, "qnorm_boxplot.pdf"), width = 16, height = 9)
oligo::boxplot(data.norm, target = "core", col = palette,
               main = "Quantile-normalized by group", cex = 0.8, las = 2)
dev.off()

# Density
pdf(paste0(plot.dir, "qnorm_density.pdf"), width = 16, height = 9)
oligo::hist(data.norm, target = "core", col = palette,
            main = "Quantile-normalized by group", cex = 0.8, lwd = 2)
dev.off()

# PCA
exp.norm <- log2(Biobase::exprs(data.norm))
PCA.norm <- prcomp(t(exp.norm), scale. = FALSE)

percentVar <- round(100*PCA.norm$sdev^2/sum(PCA.norm$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA.norm$x[,1], PC2 = PCA.norm$x[,2],
                     Group = celfiles$group)

pdf(paste0(plot.dir, "qnorm_pca.pdf"), width = 16, height = 9)
ggplot(dataGG, aes(PC1, PC2)) +
       geom_point(aes(colour = Group)) +
  ggtitle("Quantile-normalized by group") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio) +
  scale_color_brewer(palette="Dark2")
dev.off()

#----------------------------#
#### Normalization by RMA ####
#----------------------------#

data.RMA <- rma(data.norm, target = "core")
expr.mat <- exprs(data.RMA)

## Plots after RMA normalization
# Boxplot
pdf(paste0(plot.dir, "RMA_boxplot.pdf"), width = 16, height = 9)
oligo::boxplot(data.RMA, col = palette, main = "RMA-normalized data (all)",
               las = 2, ylab = "Expression levels")
dev.off()

# Density
pdf(paste0(plot.dir, "RMA_density.pdf"), width = 16, height = 9)
oligo::hist(data.RMA, col = palette, main = "RMA-normalized data (all)",
            cex = 0.8, lwd = 2)
dev.off()

# PCA
PCA.RMA <- prcomp(t(log2(expr.mat)), scale. = FALSE)

percentVar <- round(100*PCA.RMA$sdev^2/sum(PCA.RMA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA.RMA$x[,1], PC2 = PCA.RMA$x[,2],
                     Group = celfiles$group)

pdf(paste0(plot.dir, "RMA_pca.pdf"), width = 16, height = 9)
ggplot(dataGG, aes(PC1, PC2)) +
       geom_point(aes(colour = Group)) +
  ggtitle("RMA-normalized data (all)") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio) +
  scale_color_brewer(palette="Dark2")
dev.off()

#-----------------#
#### Save data ####
#-----------------#
  
saveRDS(data.norm, file = paste0(rds.dir, "data.norm.rds"))
saveRDS(data.RMA, file = paste0(rds.dir, "data.RMA.rds"))

