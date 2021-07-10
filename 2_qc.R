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
expr.RMA <- log2(exprs(data.RMA))
PCA.RMA <- prcomp(t(expr.RMA), scale. = FALSE)

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

#--------------------------#
#### Sample correlation ####
#--------------------------#

# Correlation matrix
expr.RMA.cor <- round(cor(expr.RMA),2)
expr.RMA.cor[lower.tri(expr.RMA.cor)]<- NA

# Correlation heatmap
expr.RMA.cor.long <- as.data.frame(expr.RMA.cor) %>%
  rownames_to_column(var = "Sample_1") %>%
  pivot_longer(!Sample_1, names_to = "Sample_2", values_to = "Correlation", values_drop_na = TRUE) %>%
  mutate(Sample_1 = factor(Sample_1, levels = celfiles$sample),
         Sample_2 = factor(Sample_2, levels = celfiles$sample))

pdf(paste0(plot.dir, "RMA_cor.pdf"), width = 20, height = 20)
ggplot(data = expr.RMA.cor.long, aes(Sample_1, Sample_2, fill = Correlation)) +
  geom_tile(color = "white") +
  # scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
  #                     midpoint = 0, limit = c(-1,1), space = "Lab", 
  #                     name="Pearson\nCorrelation") +
  theme_minimal()+ 
  geom_text(aes(Sample_1, Sample_2, label = Correlation), color = "black", size = 3) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, 
                               size = 12, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = "bottom",
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5)) + 
  coord_fixed()
dev.off()

#-----------------#
#### Save data ####
#-----------------#
  
saveRDS(data.norm, file = paste0(rds.dir, "data.norm.rds"))
saveRDS(data.RMA, file = paste0(rds.dir, "data.RMA.rds"))

