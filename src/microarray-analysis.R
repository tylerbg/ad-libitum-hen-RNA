# Load and set-up ----------------------------------------------------------------------------------
# Install Bioconductor and required libraries
if(!requireNamespace("BiocManager", quietly = TRUE))install.packages("BiocManager")
BiocManager::install()

BiocManager::install(c("affy", "affyPLM", "gcrma", "limma", "fifer"))

library(affy)
library(affyPLM)

# Set paths for then read in CEL files as `data`
pwd <- getwd()
celpath <- paste(pwd, "/data/raw/CEL", sep="")

data <- ReadAffy(celfile.path=celpath)

# Retrieve intensities from each CEL as `int`, annotations as `ph`, and features as `feat`
int <- intensity(data)
ph <- data@phenoData
feat <- data@featureData


# Quality control ----------------------------------------------------------------------------------
# Generate individual microarray images for all 6 samples
for (i in 1:6) {
  name <- paste0("reports/figures/microarray-images/Sample", i, ".jpg")
  jpeg(name)
  image(data[, i], main = ph@data$sample[i])
  dev.off()
}

# Create chip psuedo-images by fitting robust Probe Level linear models
Pset <- fitPLM(data)

for (i in 1:6) {
  name <- paste0("reports/figures/microarray-images/PLM", i, ".jpg")
  jpeg(name)
  image(Pset, which = i, main = ph@data$sample[i])
  dev.off()
}

# Draw a density plot of the microarray data
# Blue represents the ad-libitem animals and yellow the restricted fed hens
color = c(rep('steelblue3', 3),
          rep('goldenrod3', 3))
hist(data[, 1:6],
     lwd = 2,
     which = 'pm',
     col = color,
     ylab = 'Density',
     xlab = 'Log2 intensities',
     main = 'Density plot of the raw data')

# Generate MA plots
for (i in 1:6) {
  name = paste("reports/figures/MA-plots/MAplot", i, ".jpg", sep = "")
  jpeg(name)
  MAplot(data, which = i)
  dev.off()
}


# Normalization ------------------------------------------------------------------------------------
# Normalize using RMA method
data.rma <- rma(data)
data.matrix <- exprs(data.rma)

# Compare raw and background-corrected data
bgcorr.rma <- pm(bg.correct(data, method = "rma"))
pmexp <- pm(data)

# Pre-make vectors for loop
sampleNames <- vector()
logs <- vector()
corrlogs <- vector()

# Calculate and print images for log-background-correlations
for (i in 1:6) {
  sampleNames <- c(sampleNames, rep(ph@data[i, 1], dim(pmexp)[1]))
  logs <- c(logs, log2(pmexp[, i]))
  corrlogs <- c(corrlogs, log2(bgcorr.rma[, i]))
}

corrData <- data.frame(logInt = logs,
                       bgcorr_logInt = corrlogs,
                       sampleName = sampleNames)


# Takes some time to produce image
library(tidyverse)

jpeg('reports/figures/bgcorr.jpeg')
ggplot(corrData, aes(logInt, bgcorr_logInt)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, colour = 'tomato2') +
  facet_wrap(. ~ sampleName)
dev.off()


# EDA - PCA & clustering ---------------------------------------------------------------------------
# Get pricinpal components then plot first 3 PCs against one another
data.PC <- prcomp(t(data.matrix),
                  scale. = TRUE,
                  center = TRUE)

jpeg('reports/figures/PCA.jpeg')
par(mfrow = c(2, 2))
plot(data.PC$x[, 1:2],
     col = color,
     pch = 18,
     cex = 2)
plot(data.PC$x[, 2:3],
     col = color,
     pch = 18,
     cex = 2)
plot(data.PC$x[, c(1, 3)],
     col = color,
     pch = 18,
     cex = 2)
par(mfrow = c(1, 1))
dev.off()


# Install and load facto packages for clustering analyses
if(!require(FactoMineR)){install.packages("FactoMineR")+library(FactoMineR)}
if(!require(factoextra)){install.packages("factoextra")+library(factoextra)}

# Set up matrix
pca.matrix <- data.matrix
colnames(pca.matrix) <- c("RF1", "RF2", "RF3", "FF1", "FF2", "FF3")
pca.matrix <- t(pca.matrix)

# Run PCA and hierarchical clustering on principal components
res.pca <- PCA(pca.matrix, ncp = 10, graph = FALSE, scale.unit = TRUE)
res.hcpc <- HCPC(res.pca, min = 2, max = 2, graph = FALSE)

# Plot dendrogram and factor maps

jpeg('reports/figures/dendrogram.jpeg')
fviz_dend(res.hcpc, 
          cex = 0.7,                     # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          labels_track_height = 0.8      # Augment the room for labels
)
dev.off()

jpeg('reports/figures/factor-map.jpeg')
fviz_cluster(res.hcpc,
             repel = TRUE,            # Avoid label overlapping
             show.clust.cent = TRUE, # Show cluster centers
             palette = "jco",         # Color palette see ?ggpubr::ggpar
             ggtheme = theme_minimal(),
             main = "Factor map"
)
dev.off()

# Plot 3-dimensional hierarchical clustering on factor map
jpeg('reports/figures/3d-HC-factor-map.jpeg')
plot(res.hcpc, choice = "3D.map")
dev.off()


# Differential Expression --------------------------------------------------------------------------
library(limma)

# Set diet labels in annotation data then use them to set factored groups
ph@data[, 2] <- c(rep("RF", 3), rep("FF", 3))
colnames(ph@data)[2] <- "source"

groups <- ph@data$source
f <- factor(groups, levels = c("RF", "FF"))

# Design the model matrix then fit using a linear model for a series of arrays
design <- model.matrix(~ 0 + f)
colnames(design) <- c("RF", "FF")

data.fit <- lmFit(data.matrix, design)

# Build a contrast matrix of the differences between restricted and ad libitum feeding (full feed)
# then compute contrasts from the linear model
contrast.matrix <- makeContrasts(RF-FF, levels = design)
data.fit.con <- contrasts.fit(data.fit, contrast.matrix)

# Assess differentially expressed genes using empirical Bayes stats
data.fit.eb <- eBayes(data.fit.con)

# Create table for genes with at least 2-fold 
n <- nrow(data.matrix)
tab.2fold <- topTable(data.fit.eb,
                      lfc = 1,
                      number = n)

# Test only those with 2-fold change or greater using BH correction with FDR of 0.05
genes.2fold <- topTable(data.fit.eb,
                        coef = 1,
                        number = 1e6)
genes.2fold <- genes.2fold[abs(genes.2fold$logFC) >= 1, ]

i <- seq(1:nrow(genes.2fold))
m <- nrow(genes.2fold)
q <- 0.05
genes.2fold$imq <- i / m * q
topgenes.2fold <- genes.2fold[genes.2fold[, "P.Value"] < genes.2fold[, "imq"], ]



# Volcanoplot identifying boundaries and genes above FDR threshold and more than 2-fold change
jpeg('reports/figures/volcano-plot.jpeg')
volcanoplot(data.fit.eb, coef=1, highlight=0,
            xlim=c(-3.25,3.25))
grid()
points(data.fit.eb$coefficients, -log10(data.fit.eb$p.value),
       pch=16, cex=0.35)
points(data.fit.eb$coefficients[rownames(data.fit.eb$coefficients)
                                %in% rownames(topgenes.2fold)],
       -log10(data.fit.eb$p.value[rownames(data.fit.eb$coefficients)
                                  %in% rownames(topgenes.2fold)]),
       col="steelblue3",
       pch=16,
       cex=0.35)
abline(v=1, col="tomato", lty=2)
abline(v=-1, col="tomato", lty=2)
abline(h=-log10(max(topgenes.2fold$P.Value)), col="tomato", lty=2)
text(0, 4.5, labels="2-fold change", col="tomato")
text(-3, 2, labels=paste("FDR=", q, sep=""), col="tomato")
text(-2.9, 5.5, labels=paste(sum(topgenes.2fold$logFC <= -1), " genes\n",
                             "up in FF", sep=""), font=2, col="steelblue3")
text(2.9, 5.5, labels=paste(sum(topgenes.2fold$logFC >= 1), " genes\n",
                            "up in RF", sep=""), font=2, col="steelblue3")
dev.off()




