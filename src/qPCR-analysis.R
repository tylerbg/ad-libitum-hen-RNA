# Load and set-up ----------------------------------------------------------------------------------
library(tidyverse)

# Load dataset and modify to proper wide format
qpcr.df <- read.csv("data/raw/qpcr-data-raw.csv", header = TRUE)
rownames(qpcr.df) <- qpcr.df$X
qpcr.df <- qpcr.df[, -1]

# Convert data to long format
qpcr.long <- gather(data.frame(t(qpcr.df)),
                    Gene,
                    ddCt,
                    ANXA2:VIPR2,
                    factor_key = TRUE)

qpcr.long <- qpcr.long %>%
  mutate(Animal = factor(paste(qpcr.long$Age,
                               qpcr.long$Feed),
                         levels = c("10 Week RF",
                                    "16 Week RF",
                                    "10 Week FF",
                                    "16 Week FF")),
         Age = as.factor(Age),
         Feed = as.factor(Feed),
         ID = as.factor(ID),
         ddCt = as.numeric(ddCt),
         ddCt2 = I(2^(-ddCt)))


# Statistical Tests --------------------------------------------------------------------------------
## Two-Factor ANOVAs with type 3 SS -----
# Create a vector of gene names
qpcr.gene.list <- as.vector(rownames(qpcr.df)[4:nrow(qpcr.df)])

# Perform two-way ANOVAs with Type 3 SS and interactions for each gene
qpcr.anovas <- data.frame(sapply(qpcr.gene.list, function(x)
  Anova(aov(ddCt~Feed*Age, data=qpcr.long[qpcr.long$Gene == x,]), type=3)))

## Get all p-values and identify the minimum p-value for each test
qpcr.min.pval <- data.frame(gene=qpcr.gene.list,
                            p.val.Feed=sapply(1:46, function(x)
                              qpcr.anovas[, x]$`Pr(>F)`[2]),
                            p.val.Age=sapply(1:46, function(x)
                              qpcr.anovas[, x]$`Pr(>F)`[3]),
                            p.val.int=sapply(1:46, function(x)
                              qpcr.anovas[, x]$`Pr(>F)`[4]))

qpcr.min.pval$min <- sapply(1:46, function(x)
  min(qpcr.min.pval[x, 2:4]))

# Calculate BH Correction for 46 tests and FDR = 0.05
imq <- seq(1:46) / 46 * 0.05

# Order by the ANOVA tests by their lowest p-value result by ascending and add the BH correction
# column (imq, which is also already ordered by ascending).
qpcr.min.pval.ordered <- (qpcr.min.pval[order(qpcr.min.pval$min), ])
qpcr.min.pval.ordered$imq <- imq

# Add significance marker for interaction term and print
library(rstatix)

# Create final df of SS 3 results and write to reports as csv
qpcr.ss3.res <- qpcr.min.pval.ordered %>%
  add_significance(p.col = 'p.val.int') %>%
  select(!c(min, imq)) %>%
  arrange(gene)

write_csv(qpcr.ss3.res, 'reports/ANOVA-Type3SS-results.csv')


## Two-Factor ANOVAs with type 2 SS ----
# Use type 2 SS for genes without a significant interaction term from the type 3 SS method.

# Run anova with no interaction term and type 2 SS
qpcr.anovas <- data.frame(sapply(qpcr.gene.list, function(x)
  Anova(aov(ddCt~Feed + Age, data=qpcr.long[qpcr.long$Gene == x,]), type=2)))

# Get all p-values then minimum p-value for each test
qpcr.min.pval <- data.frame(gene=qpcr.gene.list,
                            p.val.Feed=sapply(1:46, function(x)
                              qpcr.anovas[, x]$`Pr(>F)`[1]),
                            p.val.Age=sapply(1:46, function(x)
                              qpcr.anovas[, x]$`Pr(>F)`[2])
)

# Convert to a long data frame to treat both factors in the same correction
qpcr.new <- gather(qpcr.min.pval, group, p.val, p.val.Feed:p.val.Age, factor_key = TRUE)
qpcr.new <- qpcr.new[order(qpcr.new$p.val), ]
rownames(qpcr.new) <- seq(1:nrow(qpcr.new))

# Remove HSD17B1, which had a significant interaction term, from type 2 SS analysis
qpcr.new <- qpcr.new %>%
  filter(gene != 'HSD17B1')

# Run BH adjustment
imq <- seq(1:nrow(qpcr.new))/nrow(qpcr.new) * 0.05
qpcr.new$BH.adj.Pval <- qpcr.new$p.val / imq * 0.05

# Add p-value markers separately for the Feed and Age groups
qpcr.new$Feed.sig <- if_else(qpcr.new$BH.adj.Pval <= 0.001 & qpcr.new$group == "p.val.Feed", "+++",
                             if_else(qpcr.new$BH.adj.Pval <= 0.01 & qpcr.new$group == "p.val.Feed", "++",
                                     if_else(qpcr.new$BH.adj.Pval <= 0.05 & qpcr.new$group == "p.val.Feed", "+", "")))
qpcr.new$Age.sig <- if_else(qpcr.new$BH.adj.Pval <= 0.001 & qpcr.new$group == "p.val.Age", "***",
                            if_else(qpcr.new$BH.adj.Pval <= 0.01 & qpcr.new$group == "p.val.Age", "**",
                                    if_else(qpcr.new$BH.adj.Pval <= 0.05 & qpcr.new$group == "p.val.Age", "*", "")))

# Reorder by gene name
qpcr.new <- qpcr.new[order(qpcr.new$gene), ]
rownames(qpcr.new) <- seq(1:nrow(qpcr.new))

# Write results to csv in reports
write_csv(qpcr.new, 'reports/ANOVA-Type2SS-results.csv')


