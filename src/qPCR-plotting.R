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

# Create Plots -------------------------------------------------------------------------------------
library(scales)

# Create an NA df
qpcr.long.NA <- qpcr.long
qpcr.long.NA$ddCt2 <- NA

## * Neuroactive -----
# Set list of neuro genes
neuro.genes <- list("CHRM5", "DRD4", "GABBR2", "GABRA3", "GABRA5", "GRID1",
                    "NELL2", "VIP", "VIPR1", "VIPR2")

# Create and annotate pointrange plot for neuro genes
neuro.plt <- ggplot(qpcr.long[qpcr.long$Gene %in% neuro.genes, ],
       aes(x=Gene,
           y=ddCt2,
           fill=Animal)) +
  theme_bw(base_size = 12) +
  scale_fill_brewer(palette = "Paired") +
  scale_color_brewer(palette = "Paired") +
  stat_summary(aes(color = Animal,
                   shape = Animal),
               fun.data = mean_sdl,
               fun.args=list(mult = 1),
               geom = "pointrange",
               size = 1.25,
               position=position_dodge(0.75)) +
  scale_shape_manual(values = rep("_", 4)) +
  labs(x = NULL,
       y = expression(2^("-"*Delta*Delta*C*"T"))) +
  annotate("text", x=1, y=8, label="+", ) +
  annotate("text", x=2, y=8, label="+++", ) +
  annotate("text", x=4, y=8, label="++", ) +
  annotate("text", x=5, y=8, label="+", ) +
  annotate("text", x=6, y=9, label="****", ) +
  annotate("text", x=7, y=8, label="++", ) +
  annotate("text", x=8, y=9, label="*", ) +
  annotate("text", x=10, y=8, label="+++", ) +
  annotate("text", x=10, y=9, label="***", ) +
  scale_y_continuous(trans="log2",
                     breaks=trans_breaks("log2", function(x) 2^x, n=10),
                     minor_breaks = NULL,
                     labels=trans_format("log2", math_format(2^.x)),
                     limits=c(NA, 10.5)) +
  geom_bar(data = qpcr.long.NA[qpcr.long.NA$Gene %in% neuro.genes, ],
           aes(x = Gene,
               y = ddCt2,
               fill = Animal),
           stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        legend.title = element_blank(),
        text = element_text())

# Write plot to reports
png("reports/figures/qPCR-plots/neuroactive.png")
neuro.plt
dev.off()


## * Steroidogenic -----
# Set list of neuro genes
steroid.genes <- list("CYP11A1", "CYP19A1", "DHCR24", "HMGCR", "HSD17B1",
                      "HSD3B2", "StAR", "StARD4")

steroid.plt <- ggplot(qpcr.long[qpcr.long$Gene %in% steroid.genes, ],
       aes(x=Gene,
           y=ddCt2,
           fill=Animal)) +
  theme_bw(base_size = 12) +
  scale_fill_brewer(palette = "Paired") +
  scale_color_brewer(palette = "Paired") +
  stat_summary(aes(color = Animal,
                   shape = Animal),
               fun.data = mean_sdl,
               fun.args=list(mult = 1),
               geom = "pointrange",
               size = 1.25,
               position=position_dodge(0.75)) +
  scale_shape_manual(values = rep("_", 4)) +
  labs(x = NULL,
       y = expression(2^("-"*Delta*Delta*C*"T"))) +
  annotate("text", x=1, y=8, label="+++", ) +
  annotate("text", x=2, y=8, label="++", ) +
  annotate("text", x=3, y=8, label="++", ) +
  annotate("text", x=4, y=9, label="****", ) +
  annotate("text", x=5, y=8, label="++", ) +
  annotate("text", x=5, y=10.5, label="&", ) +
  annotate("text", x=6, y=8, label="+", ) +
  annotate("text", x=6, y=9, label="*", ) +
  annotate("text", x=7, y=8, label="+++", ) +
  annotate("text", x=7, y=9, label="***", ) +
  annotate("text", x=8, y=8, label="++", ) +
  annotate("text", x=8, y=9, label="**", ) +
  scale_y_continuous(trans="log2",
                     breaks=trans_breaks("log2", function(x) 2^x, n=6),
                     minor_breaks = NULL,
                     labels=trans_format("log2", math_format(2^.x)),
                     limits=c(NA, 10.5)) +
  geom_bar(data = qpcr.long.NA[qpcr.long.NA$Gene %in% steroid.genes, ],
           aes(x = Gene,
               y = ddCt2,
               fill = Animal),
           stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        legend.title = element_blank(),
        text = element_text())

# Write plot to reports
png("reports/figures/qPCR-plots/steroidogenic.png")
steroid.plt
dev.off()


## * ECM -----
# Set list of ECM genes
ecm.genes <- list("ANXA2", "COL1A2", "COL3A1", "COL4A2", "COL6A1", "FHL2",
                  "FN1", "LAMA2", "MMP9", "MMP10", "VCAN")

ecm.plt <- ggplot(qpcr.long[qpcr.long$Gene %in% ecm.genes, ],
       aes(x=Gene,
           y=ddCt2,
           fill=Animal)) +
  theme_bw(base_size = 12) +
  scale_fill_brewer(palette = "Paired") +
  scale_color_brewer(palette = "Paired") +
  stat_summary(aes(color = Animal,
                   shape = Animal),
               fun.data = mean_sdl,
               fun.args=list(mult = 1),
               geom = "pointrange",
               size = 1.25,
               position=position_dodge(0.75)) +
  scale_shape_manual(values = rep("_", 4)) +
  labs(x = NULL,
       y = expression(2^("-"*Delta*Delta*C*"T"))) +
  annotate("text", x=1, y=35, label="+++") +
  annotate("text", x=5, y=35, label="+") +
  annotate("text", x=6, y=35, label="++") +
  annotate("text", x=7, y=35, label="++") +
  annotate("text", x=8, y=35, label="+") +
  annotate("text", x=9, y=35, label="++") +
  annotate("text", x=11, y=35, label="+++") +
  annotate("text", x=1, y=39.5, label=sprintf('**')) +
  annotate("text", x=6, y=39.5, label=sprintf('**')) +
  annotate("text", x=8, y=39.5, label=sprintf('**')) +
  annotate("text", x=11, y=39.5, label=sprintf('***')) +
  scale_y_continuous(trans="log2",
                     breaks=trans_breaks("log2", function(x) 2^x, n=8),
                     minor_breaks = NULL,
                     labels=trans_format("log2", math_format(2^.x)),
                     limits=c(NA, 40)) +
  geom_bar(data = qpcr.long.NA[qpcr.long.NA$Gene %in% ecm.genes, ],
           aes(x = Gene,
               y = ddCt2,
               fill = Animal),
           stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        legend.title = element_blank(),
        text = element_text())

# Write plot to reports
png("reports/figures/qPCR-plots/ecm.png")
ecm.plt
dev.off()


## * Cell Adhesion -----
# Set list of cell adhesion genes
adhesion.genes <- list("CCDC80", "CD36", "CDH13", "CLDN18", "MGLL", "MYH11",
                       "NEGR1", "POSTN", "PTEN")

adhesion.plt <- ggplot(qpcr.long[qpcr.long$Gene %in% adhesion.genes, ],
       aes(x=Gene,
           y=ddCt2,
           fill=Animal)) +
  theme_bw(base_size = 12) +
  scale_fill_brewer(palette = "Paired") +
  scale_color_brewer(palette = "Paired") +
  stat_summary(aes(color = Animal,
                   shape = Animal),
               fun.data = mean_sdl,
               fun.args=list(mult = 1),
               geom = "pointrange",
               size = 1.25,
               position=position_dodge(0.75)) +
  scale_shape_manual(values = rep("_", 4)) +
  labs(x = NULL,
       y = expression(2^("-"*Delta*Delta*C*"T"))) +
  annotate("text", x=1, y=30, label="++++") +
  annotate("text", x=3, y=30, label="+") +
  annotate("text", x=5, y=30, label="+++") +
  annotate("text", x=6, y=30, label="++") +
  annotate("text", x=8, y=30, label="+") +
  annotate("text", x=9, y=30, label="+") +
  annotate("text", x=1, y=35, label=sprintf('****')) +
  annotate("text", x=3, y=35, label=sprintf('*')) +
  annotate("text", x=4, y=35, label=sprintf('*')) +
  annotate("text", x=5, y=35, label=sprintf('**')) +
  annotate("text", x=6, y=35, label=sprintf('***')) +
  annotate("text", x=7, y=35, label=sprintf('***')) +
  annotate("text", x=9, y=35, label=sprintf('****')) +
  scale_y_continuous(trans="log2",
                     breaks=trans_breaks("log2", function(x) 2^x, n=10),
                     minor_breaks = NULL,
                     labels=trans_format("log2", math_format(2^.x)),
                     limits=c(NA, 40)) +
  geom_bar(data = qpcr.long.NA[qpcr.long.NA$Gene %in% adhesion.genes, ],
           aes(x = Gene,
               y = ddCt2,
               fill = Animal),
           stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        legend.title = element_blank(),
        text = element_text())

# Write plot to reports
png("reports/figures/qPCR-plots/cell-adhesion.png")
adhesion.plt
dev.off()


## * Transcription and defense -----
# Set list of transcription and defense genes
trandef.genes <- list("CATHL3", "DLX5", "GAL7", "LECT2", "RBP4",
                      "SATB2", "SUPT3H", "TEAD1")

trandef.plt <- ggplot(qpcr.long[qpcr.long$Gene %in% trandef.genes, ],
       aes(x=Gene,
           y=ddCt2,
           fill=Animal)) +
  theme_bw(base_size = 12) +
  scale_fill_brewer(palette = "Paired") +
  scale_color_brewer(palette = "Paired") +
  stat_summary(aes(color = Animal,
                   shape = Animal),
               fun.data = mean_sdl,
               fun.args=list(mult = 1),
               geom = "pointrange",
               size = 1.25,
               position=position_dodge(0.75)) +
  scale_shape_manual(values = rep("_", 4)) +
  labs(x = NULL,
       y = expression(2^("-"*Delta*Delta*C*"T"))) +
  annotate("text", x=2, y=15, label=sprintf('*')) +
  annotate("text", x=6, y=15, label=sprintf('**')) +
  annotate("text", x=7, y=15, label=sprintf('**')) +
  scale_y_continuous(trans="log2",
                     breaks=trans_breaks("log2", function(x) 2^x, n=10),
                     minor_breaks = NULL,
                     labels=trans_format("log2", math_format(2^.x)),
                     limits=c(NA, 16)) +
  geom_bar(data = qpcr.long.NA[qpcr.long.NA$Gene %in% trandef.genes, ],
           aes(x = Gene,
               y = ddCt2,
               fill = Animal),
           stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        legend.title = element_blank(),
        text = element_text())

# Write plot to reports
png("reports/figures/qPCR-plots/transcription-defence.png")
trandef.plt
dev.off()
