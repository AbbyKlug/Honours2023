### libraries ###
library(ggplot2)
library(DESeq2)
library(tximport)
library(tidyverse)
library(ggsignif)
library(ggrepel)
library(clusterProfiler)
library(VennDiagram)
library(enrichplot)
library(stringr)
library(forcats)
library(org.Hs.eg.db)
library(ggpubr)

install.packages("VennDiagram")
install.packages("ggsignif")
install.packages("ggrepel")
install.packages("egg")
install.packages("ggpubr")
library(egg)
## Yi's tx2gene
tx2gene = readRDS("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/DIFF-EXP/tx2gene.rds")
setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/FASTQ.files/")

tx2gene = readRDS("D://ABBY.windows_surface_pro/DIFF-EXP/tx2gene.rds")
setwd("D://ABBY.windows_surface_pro/FASTQ.files/")

## MALE
# files = c('./AD/20-1T-AD.fastq.gz_counts/quant.sf',
#           './AD/21-1A-AD.fastq.gz_counts/quant.sf',
#           './AD/22-2T-AD.fastq.gz_counts/quant.sf',
#           './AD/23-2A-AD.fastq.gz_counts/quant.sf',
#           './AD/24-3T-AD.fastq.gz_counts/quant.sf',
#           './AD/25-5T-AD.fastq.gz_counts/quant.sf',
#           './AD/26-3A-AD.fastq.gz_counts/quant.sf',
#           './AD/28-8T-AD.fastq.gz_counts/quant.sf',
#           './AD/29-6T-AD.fastq.gz_counts/quant.sf',
#           './AD/30-9T-AD.fastq.gz_counts/quant.sf',
#           './AD/31-7T-AD.fastq.gz_counts/quant.sf',
#           './Old/10-8A-Old.fastq.gz_counts/quant.sf',
#           './Old/11-10T-Old.fastq.gz_counts/quant.sf',
#           './Old/12-6A-Old.fastq.gz_counts/quant.sf',
#           './Old/13-11T-Old.fastq.gz_counts/quant.sf',
#           './Old/14-7A-Old.fastq.gz_counts/quant.sf',
#           './Old/15-13T-Old.fastq.gz_counts/quant.sf',
#           './Old/16-14T-Old.fastq.gz_counts/quant.sf',
#           './Old/17-9A-Old.fastq.gz_counts/quant.sf',
#           './Old/18-10A-Old.fastq.gz_counts/quant.sf',
#           './Young/2-12A-Young.fastq.gz_counts/quant.sf',
#           './Young/3-17T-Young.fastq.gz_counts/quant.sf',
#           './Young/4-13A-Young.fastq.gz_counts/quant.sf',
#           './Young/5-18T-Young.fastq.gz_counts/quant.sf',
#           './Young/6-14A-Young.fastq.gz_counts/quant.sf',
#           './Young/7-19T-Young.fastq.gz_counts/quant.sf',
#           './Young/8-15A-Young.fastq.gz_counts/quant.sf')

###### AD vs Old ######

files = c('./AD/20-1T-AD.fastq.gz_counts/quant.sf',
          './AD/21-1A-AD.fastq.gz_counts/quant.sf',
          './AD/22-2T-AD.fastq.gz_counts/quant.sf',
          './AD/23-2A-AD.fastq.gz_counts/quant.sf',
          './AD/24-3T-AD.fastq.gz_counts/quant.sf',
          './AD/25-5T-AD.fastq.gz_counts/quant.sf',
          './AD/26-3A-AD.fastq.gz_counts/quant.sf',
          './AD/28-8T-AD.fastq.gz_counts/quant.sf',
          './AD/29-6T-AD.fastq.gz_counts/quant.sf',
          './AD/30-9T-AD.fastq.gz_counts/quant.sf',
          './AD/31-7T-AD.fastq.gz_counts/quant.sf',
          './Old/10-8A-Old.fastq.gz_counts/quant.sf',
          './Old/11-10T-Old.fastq.gz_counts/quant.sf',
          './Old/12-6A-Old.fastq.gz_counts/quant.sf',
          './Old/13-11T-Old.fastq.gz_counts/quant.sf',
          './Old/14-7A-Old.fastq.gz_counts/quant.sf',
          './Old/15-13T-Old.fastq.gz_counts/quant.sf',
          './Old/16-14T-Old.fastq.gz_counts/quant.sf',
          './Old/17-9A-Old.fastq.gz_counts/quant.sf',
          './Old/18-10A-Old.fastq.gz_counts/quant.sf')

txi.salmon <- tximport(files = files, 
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)

sampleTable <- data.frame(condition = factor(c("AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","Old","Old","Old","Old","Old","Old","Old","Old","Old")))
names <-  (c("AD1", "AD2","AD3","AD4","AD5","AD6","AD7","AD8","AD9","AD10","AD11",
             "Old1","Old2","Old3","Old4","Old5","Old6","Old7","Old8","Old9"))

rownames(sampleTable) = names
colnames(txi.salmon$counts) <- names

dds_ADvsOld <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~ condition)
dds_ADvsOld$condition <- relevel(dds_ADvsOld$condition, ref = "Old")
dds_ADvsOld <- DESeq(dds_ADvsOld)
resultsNames(dds_ADvsOld)

## PCA ##

vsd <- vst(dds_ADvsOld, blind = FALSE) # or varianceStabilizingTransformation instead of vst

# Create the PCA plot with custom colors

PCA_ADOld <- plotPCA(vsd, intgroup="condition") +
  geom_point(shape = 21, aes(fill = condition, colour=condition), size = 4) +  # Use points (dots) for labels
  geom_text(aes(label = name, fontface="bold"), colour = "black", nudge_y = 1.35, nudge_x = 0.35, size = 4) +
  scale_fill_manual(values = c("AD" = "orangered", "Old" = "darkslateblue")) +
  scale_colour_manual(values = c("darkslateblue", "orangered")) +
  labs(fill="", colour="") +
  theme(legend.position = "top")
PCA_ADOld  

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}

PCA_ADOld <- edit_plots(PCA_ADOld)
PCA_ADOld 

?plotPCA
setwd("D://ABBY.windows_surface_pro/diffexp_results/")
saveRDS(PCA_ADOld, "./PCA_ADOld_all_11sept.rds")

#geom_label(aes(label = name, colour = condition))


setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp_update/")
saveRDS(PCA_ADOld, "./PCA_ADOld_all_b&w.rds")

###### Old vs Young ######

files = c('./Old/10-8A-Old.fastq.gz_counts/quant.sf',
          './Old/11-10T-Old.fastq.gz_counts/quant.sf',
          './Old/12-6A-Old.fastq.gz_counts/quant.sf',
          './Old/13-11T-Old.fastq.gz_counts/quant.sf',
          './Old/14-7A-Old.fastq.gz_counts/quant.sf',
          './Old/15-13T-Old.fastq.gz_counts/quant.sf',
          './Old/16-14T-Old.fastq.gz_counts/quant.sf',
          './Old/17-9A-Old.fastq.gz_counts/quant.sf',
          './Old/18-10A-Old.fastq.gz_counts/quant.sf',
          './Young/2-12A-Young.fastq.gz_counts/quant.sf',
          './Young/3-17T-Young.fastq.gz_counts/quant.sf',
          './Young/4-13A-Young.fastq.gz_counts/quant.sf',
          './Young/5-18T-Young.fastq.gz_counts/quant.sf',
          './Young/6-14A-Young.fastq.gz_counts/quant.sf',
          './Young/7-19T-Young.fastq.gz_counts/quant.sf',
          './Young/8-15A-Young.fastq.gz_counts/quant.sf')

txi.salmon <- tximport(files = files, 
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)

sampleTable <- data.frame(condition = factor(c("Old","Old","Old","Old","Old","Old","Old","Old","Old", "Young","Young","Young","Young","Young","Young","Young")))
names <-  (c("Old1","Old2","Old3","Old4","Old5","Old6","Old7","Old8","Old9",
             "Young1","Young2","Young3","Young4","Young5","Young6","Young7"))

rownames(sampleTable) = names
colnames(txi.salmon$counts) <- names

dds_OldvsYoung <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~ condition)
dds_OldvsYoung$condition <- relevel(dds_OldvsYoung$condition, ref = "Young")
dds_OldvsYoung <- DESeq(dds_OldvsYoung)
resultsNames(dds_OldvsYoung)

## PCA ##

vsd <- vst(dds_OldvsYoung, blind = FALSE) # or varianceStabilizingTransformation instead of vst

PCA_OldYoung <- plotPCA(vsd, intgroup="condition") +
  geom_point(shape = 21, aes(fill = condition, colour=condition), size = 4) +  # Use points (dots) for labels
  geom_text(aes(label = name, fontface="bold"), colour = "black", nudge_y = 1.35, nudge_x = 0.35, size = 4) +
  scale_fill_manual(values = c("Old" = "darkslateblue", "Young" = "gold")) +
  scale_colour_manual(values = c("gold", "darkslateblue")) +
  labs(fill="", colour="") +
  theme(legend.position = "top")
PCA_OldYoung 

# PCA_OldYoung <- plotPCA(vsd, intgroup="condition") +
#   geom_point(aes(color = condition), size = 1) +  # Use points (dots) for labels
#   geom_text(aes(label = name, fontface="bold"), colour = "black", nudge_y = 1.85, nudge_x = 0.35, size = 4) +
#   scale_colour_manual(values = c("Old" = "plum", "Young"="orange")) 

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}

PCA_OldYoung <- edit_plots(PCA_OldYoung)
PCA_OldYoung

setwd("D://ABBY.windows_surface_pro/diffexp_results/")
saveRDS(PCA_OldYoung, "./PCA_OldYoung_all_11sept.rds")

setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp_update/")
saveRDS(PCA_OldYoung, "./PCA_OldYoung_all_b&w.rds")

#plotPCA(vsd, intgroup="condition") +
#  geom_label(aes(label = name))

###### AD vs Young ######

files = c('./AD/20-1T-AD.fastq.gz_counts/quant.sf',
          './AD/21-1A-AD.fastq.gz_counts/quant.sf',
          './AD/22-2T-AD.fastq.gz_counts/quant.sf',
          './AD/23-2A-AD.fastq.gz_counts/quant.sf',
          './AD/24-3T-AD.fastq.gz_counts/quant.sf',
          './AD/25-5T-AD.fastq.gz_counts/quant.sf',
          './AD/26-3A-AD.fastq.gz_counts/quant.sf',
          './AD/28-8T-AD.fastq.gz_counts/quant.sf',
          './AD/29-6T-AD.fastq.gz_counts/quant.sf',
          './AD/30-9T-AD.fastq.gz_counts/quant.sf',
          './AD/31-7T-AD.fastq.gz_counts/quant.sf',
          './Young/2-12A-Young.fastq.gz_counts/quant.sf',
          './Young/3-17T-Young.fastq.gz_counts/quant.sf',
          './Young/4-13A-Young.fastq.gz_counts/quant.sf',
          './Young/5-18T-Young.fastq.gz_counts/quant.sf',
          './Young/6-14A-Young.fastq.gz_counts/quant.sf',
          './Young/7-19T-Young.fastq.gz_counts/quant.sf',
          './Young/8-15A-Young.fastq.gz_counts/quant.sf')

txi.salmon <- tximport(files = files, 
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)

sampleTable <- data.frame(condition = factor(c("AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","Young","Young","Young","Young","Young","Young","Young")))
names <-  (c("AD1", "AD2","AD3","AD4","AD5","AD6","AD7","AD8","AD9","AD10","AD11",
             "Young1","Young2","Young3","Young4","Young5","Young6","Young7"))


rownames(sampleTable) = names
colnames(txi.salmon$counts) <- names

dds_ADvsYoung <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~ condition)
dds_ADvsYoung$condition <- relevel(dds_ADvsYoung$condition, ref = "Young")
dds_ADvsYoung <- DESeq(dds_ADvsYoung)
resultsNames(dds_ADvsYoung)

## PCA ##

vsd <- vst(dds_ADvsYoung, blind = FALSE) # or varianceStabilizingTransformation instead of vst

PCA_ADYoung <- plotPCA(vsd, intgroup="condition") +
  geom_point(shape = 21, aes(fill = condition, colour=condition), size = 4) +  # Use points (dots) for labels
  geom_text(aes(label = name, fontface="bold"), colour = "black", nudge_y = 1.35, nudge_x = 0.35, size = 4) +
  scale_fill_manual(values = c("AD" = "orangered", "Young" = "gold")) +
  scale_colour_manual(values = c("gold", "orangered")) +
  labs(fill="", colour="") +
  theme(legend.position = "top")
PCA_ADYoung

# PCA_ADYoung <- plotPCA(vsd, intgroup="condition") +
#   geom_point(aes(color = condition), size = 1) +  # Use points (dots) for labels
#   geom_text(aes(label = name, fontface="bold"), colour = "black", nudge_y = 1.85, nudge_x = 0.35, size = 4) +
#   scale_colour_manual(values = c("AD" = "slategray3", "Young"="orange")) 

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}

PCA_ADYoung <- edit_plots(PCA_ADYoung)
PCA_ADYoung

setwd("D://ABBY.windows_surface_pro/diffexp_results/")
saveRDS(PCA_ADYoung, "./PCA_ADYoung_all_11sept.rds")

setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp_update/")
saveRDS(PCA_ADYoung, "./PCA_ADYoung_all_b&w.rds")

#plotPCA(vsd, intgroup="condition") +
#  geom_label(aes(label = name))

### join 3 PCA plots (all samples) ###
library(ggpubr)

setwd("D://ABBY.windows_surface_pro/differentialexp_update/")

plot1_all <- readRDS("./PCA_ADOld_all_9sept.rds")
plot2_all <- readRDS("./PCA_OldYoung_all_9sept.rds")
plot3_all <- readRDS("./PCA_ADYoung_all_9sept.rds")
plot4 <- readRDS("./PCA_ADOld_3each_9sept.rds")
plot5 <- readRDS("./PCA_OldYoung_3each_9sept.rds")
plot6 <- readRDS("./PCA_ADYoung_3each_9sept.rds")

multi_panel_figure_all <- ggarrange(plot1_all, plot2_all, plot3_all, ncol = 3,
                                    common.legend = FALSE,   # Add a single legend for all plots
                                    legend = "top", # Adjust the legend position as needed
                                    labels = "AUTO",
                                    label.x = 0.0,
                                    align = "hv",
                                    vjust = 1.8
)

multi_panel_figure_3s <- ggarrange(plot4, plot5, plot6, ncol = 3,
                                    common.legend = FALSE,   # Add a single legend for all plots
                                    legend = "top", # Adjust the legend position as needed
                                    labels = "AUTO",
                                    label.x = 0.0,
                                    align = "hv",
                                    vjust = 1.8
)

# Display the multi-panel figure
multi_panel_figure_all
multi_panel_figure_3s
saveRDS(multi_panel_figure_all, "./multipanel_PCA_all_9sept.rds")
saveRDS(multi_panel_figure_3s, "./multipanel_PCA_3each_9sept.rds")
####### DGE SELECTING SAMPLES PER CONDITION #######
# files = c('./AD/22-2T-AD.fastq.gz_counts/quant.sf', #3
#           './AD/25-5T-AD.fastq.gz_counts/quant.sf', #6
#           './AD/26-3A-AD.fastq.gz_counts/quant.sf', #7
#           './AD/27-5A-AD.fastq.gz_counts/quant.sf', #8
#           './AD/28-8T-AD.fastq.gz_counts/quant.sf', #9
#           './Old/10-8A-Old.fastq.gz_counts/quant.sf', #1
#           './Old/11-10T-Old.fastq.gz_counts/quant.sf',#2
#           './Old/12-6A-Old.fastq.gz_counts/quant.sf', #3
#           './Old/13-11T-Old.fastq.gz_counts/quant.sf', #4
#           './Old/14-7A-Old.fastq.gz_counts/quant.sf', #5
#           './Old/15-13T-Old.fastq.gz_counts/quant.sf', #6
#           './Old/16-14T-Old.fastq.gz_counts/quant.sf', #7
#           './Young/2-12A-Young.fastq.gz_counts/quant.sf', #1
#           './Young/3-17T-Young.fastq.gz_counts/quant.sf', #2
#           './Young/4-13A-Young.fastq.gz_counts/quant.sf', #3
#           './Young/6-14A-Young.fastq.gz_counts/quant.sf') #5

## Yi's tx2gene
# tx2gene = readRDS("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/DIFF-EXP/tx2gene.rds")
# setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/FASTQ.files/")
# ###### AD vs Old ######
# 
# files = c('./AD/22-2T-AD.fastq.gz_counts/quant.sf', #3
#           './AD/25-5T-AD.fastq.gz_counts/quant.sf', #6
#           './AD/26-3A-AD.fastq.gz_counts/quant.sf', #7
#           './AD/27-5A-AD.fastq.gz_counts/quant.sf', #8
#           './AD/28-8T-AD.fastq.gz_counts/quant.sf', #9
#           './Old/10-8A-Old.fastq.gz_counts/quant.sf', #1
#           './Old/11-10T-Old.fastq.gz_counts/quant.sf',#2
#           './Old/12-6A-Old.fastq.gz_counts/quant.sf', #3
#           './Old/13-11T-Old.fastq.gz_counts/quant.sf', #4
#           './Old/14-7A-Old.fastq.gz_counts/quant.sf', #5
#           './Old/15-13T-Old.fastq.gz_counts/quant.sf', #6
#           './Old/16-14T-Old.fastq.gz_counts/quant.sf') #7
# 
# txi.salmon <- tximport(files = files, 
#                        type = "salmon",
#                        txOut = FALSE,
#                        tx2gene = tx2gene,
#                        ignoreAfterBar = T)
# 
# sampleTable <- data.frame(condition = factor(c("AD", "AD","AD","AD","AD","Old","Old","Old","Old","Old","Old","Old")))
# names <-  c("AD3","AD6","AD7", "AD8","AD9",
#              "Old1","Old2","Old3","Old4","Old5","Old6","Old7")
# rownames(sampleTable) = names
# colnames(txi.salmon$counts) <- names
# 
# dds_ADvsOld <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~ condition)
# dds_ADvsOld$condition <- relevel(dds_ADvsOld$condition, ref = "Old")
# dds_ADvsOld <- DESeq(dds_ADvsOld)
# resultsNames(dds_ADvsOld)
# 
# 
# ## PCA ##
# 
# vsd <- vst(dds_ADvsOld, blind = FALSE) # or varianceStabilizingTransformation instead of vst
# 
# plotPCA(vsd, intgroup="condition") +
#   geom_label(aes(label = name))
# 
# 
# ###### Old vs Young ######
# files = c('./Old/10-8A-Old.fastq.gz_counts/quant.sf', #1
#           './Old/11-10T-Old.fastq.gz_counts/quant.sf',#2
#           './Old/12-6A-Old.fastq.gz_counts/quant.sf', #3
#           './Old/13-11T-Old.fastq.gz_counts/quant.sf', #4
#           './Old/14-7A-Old.fastq.gz_counts/quant.sf', #5
#           './Old/15-13T-Old.fastq.gz_counts/quant.sf', #6
#           './Old/16-14T-Old.fastq.gz_counts/quant.sf', #7
#           './Young/2-12A-Young.fastq.gz_counts/quant.sf', #1
#           './Young/3-17T-Young.fastq.gz_counts/quant.sf', #2
#           './Young/4-13A-Young.fastq.gz_counts/quant.sf', #3
#           './Young/6-14A-Young.fastq.gz_counts/quant.sf') #5
# txi.salmon <- tximport(files = files, 
#                        type = "salmon",
#                        txOut = FALSE,
#                        tx2gene = tx2gene,
#                        ignoreAfterBar = T)
# 
# sampleTable <- data.frame(condition = factor(c("Old","Old","Old","Old","Old","Old","Old", "Young","Young","Young","Young")))
# names <-  c("Old1","Old2","Old3","Old4","Old5","Old6","Old7", 
#              "Young1","Young2","Young3","Young5")
# rownames(sampleTable) = names
# colnames(txi.salmon$counts) <- names
# 
# dds_OldvsYoung <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~ condition)
# dds_OldvsYoung$condition <- relevel(dds_OldvsYoung$condition, ref = "Young")
# dds_OldvsYoung <- DESeq(dds_OldvsYoung)
# resultsNames(dds_OldvsYoung)
# 
# ## PCA ##
# 
# vsd <- vst(dds_OldvsYoung, blind = FALSE) # or varianceStabilizingTransformation instead of vst
# 
# plotPCA(vsd, intgroup="condition") +
#   geom_label(aes(label = name))
# 
# 
# ###### AD vs Young ######
# 
# files = c('./AD/22-2T-AD.fastq.gz_counts/quant.sf', #3
#           './AD/25-5T-AD.fastq.gz_counts/quant.sf', #6
#           './AD/26-3A-AD.fastq.gz_counts/quant.sf', #7
#           './AD/27-5A-AD.fastq.gz_counts/quant.sf', #8
#           './AD/28-8T-AD.fastq.gz_counts/quant.sf', #9
#           './Young/2-12A-Young.fastq.gz_counts/quant.sf', #1
#           './Young/3-17T-Young.fastq.gz_counts/quant.sf', #2
#           './Young/4-13A-Young.fastq.gz_counts/quant.sf', #3
#           './Young/6-14A-Young.fastq.gz_counts/quant.sf') #5
# 
# txi.salmon <- tximport(files = files, 
#                        type = "salmon",
#                        txOut = FALSE,
#                        tx2gene = tx2gene,
#                        ignoreAfterBar = T)
# 
# 
# sampleTable <- data.frame(condition = factor(c("AD", "AD","AD","AD","AD", "Young","Young","Young","Young")))
# names <-  c("AD3","AD6","AD7", "AD8","AD9", 
#             "Young1","Young2","Young3","Young5")
# rownames(sampleTable) = names
# colnames(txi.salmon$counts) <- names
# 
# dds_ADvsYoung <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~ condition)
# dds_ADvsYoung$condition <- relevel(dds_ADvsYoung$condition, ref = "Young")
# dds_ADvsYoung <- DESeq(dds_ADvsYoung)
# resultsNames(dds_ADvsYoung)
# 
# ## PCA ##
# 
# vsd <- vst(dds_ADvsYoung, blind = FALSE) # or varianceStabilizingTransformation instead of vst
# 
# plotPCA(vsd, intgroup="condition") +
#   geom_label(aes(label = name))
# 
# 
# 
# ####### 2nd round selection: DGE SELECTING SAMPLES PER CONDITION #######
# 
# # files = c('./AD/22-2T-AD.fastq.gz_counts/quant.sf', #3
# #           './AD/25-5T-AD.fastq.gz_counts/quant.sf', #6
# #           './AD/26-3A-AD.fastq.gz_counts/quant.sf', #7
# #           './AD/28-8T-AD.fastq.gz_counts/quant.sf', #9
# #           './Old/13-11T-Old.fastq.gz_counts/quant.sf', #4
# #           './Old/14-7A-Old.fastq.gz_counts/quant.sf', #5
# #           './Old/15-13T-Old.fastq.gz_counts/quant.sf', #6
# #           './Old/16-14T-Old.fastq.gz_counts/quant.sf', #7
# #           './Young/2-12A-Young.fastq.gz_counts/quant.sf', #1
# #           './Young/3-17T-Young.fastq.gz_counts/quant.sf', #2
# #           './Young/4-13A-Young.fastq.gz_counts/quant.sf', #3
# #           './Young/6-14A-Young.fastq.gz_counts/quant.sf') #5
# 
# ## Yi's tx2gene
# 
# tx2gene = readRDS("D://ABBY.windows_surface_pro/DIFF-EXP/tx2gene.rds")
# setwd("D://ABBY.windows_surface_pro/FASTQ.files/")
# 
# ###### AD vs Old ######
# 
# files = c('./AD/22-2T-AD.fastq.gz_counts/quant.sf', #3
#           './AD/25-5T-AD.fastq.gz_counts/quant.sf', #6
#           './AD/26-3A-AD.fastq.gz_counts/quant.sf', #7
#           './AD/28-8T-AD.fastq.gz_counts/quant.sf', #9
#           './Old/13-11T-Old.fastq.gz_counts/quant.sf', #4
#           './Old/14-7A-Old.fastq.gz_counts/quant.sf', #5
#           './Old/15-13T-Old.fastq.gz_counts/quant.sf', #6
#           './Old/16-14T-Old.fastq.gz_counts/quant.sf') #7
# 
# txi.salmon <- tximport(files = files, 
#                        type = "salmon",
#                        txOut = FALSE,
#                        tx2gene = tx2gene,
#                        ignoreAfterBar = T)
# 
# sampleTable <- data.frame(condition = factor(c("AD", "AD","AD","AD","Old","Old","Old","Old")))
# names <-  c("AD3","AD6","AD7","AD9",
#             "Old4","Old5","Old6","Old7")
# rownames(sampleTable) = names
# colnames(txi.salmon$counts) <- names
# 
# dds_ADvsOld <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~ condition)
# dds_ADvsOld$condition <- relevel(dds_ADvsOld$condition, ref = "Old")
# dds_ADvsOld <- DESeq(dds_ADvsOld)
# resultsNames(dds_ADvsOld)
# 
# 
# ## PCA ##
# 
# vsd <- vst(dds_ADvsOld, blind = FALSE) # or varianceStabilizingTransformation instead of vst
# 
# plotPCA(vsd, intgroup="condition") +
#   geom_label(aes(label = name))
# 
# 
# ###### Old vs Young ######
# files = c('./Old/13-11T-Old.fastq.gz_counts/quant.sf', #4
#           './Old/14-7A-Old.fastq.gz_counts/quant.sf', #5
#           './Old/15-13T-Old.fastq.gz_counts/quant.sf', #6
#           './Old/16-14T-Old.fastq.gz_counts/quant.sf', #7
#           './Young/2-12A-Young.fastq.gz_counts/quant.sf', #1
#           './Young/3-17T-Young.fastq.gz_counts/quant.sf', #2
#           './Young/4-13A-Young.fastq.gz_counts/quant.sf', #3
#           './Young/6-14A-Young.fastq.gz_counts/quant.sf') #5
# txi.salmon <- tximport(files = files, 
#                        type = "salmon",
#                        txOut = FALSE,
#                        tx2gene = tx2gene,
#                        ignoreAfterBar = T)
# 
# sampleTable <- data.frame(condition = factor(c("Old","Old","Old","Old", "Young","Young","Young","Young")))
# names <-  c("Old4","Old5","Old6","Old7", 
#             "Young1","Young2","Young3","Young5")
# rownames(sampleTable) = names
# colnames(txi.salmon$counts) <- names
# 
# dds_OldvsYoung <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~ condition)
# dds_OldvsYoung$condition <- relevel(dds_OldvsYoung$condition, ref = "Young")
# dds_OldvsYoung <- DESeq(dds_OldvsYoung)
# resultsNames(dds_OldvsYoung)
# 
# ## PCA ##
# 
# vsd <- vst(dds_OldvsYoung, blind = FALSE) # or varianceStabilizingTransformation instead of vst
# 
# plotPCA(vsd, intgroup="condition") +
#   geom_label(aes(label = name))
# 
# 
# ###### AD vs Young ######
# 
# files = c('./AD/22-2T-AD.fastq.gz_counts/quant.sf', #3
#           './AD/25-5T-AD.fastq.gz_counts/quant.sf', #6
#           './AD/26-3A-AD.fastq.gz_counts/quant.sf', #7
#           './AD/28-8T-AD.fastq.gz_counts/quant.sf', #9
#           './Young/2-12A-Young.fastq.gz_counts/quant.sf', #1
#           './Young/3-17T-Young.fastq.gz_counts/quant.sf', #2
#           './Young/4-13A-Young.fastq.gz_counts/quant.sf', #3
#           './Young/6-14A-Young.fastq.gz_counts/quant.sf') #5
# 
# txi.salmon <- tximport(files = files, 
#                        type = "salmon",
#                        txOut = FALSE,
#                        tx2gene = tx2gene,
#                        ignoreAfterBar = T)
# 
# 
# sampleTable <- data.frame(condition = factor(c("AD","AD","AD","AD", "Young","Young","Young","Young")))
# names <-  c("AD3","AD6","AD7", "AD9", 
#             "Young1","Young2","Young3","Young5")
# rownames(sampleTable) = names
# colnames(txi.salmon$counts) <- names
# 
# dds_ADvsYoung <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~ condition)
# dds_ADvsYoung$condition <- relevel(dds_ADvsYoung$condition, ref = "Young")
# dds_ADvsYoung <- DESeq(dds_ADvsYoung)
# resultsNames(dds_ADvsYoung)
# 
# ## PCA ##
# 
# vsd <- vst(dds_ADvsYoung, blind = FALSE) # or varianceStabilizingTransformation instead of vst
# 
# plotPCA(vsd, intgroup="condition") +
#   geom_label(aes(label = name))
# 
# 
# 
# #### 3 per condition ####
# 
# 
# ## Yi's tx2gene
# 
tx2gene = readRDS("D://ABBY.windows_surface_pro/DIFF-EXP/tx2gene.rds")
setwd("D://ABBY.windows_surface_pro/FASTQ.files/")
# 
# ###### AD vs Old ######
# 
# files = c('./AD/22-2T-AD.fastq.gz_counts/quant.sf', #3
#           './AD/26-3A-AD.fastq.gz_counts/quant.sf', #7
#           './AD/28-8T-AD.fastq.gz_counts/quant.sf', #9
#           './Old/13-11T-Old.fastq.gz_counts/quant.sf', #4
#           './Old/14-7A-Old.fastq.gz_counts/quant.sf', #5
#           './Old/15-13T-Old.fastq.gz_counts/quant.sf') #6 
# 
# txi.salmon <- tximport(files = files, 
#                        type = "salmon",
#                        txOut = FALSE,
#                        tx2gene = tx2gene,
#                        ignoreAfterBar = T)
# 
# sampleTable <- data.frame(condition = factor(c("AD", "AD","AD","Old","Old","Old")))
# names <-  c("AD3","AD7","AD9",
#             "Old4","Old5","Old6")
# rownames(sampleTable) = names
# colnames(txi.salmon$counts) <- names
# 
# dds_ADvsOld <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~ condition)
# dds_ADvsOld$condition <- relevel(dds_ADvsOld$condition, ref = "Old")
# dds_ADvsOld <- DESeq(dds_ADvsOld)
# resultsNames(dds_ADvsOld)
# 
# 
# ## PCA ##
# 
# vsd <- vst(dds_ADvsOld, blind = FALSE) # or varianceStabilizingTransformation instead of vst
# 
# plotPCA(vsd, intgroup="condition") +
#   geom_label(aes(label = name))
# 
# 
# ###### Old vs Young ######
# files = c('./Old/13-11T-Old.fastq.gz_counts/quant.sf', #4
#           './Old/14-7A-Old.fastq.gz_counts/quant.sf', #5
#           './Old/15-13T-Old.fastq.gz_counts/quant.sf', #6
#           './Young/2-12A-Young.fastq.gz_counts/quant.sf', #1
#           './Young/4-13A-Young.fastq.gz_counts/quant.sf', #3
#           './Young/6-14A-Young.fastq.gz_counts/quant.sf') #5
# txi.salmon <- tximport(files = files, 
#                        type = "salmon",
#                        txOut = FALSE,
#                        tx2gene = tx2gene,
#                        ignoreAfterBar = T)
# 
# sampleTable <- data.frame(condition = factor(c("Old","Old","Old","Young","Young","Young")))
# names <-  c("Old4","Old5","Old6", 
#             "Young1","Young3","Young5")
# rownames(sampleTable) = names
# colnames(txi.salmon$counts) <- names
# 
# dds_OldvsYoung <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~ condition)
# dds_OldvsYoung$condition <- relevel(dds_OldvsYoung$condition, ref = "Young")
# dds_OldvsYoung <- DESeq(dds_OldvsYoung)
# resultsNames(dds_OldvsYoung)
# 
# ## PCA ##
# 
# vsd <- vst(dds_OldvsYoung, blind = FALSE) # or varianceStabilizingTransformation instead of vst
# 
# plotPCA(vsd, intgroup="condition") +
#   geom_label(aes(label = name))
# 
# 
# ###### AD vs Young ######
# 
# files = c('./AD/22-2T-AD.fastq.gz_counts/quant.sf', #3
#           './AD/26-3A-AD.fastq.gz_counts/quant.sf', #7
#           './AD/28-8T-AD.fastq.gz_counts/quant.sf', #9
#           './Young/2-12A-Young.fastq.gz_counts/quant.sf', #1
#           './Young/4-13A-Young.fastq.gz_counts/quant.sf', #3
#           './Young/6-14A-Young.fastq.gz_counts/quant.sf') #5
# 
# txi.salmon <- tximport(files = files, 
#                        type = "salmon",
#                        txOut = FALSE,
#                        tx2gene = tx2gene,
#                        ignoreAfterBar = T)
# 
# 
# sampleTable <- data.frame(condition = factor(c("AD","AD","AD","Young","Young","Young")))
# names <-  c("AD3","AD7", "AD9", 
#             "Young1","Young3","Young5")
# rownames(sampleTable) = names
# colnames(txi.salmon$counts) <- names
# 
# dds_ADvsYoung <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~ condition)
# dds_ADvsYoung$condition <- relevel(dds_ADvsYoung$condition, ref = "Young")
# dds_ADvsYoung <- DESeq(dds_ADvsYoung)
# resultsNames(dds_ADvsYoung)
# 
# ## PCA ##
# 
# vsd <- vst(dds_ADvsYoung, blind = FALSE) # or varianceStabilizingTransformation instead of vst
# 
# plotPCA(vsd, intgroup="condition") +
#   geom_label(aes(label = name))


#### AD 3,6,9 

#### 3 per condition ####


## Yi's tx2gene

tx2gene = readRDS("D://ABBY.windows_surface_pro/DIFF-EXP/tx2gene.rds")
setwd("D://ABBY.windows_surface_pro/FASTQ.files/")

tx2gene = readRDS("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/DIFF-EXP/tx2gene.rds")
setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/FASTQ.files/")

###### AD vs Old ######

files = c('./AD/22-2T-AD.fastq.gz_counts/quant.sf', #3
          './AD/25-5T-AD.fastq.gz_counts/quant.sf', #6
          './AD/28-8T-AD.fastq.gz_counts/quant.sf', #9
          './Old/13-11T-Old.fastq.gz_counts/quant.sf', #4
          './Old/14-7A-Old.fastq.gz_counts/quant.sf', #5
          './Old/15-13T-Old.fastq.gz_counts/quant.sf') #6 

txi.salmon <- tximport(files = files, 
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)

sampleTable <- data.frame(condition = factor(c("AD", "AD","AD","Old","Old","Old")))
names <-  c("AD3","AD6","AD9",
            "Old4","Old5","Old6")
rownames(sampleTable) = names
colnames(txi.salmon$counts) <- names

dds_ADvsOld <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~ condition)
dds_ADvsOld$condition <- relevel(dds_ADvsOld$condition, ref = "Old")
dds_ADvsOld <- DESeq(dds_ADvsOld)
resultsNames(dds_ADvsOld)
res_ADvsOld <- lfcShrink(dds_ADvsOld, coef="condition_AD_vs_Old", type="apeglm")
summary(res_ADvsOld)
sum(res_ADvsOld$padj < 0.05, na.rm = TRUE)

## saving ##

setwd("D://ABBY.windows_surface_pro/differentialexp_update//")
saveRDS(res_ADvsOld, file = "./DESEQres_pairwise_ADvsOld(3each-bulk).rds")

## PCA ##

vsd <- vst(dds_ADvsOld, blind = FALSE) # or varianceStabilizingTransformation instead of vst

# PCA_ADOld<- plotPCA(vsd, intgroup="condition") +
#   geom_point(aes(color = condition), size = 1) +  # Use points (dots) for labels
#   geom_text(aes(label = name, fontface="bold"), colour = "black", nudge_y = 1.4, nudge_x = 0.35, size = 4) +
#   scale_colour_manual(values = c("AD" = "slategray3", "Old"="plum")) 

PCA_ADOld <- plotPCA(vsd, intgroup="condition") +
  geom_point(shape = 21, aes(fill = condition, colour=condition), size = 4) +  # Use points (dots) for labels
  geom_text(aes(label = name, fontface="bold"), colour = "black", nudge_y = 1.35, nudge_x = 0.35, size = 4) +
  scale_fill_manual(values = c("AD" = "orangered", "Old" = "darkslateblue")) +
  scale_colour_manual(values = c("darkslateblue", "orangered")) +
  labs(fill="", colour="") +
  theme(legend.position = "top")
PCA_ADOld 

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}

PCA_ADOld <- edit_plots(PCA_ADOld)
PCA_ADOld

setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp_update/")
saveRDS(PCA_ADOld, "./PCA_ADOld_3each_b&w.rds")


setwd("D://ABBY.windows_surface_pro/diffexp_results/")
saveRDS(PCA_ADOld, "./PCA_ADOld_3each_11sept.rds")

# plotPCA(vsd, intgroup="condition") +
#   geom_label(aes(label = name))


###### Old vs Young ######
files = c('./Old/13-11T-Old.fastq.gz_counts/quant.sf', #4
          './Old/14-7A-Old.fastq.gz_counts/quant.sf', #5
          './Old/15-13T-Old.fastq.gz_counts/quant.sf', #6
          './Young/2-12A-Young.fastq.gz_counts/quant.sf', #1
          './Young/4-13A-Young.fastq.gz_counts/quant.sf', #3
          './Young/6-14A-Young.fastq.gz_counts/quant.sf') #5
txi.salmon <- tximport(files = files, 
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)

sampleTable <- data.frame(condition = factor(c("Old","Old","Old","Young","Young","Young")))
names <-  c("Old4","Old5","Old6", 
            "Young1","Young3","Young5")
rownames(sampleTable) = names
colnames(txi.salmon$counts) <- names

dds_OldvsYoung <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~ condition)
dds_OldvsYoung$condition <- relevel(dds_OldvsYoung$condition, ref = "Young")
dds_OldvsYoung <- DESeq(dds_OldvsYoung)
resultsNames(dds_OldvsYoung)
res_OldvsYoung <- lfcShrink(dds_OldvsYoung, coef="condition_Old_vs_Young", type="apeglm")
summary(res_OldvsYoung)
sum(res_OldvsYoung$padj < 0.05, na.rm = TRUE)

## saving ##

setwd("D://ABBY.windows_surface_pro/differentialexp_update//")
saveRDS(res_OldvsYoung, file = "./DESEQres_pairwise_OldvsYoung(3each-bulk).rds")

## PCA ##

vsd <- vst(dds_OldvsYoung, blind = FALSE) # or varianceStabilizingTransformation instead of vst

# PCA_OldYoung<- plotPCA(vsd, intgroup="condition") +
#   geom_point(aes(color = condition), size = 1) +  # Use points (dots) for labels
#   geom_text(aes(label = name, fontface="bold"), colour = "black", nudge_y = 1.2, nudge_x = 0.35, size = 4) +
#   scale_colour_manual(values = c("Old"="plum", "Young"="orange")) 

PCA_OldYoung <- plotPCA(vsd, intgroup="condition") +
  geom_point(shape = 21, aes(fill = condition, colour=condition), size = 4) +  # Use points (dots) for labels
  geom_text(aes(label = name, fontface="bold"), colour = "black", nudge_y = 1.35, nudge_x = 0.35, size = 4) +
  scale_fill_manual(values = c("Old" = "darkslateblue", "Young" = "gold")) +
  scale_colour_manual(values = c("gold", "darkslateblue")) +
  labs(fill="", colour="") +
  theme(legend.position = "top")
PCA_OldYoung 


edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}

PCA_OldYoung <- edit_plots(PCA_OldYoung)
PCA_OldYoung

setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp_update/")
saveRDS(PCA_OldYoung, "./PCA_OldYoung_3each_b&w.rds")


setwd("D://ABBY.windows_surface_pro/diffexp_results/")
saveRDS(PCA_OldYoung, "./PCA_OldYoung_3each_11sept.rds")

# plotPCA(vsd, intgroup="condition") +
#   geom_label(aes(label = name))


###### AD vs Young ######

files = c('./AD/22-2T-AD.fastq.gz_counts/quant.sf', #3
          './AD/25-5T-AD.fastq.gz_counts/quant.sf', #6
          './AD/28-8T-AD.fastq.gz_counts/quant.sf', #9
          './Young/2-12A-Young.fastq.gz_counts/quant.sf', #1
          './Young/4-13A-Young.fastq.gz_counts/quant.sf', #3
          './Young/6-14A-Young.fastq.gz_counts/quant.sf') #5

txi.salmon <- tximport(files = files, 
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)


sampleTable <- data.frame(condition = factor(c("AD","AD","AD","Young","Young","Young")))
names <-  c("AD3","AD6", "AD9", 
            "Young1","Young3","Young5")
rownames(sampleTable) = names
colnames(txi.salmon$counts) <- names

dds_ADvsYoung <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~ condition)
dds_ADvsYoung$condition <- relevel(dds_ADvsYoung$condition, ref = "Young")
dds_ADvsYoung <- DESeq(dds_ADvsYoung)
resultsNames(dds_ADvsYoung)
res_ADvsYoung <- lfcShrink(dds_ADvsYoung, coef="condition_AD_vs_Young", type="apeglm")
summary(res_ADvsYoung)
sum(res_ADvsYoung$padj < 0.05, na.rm = TRUE)

## saving ##

setwd("D://ABBY.windows_surface_pro/differentialexp_update//")
saveRDS(res_ADvsYoung, file = "./DESEQres_pairwise_ADvsYoung(3each-bulk).rds")

## PCA ##

vsd <- vst(dds_ADvsYoung, blind = FALSE) # or varianceStabilizingTransformation instead of vst

# PCA_ADYoung<- plotPCA(vsd, intgroup="condition") +
#   geom_point(aes(color = condition), size = 1) +  # Use points (dots) for labels
#   geom_text(aes(label = name, fontface="bold"), colour = "black", nudge_y = 1.7, nudge_x = 0.35, size = 4) +
#   scale_colour_manual(values = c("AD"="slategray3", "Young"="orange")) 

PCA_ADYoung <- plotPCA(vsd, intgroup="condition") +
  geom_point(shape = 21, aes(fill = condition, colour=condition), size = 4) +  # Use points (dots) for labels
  geom_text(aes(label = name, fontface="bold"), colour = "black", nudge_y = 1.35, nudge_x = 0.35, size = 4) +
  scale_fill_manual(values = c("AD" = "orangered", "Young" = "gold")) +
  scale_colour_manual(values = c("gold", "orangered")) +
  labs(fill="", colour="") +
  theme(legend.position = "top")
PCA_ADYoung

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}

PCA_ADYoung <- edit_plots(PCA_ADYoung)
PCA_ADYoung

setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp_update/")
saveRDS(PCA_ADYoung, "./PCA_ADYoung_3each_b&w.rds")

setwd("D://ABBY.windows_surface_pro/diffexp_results/")
saveRDS(PCA_ADYoung, "./PCA_ADYoung_3each_11sept.rds")

# plotPCA(vsd, intgroup="condition") +
#   geom_label(aes(label = name))



###### Looking at DESeq2 results ######

setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp/")
setwd("D://ABBY.windows_surface_pro/differentialexp_update//")
## AD vs Old ## 

res_ADvsOld <- readRDS("./DESEQres_pairwise_ADvsOld(3each-bulk).rds")
ADvsOld_df <- data.frame(Gene = res_ADvsOld@rownames, 
                         LFC = res_ADvsOld@listData[["log2FoldChange"]],
                         pvalue = res_ADvsOld@listData[["pvalue"]],
                         padj = res_ADvsOld@listData[["padj"]],
                         basemean = res_ADvsOld@listData[["baseMean"]])

base_mean_values = res_ADvsOld@listData[["baseMean"]]
hist(log10(base_mean_values), 
     breaks = 250,  # Number of bins
     main = "Base Mean Distribution ADvsOld",
     xlab = "log10(Base Mean)",
     ylab = "Frequency")

filtered_ADvsOld_basemean <- ADvsOld_df[ADvsOld_df$basemean > 10, ]
filtered_ADvsOld_basemean <- filtered_ADvsOld_basemean[filtered_ADvsOld_basemean$padj < 0.05, ]
filtered_ADvsOld_LFC1 <- filtered_ADvsOld_basemean[filtered_ADvsOld_basemean$LFC > 1 | filtered_ADvsOld_basemean$LFC < -1, ]
filtered_ADvsOld_LFC1 <- na.omit(filtered_ADvsOld_LFC1)
saveRDS(filtered_ADvsOld_LFC1, "./dataframe_filtered_ADvsOld_LFC1_basemean10_padj0.05.rds")

## Old vs Young
res_OldvsYoung <- readRDS("./DESEQres_pairwise_OldvsYoung(3each-bulk).rds")
OldvsYoung_df <- data.frame(Gene = res_OldvsYoung@rownames, 
                            LFC = res_OldvsYoung@listData[["log2FoldChange"]],
                            pvalue = res_OldvsYoung@listData[["pvalue"]],
                            padj = res_OldvsYoung@listData[["padj"]],
                            basemean = res_OldvsYoung@listData[["baseMean"]])

base_mean_values = res_OldvsYoung@listData[["baseMean"]]
hist(log10(base_mean_values), 
     breaks = 250,  # Number of bins
     main = "Base Mean Distribution OldvsYoung",
     xlab = "log10(Base Mean)",
     ylab = "Frequency")


filtered_OldvsYoung_basemean <- OldvsYoung_df[OldvsYoung_df$basemean > 10, ]
filtered_OldvsYoung_basemean <- filtered_OldvsYoung_basemean[filtered_OldvsYoung_basemean$padj < 0.05, ]
filtered_OldvsYoung_LFC1 <- filtered_OldvsYoung_basemean[filtered_OldvsYoung_basemean$LFC > 1 | filtered_OldvsYoung_basemean$LFC < -1, ]
filtered_OldvsYoung_LFC1 <- na.omit(filtered_OldvsYoung_LFC1)
saveRDS(filtered_OldvsYoung_LFC1, "./dataframe_filtered_OldvsYoung_LFC1_basemean10_padj0.05.rds")


## AD vs YOUNG
res_ADvsYoung <- readRDS("./DESEQres_pairwise_ADvsYoung(3each-bulk).rds")
ADvsYoung_df <- data.frame(Gene = res_ADvsYoung@rownames, 
                           LFC = res_ADvsYoung@listData[["log2FoldChange"]],
                           pvalue = res_ADvsYoung@listData[["pvalue"]],
                           padj = res_ADvsYoung@listData[["padj"]],
                           basemean = res_ADvsYoung@listData[["baseMean"]])

base_mean_values = res_ADvsYoung@listData[["baseMean"]]
hist(log10(base_mean_values), 
     breaks = 250,  # Number of bins
     main = "Base Mean Distribution ADvsYoung",
     xlab = "log10(Base Mean)",
     ylab = "Frequency")

filtered_ADvsYoung_basemean <- ADvsYoung_df[ADvsYoung_df$basemean > 10, ]
filtered_ADvsYoung_basemean <- filtered_ADvsYoung_basemean[filtered_ADvsYoung_basemean$padj < 0.05, ]
filtered_ADvsYoung_LFC1 <- filtered_ADvsYoung_basemean[filtered_ADvsYoung_basemean$LFC > 1 | filtered_ADvsYoung_basemean$LFC < -1, ]
filtered_ADvsYoung_LFC1 <- na.omit(filtered_ADvsYoung_LFC1)
saveRDS(filtered_ADvsYoung_LFC1, "./dataframe_filtered_ADvsYoung_LFC1_basemean10_padj0.05.rds")



### filtered DESeq2 genes ###

setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp_update/")
setwd("D://ABBY.windows_surface_pro/differentialexp_update/")
ADvsOld <- readRDS("./dataframe_filtered_ADvsOld_LFC1_basemean10_padj0.05.rds")
OldvsYoung <- readRDS("./dataframe_filtered_OldvsYoung_LFC1_basemean10_padj0.05.rds")
ADvsYoung <- readRDS("./dataframe_filtered_ADvsYoung_LFC1_basemean10_padj0.05.rds")

## AD vs OLD - UP- & DOWN-REGULATED
ADvsOld_up <- ADvsOld[ADvsOld$LFC > 0, ]
ADvsOld_down <- ADvsOld[ADvsOld$LFC < 0, ]
saveRDS(ADvsOld_up, "./ADvsOld_up.rds")
saveRDS(ADvsOld_down, "./ADvsOld_down.rds")

## OLD vs YOUNG - UP- & DOWN-REGULATED
OldvsYoung_up <- OldvsYoung[OldvsYoung$LFC > 0, ]
OldvsYoung_down <- OldvsYoung[OldvsYoung$LFC < 0, ]
saveRDS(OldvsYoung_up, "./OldvsYoung_up.rds")
saveRDS(OldvsYoung_down, "./OldvsYoung_down.rds")

## AD vs YOUNG - UP- & DOWN-REGULATED
ADvsYoung_up <- ADvsYoung[ADvsYoung$LFC > 0, ]
ADvsYoung_down <- ADvsYoung[ADvsYoung$LFC < 0, ]
saveRDS(ADvsYoung_up, "./ADvsYoung_up.rds")
saveRDS(ADvsYoung_down, "./ADvsYoung_down.rds")



ADOld_genes_up <- ADvsOld_up$Gene
ADOld_genes_down <- ADvsOld_down$Gene

OldYoung_genes_up <- OldvsYoung_up$Gene
OldYoung_genes_down <- OldvsYoung_down$Gene

ADYoung_genes_up <- ADvsYoung_up$Gene
ADYoung_genes_down <- ADvsYoung_down$Gene


## Venn diagrams ##
#### UPREGULATED ####


gene_lists_up <- list(
  ADvsOld = ADOld_genes_up,
  OldvsYoung = OldYoung_genes_up,
  ADvsYoung = ADYoung_genes_up)

venn_result_up <- venn.diagram(
  x = gene_lists_up,
  category.names = c("AD vs Old", "Old vs Young", "AD vs Young"),
  filename = NULL,
  cat.default.pos = "text",
  cat.pos = c(0, 0, 0),  # Place labels at default position (inside the sets)
  cat.dist = c(-0.23,-0.235,0.305),  # Distance of labels from the sets
  cat.cex = 1.2,   # Label text size
  ext.text = TRUE,
  label.col = "black",
  fill = c("mediumpurple", "lightgreen", "orange"),  # Colors for the sets
  alpha = 0.5,
  col = c("mediumpurple", "lightgreen", "orange"),
  cat.fontface = "bold",
  fontface = "bold"
)

grid.newpage()
grid.draw(venn_result_up)

saveRDS(venn_result_up, "D://ABBY.windows_surface_pro/diffexp_results/venn_up_11sept.rds")

#### DOWNREGULATED ####
gene_lists_down <- list(
  ADvsOld = ADOld_genes_down,
  OldvsYoung = OldYoung_genes_down,
  ADvsYoung = ADYoung_genes_down)

venn_result_down <- venn.diagram(
  x = gene_lists_down,
  category.names = c("AD vs Old", "Old vs Young", "AD vs Young"),
  filename = NULL,
  cat.default.pos = "text",
  cat.pos = c(0, 0, 0),  
  cat.dist = c(-0.23,-0.235,0.305),  # Distance of labels from the sets
  cat.cex = 1.2,   # Label text size
  ext.text = TRUE,
  label.col = "black",
  fill = c("mediumpurple", "lightgreen", "orange"),  # Colors for the sets
  alpha = 0.5,
  col = c("mediumpurple", "lightgreen", "orange"),
  cat.fontface = "bold",
  fontface = "bold" # Make category names bold
  #fontfamily = "sans"  # Font family (adjust as needed)
)

grid.newpage()
grid.draw(venn_result_down)

saveRDS(venn_result_down, "D://ABBY.windows_surface_pro/diffexp_results/venn_down_11sept.rds")


multi_panel_venn <- ggarrange(venn_result_up, venn_result_down, ncol = 2,
                                    labels = "AUTO",
                                    label.x = 0,
                                    align = "hv",
                                    vjust = 4.5)
multi_panel_venn

?venn.diagram
?ggarrange

saveRDS(multi_panel_venn, "./multi_venn_8sept.rds")

###
# install.packages("ggvenn")
# library(ggvenn)
# 
# ggvenn(data=gene_lists_down,
#        fill_color = c("deeppink1", "hotpink1", "lightpink"),
#        show_percentage = FALSE,
#        set_name_size = 6,
#        text_size = 5
#        )
# ?ggvenn



# COMMON UPREGULATED
common_ADOld_ADYoung_up <- intersect(ADOld_genes_up, ADYoung_genes_up)
common_OldYoung_ADYoung_up <- intersect(OldYoung_genes_up, ADYoung_genes_up)
common_ADOld_OldYoung_up <- intersect(ADOld_genes_up, OldYoung_genes_up)

saveRDS(common_ADOld_ADYoung_up, "./common_ADOld_ADYoung_up.rds")
saveRDS(common_OldYoung_ADYoung_up, "./common_OldYoung_ADYoung_up.rds")
saveRDS(common_ADOld_OldYoung_up, "./common_ADOld_OldYoung_up.rds")
## DEGs IN OLD/YOUNG+AD/YOUNG TO EXCLUDE FROM DEGs IN AD/OLD+AD/YOUNG UPregulated
genes_exclude_up <- intersect(common_ADOld_ADYoung_up, common_OldYoung_ADYoung_up)
common_OldYoung_ADYoung_up %in% common_ADOld_ADYoung_up

# COMMON DOWNREGULATED
common_ADOld_ADYoung_down <- intersect(ADOld_genes_down, ADYoung_genes_down)
common_OldYoung_ADYoung_down <- intersect(OldYoung_genes_down, ADYoung_genes_down)
common_ADOld_OldYoung_down <- intersect(ADOld_genes_down, OldYoung_genes_down)

saveRDS(common_ADOld_ADYoung_down, "./common_ADOld_ADYoung_down.rds")
saveRDS(common_OldYoung_ADYoung_down, "./common_OldYoung_ADYoung_down.rds")
saveRDS(common_ADOld_OldYoung_down, "./common_ADOld_OldYoung_down.rds")


## DEGs IN OLD/YOUNG+AD/YOUNG TO EXCLUDE FROM DEGs IN AD/OLD+AD/YOUNG DOWNregulated
genes_exclude_down <- intersect(common_ADOld_ADYoung_down, common_OldYoung_ADYoung_down)

common_OldYoung_ADYoung_down %in% common_ADOld_ADYoung_down

saveRDS(genes_exclude_up, "./genes_exclude_up.rds")
saveRDS(genes_exclude_down, "./genes_exclude_down.rds")





###### PLOT TOTAL GENE COUNTS PER SAMPLE (selected samples) ######

#### boxlot ####

tx2gene = readRDS("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/DIFF-EXP/tx2gene.rds")
setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/FASTQ.files/")

files = c('./AD/22-2T-AD.fastq.gz_counts/quant.sf', #3
          './AD/25-5T-AD.fastq.gz_counts/quant.sf', #6
          './AD/28-8T-AD.fastq.gz_counts/quant.sf', #9
          './Old/13-11T-Old.fastq.gz_counts/quant.sf', #4
          './Old/14-7A-Old.fastq.gz_counts/quant.sf', #5
          './Old/15-13T-Old.fastq.gz_counts/quant.sf', #6
          './Young/2-12A-Young.fastq.gz_counts/quant.sf', #1
          './Young/4-13A-Young.fastq.gz_counts/quant.sf', #3
          './Young/6-14A-Young.fastq.gz_counts/quant.sf') #5

txi.salmon <- tximport(files = files, 
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)

sampleTable <- data.frame(condition = factor(c("AD", "AD","AD","Old","Old","Old","Young","Young","Young")))
names <-  c("AD3","AD6","AD9",
            "Old4","Old5","Old6",
            "Young1","Young3","Young5")
rownames(sampleTable) = names
colnames(txi.salmon$counts) <- names
colnames(txi.salmon$abundance) <- names

counts_df <- as.data.frame(txi.salmon$counts)

# Calculate total gene counts per sample
totals <- colSums(counts_df)
samples <- names(totals)
condition <- c("AD","AD","AD",
               "Old", "Old","Old",
               "Young", "Young","Young")

# Convert the result to a data frame
total_counts_df <- data.frame(Sample = samples, TotalCount = totals, condition)

boxplot_counts <- ggplot(total_counts_df, aes(x = condition, y = TotalCount)) +
  geom_boxplot(fill="slategrey") +
  labs(x = "Sample Group",
       y = "Total Gene Count") +
  scale_fill_discrete(name = "Sample") +
  theme_minimal() +
  geom_signif(comparisons = list(c("AD","Old"),c("Old", "Young"), c("AD", "Young")), 
              map_signif_level=TRUE)

#Darisia Edit_plot
edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          # panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
} 
boxplot_counts <- edit_plots(boxplot_counts)
boxplot_counts

## adding individual samples onto plot
boxplot_counts <- boxplot_counts +
  geom_point(data = total_counts_df, aes(x = condition, y = TotalCount, colour = condition), ## change shape: shape = condition
              size = 3) +
  geom_text_repel(label= total_counts_df$Sample) +
  scale_colour_manual(values = c("AD" = "slategray3", "Old" = "plum", "Young" = "orange"))

boxplot_counts


### PLOTTING TOTAL DEGs PER PAIRWISE COMPARISON ###
### filtered DESeq2 genes ###
setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp/")

ADvsOld <- readRDS("./dataframe_filtered_ADvsOld_LFC1_basemean10_padj0.05.rds")
OldvsYoung <- readRDS("./dataframe_filtered_OldvsYoung_LFC1_basemean10_padj0.05.rds")
ADvsYoung <- readRDS("./dataframe_filtered_ADvsYoung_LFC1_basemean10_padj0.05.rds")

nrow(ADvsOld)
nrow(OldvsYoung)
nrow(ADvsYoung)

totals <- c(nrow(ADvsOld), nrow(OldvsYoung), nrow(ADvsYoung))
comparisons <- factor(c("ADvsOld", "OldvsYoung", "ADvsYoung"), levels = c("ADvsOld", "OldvsYoung", "ADvsYoung")) ### setting levels keeps order (not alphabetically)
total_counts_df <- data.frame(Comaprison = comparisons, TotalCount = totals)
levels(comparisons)
comp <- c("ADvsOld", "OldvsYoung", "ADvsYoung")
plot_counts <- ggplot(total_counts_df, aes(x = Comaprison, y = TotalCount, fill = comp)) +
  geom_bar(stat = "identity") +
  labs(title = "Total DEGs per Comparison", x = "Comparison", y = "Total DEGs") +
  scale_fill_manual(values = c("ADvsOld" = "slategray4", "OldvsYoung" = "slategray3", "ADvsYoung" = "slategray1"))

plot_counts


### PLOTTING UP- AND DOWN-REGULATED DEGs PER PAIRWISE COMPARISON ###
### filtered DESeq2 genes ###
setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp_update/")

ADvsOld <- readRDS("./dataframe_filtered_ADvsOld_LFC1_basemean10_padj0.05.rds")
OldvsYoung <- readRDS("./dataframe_filtered_OldvsYoung_LFC1_basemean10_padj0.05.rds")
ADvsYoung <- readRDS("./dataframe_filtered_ADvsYoung_LFC1_basemean10_padj0.05.rds")

## AD vs OLD - UP- & DOWN-REGULATED
ADvsOld_up <- ADvsOld[ADvsOld$LFC > 0, ]
ADvsOld_down <- ADvsOld[ADvsOld$LFC < 0, ]
## OLD vs YOUNG - UP- & DOWN-REGULATED
OldvsYoung_up <- OldvsYoung[OldvsYoung$LFC > 0, ]
OldvsYoung_down <- OldvsYoung[OldvsYoung$LFC < 0, ]
## AD vs YOUNG - UP- & DOWN-REGULATED
ADvsYoung_up <- ADvsYoung[ADvsYoung$LFC > 0, ]
ADvsYoung_down <- ADvsYoung[ADvsYoung$LFC < 0, ]

totals <- c(nrow(ADvsOld_up), nrow(ADvsOld_down), nrow(OldvsYoung_up),nrow(OldvsYoung_down), nrow(ADvsYoung_up),nrow(ADvsYoung_down))
comparisons <- factor(c("ADvsOld", "ADvsOld", "OldvsYoung", "OldvsYoung", "ADvsYoung", "ADvsYoung"), levels = c("ADvsOld", "OldvsYoung", "ADvsYoung")) 
direction = factor(c("UP", "DOWN","UP", "DOWN","UP", "DOWN"), levels = c("UP", "DOWN"))
total_counts_df <- data.frame(Comaprison = comparisons, TotalCount = totals, Direction = direction)
levels(direction)


plot_counts <- ggplot(total_counts_df, aes(x = Comaprison, y = TotalCount, fill = Direction)) +
  geom_bar(aes(fill = direction), position = "dodge", stat="identity") +
  geom_text(aes(label = TotalCount), vjust = -0.5, position = position_dodge(width = 0.9)) + 
  labs(title = "Up- and Down-Regulated DEGs per Comparison", x = "Comparison", y = "Total DEGs") +
  scale_fill_manual(values = c("UP" = "springgreen", "DOWN" = "firebrick2"))
plot_counts

plot_counts <- edit_plots(plot_counts)
plot_counts

saveRDS(total_counts_df, "./dataframe_UP-DOWN_DEGs_percomparison_forplotting.rds")



### PLOTTING TOTAL DEGs PER PAIRWISE COMPARISON ###
### filtered DESeq2 genes ###
setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp_update/")

ADvsOld <- readRDS("./dataframe_filtered_ADvsOld_LFC1_basemean10_padj0.05.rds")
OldvsYoung <- readRDS("./dataframe_filtered_OldvsYoung_LFC1_basemean10_padj0.05.rds")
ADvsYoung <- readRDS("./dataframe_filtered_ADvsYoung_LFC1_basemean10_padj0.05.rds")

nrow(ADvsOld)
nrow(OldvsYoung)
nrow(ADvsYoung)

totals <- c(nrow(ADvsOld), nrow(OldvsYoung), nrow(ADvsYoung))
comparisons <- factor(c("ADvsOld", "OldvsYoung", "ADvsYoung"), levels = c("ADvsOld", "OldvsYoung", "ADvsYoung")) ### setting levels keeps order (not alphabetically)
total_counts_df <- data.frame(Comaprison = comparisons, TotalCount = totals)
levels(comparisons)
comp <- c("ADvsOld", "OldvsYoung", "ADvsYoung")
plot_counts <- ggplot(total_counts_df, aes(x = Comaprison, y = TotalCount, fill = comp)) +
  geom_bar(stat = "identity") +
  labs(title = "Total DEGs per Comparison", x = "Comparison", y = "Total DEGs") +
  scale_fill_manual(values = c("ADvsOld" = "slategray4", "OldvsYoung" = "slategray3", "ADvsYoung" = "slategray1"))

plot_counts

#Darisia Edit_plot
edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          # panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
} 
plot_counts <- edit_plots(plot_counts)
plot_counts




###### Looking at DESeq2 results ######
### label volcano plots with genes LFC > 4 or <-4 ###
setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp_update/")

## AD vs Old ## 
res_ADvsOld <- readRDS("./DESEQres_pairwise_ADvsOld(3each-bulk).rds")
ADvsOld_df <- data.frame(Gene = res_ADvsOld@rownames, 
                         LFC = res_ADvsOld@listData[["log2FoldChange"]],
                         pvalue = res_ADvsOld@listData[["pvalue"]],
                         padj = res_ADvsOld@listData[["padj"]],
                         basemean = res_ADvsOld@listData[["baseMean"]])

filtered_ADvsOld_basemean <- ADvsOld_df[ADvsOld_df$basemean > 10, ]
filtered_ADvsOld_basemean <- na.omit(filtered_ADvsOld_basemean)

## volcano plot no LFC cutoff
volcano <- ggplot(filtered_ADvsOld_basemean, aes(x = LFC, y = -log10(padj), color = filtered_ADvsOld_basemean$padj < 0.05 & abs(filtered_ADvsOld_basemean$LFC) > 1)) +
  geom_point() +
  scale_color_manual(values = c("slategrey", "springgreen")) +  # Color for non-significant and significant genes
  theme_minimal() +
  labs(
    x = "Log Fold Change (LFC)",
    y = "-log10(Adjusted p-value)",
    colour = "Differentially Expressed Genes"
  ) +
  geom_text_repel(label= ifelse(abs(filtered_ADvsOld_basemean$padj) < 0.05 & abs(filtered_ADvsOld_basemean$LFC) > 4, filtered_ADvsOld_basemean$Gene, NA))

sum(abs(filtered_ADvsOld_basemean$padj) < 0.05 & abs(filtered_ADvsOld_basemean$LFC) > 4, na.rm = TRUE)

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          # panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
} 
volcano <- edit_plots(volcano)
volcano

## Old vs Young
res_OldvsYoung <- readRDS("./DESEQres_pairwise_OldvsYoung(3each-bulk).rds")
OldvsYoung_df <- data.frame(Gene = res_OldvsYoung@rownames, 
                            LFC = res_OldvsYoung@listData[["log2FoldChange"]],
                            pvalue = res_OldvsYoung@listData[["pvalue"]],
                            padj = res_OldvsYoung@listData[["padj"]],
                            basemean = res_OldvsYoung@listData[["baseMean"]])
filtered_OldvsYoung_basemean <- OldvsYoung_df[OldvsYoung_df$basemean > 10, ]
filtered_OldvsYoung_basemean <- na.omit(filtered_OldvsYoung_basemean)

## volcano plot no LFC cutoff
volcano <- ggplot(filtered_OldvsYoung_basemean, aes(x = LFC, y = -log10(padj), color = filtered_OldvsYoung_basemean$padj < 0.05 & abs(filtered_OldvsYoung_basemean$LFC) > 1)) +
  geom_point() +
  scale_color_manual(values = c("slategrey", "springgreen")) +  # Color for non-significant and significant genes
  theme_minimal() +
  labs(
    x = "Log Fold Change (LFC)",
    y = "-log10(Adjusted p-value)",
    colour = "Differentially Expressed Genes"
  ) +
  geom_text_repel(label= ifelse(abs(filtered_OldvsYoung_basemean$padj) < 0.05 & abs(filtered_OldvsYoung_basemean$LFC) > 4, filtered_OldvsYoung_basemean$Gene, NA))

sum(abs(filtered_OldvsYoung_basemean$padj) < 0.05 & abs(filtered_OldvsYoung_basemean$LFC) > 4, na.rm = TRUE)
edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          # panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
} 
volcano <- edit_plots(volcano)
volcano

## AD vs Young
res_ADvsYoung <- readRDS("./DESEQres_pairwise_ADvsYoung(3each-bulk).rds")
ADvsYoung_df <- data.frame(Gene = res_ADvsYoung@rownames, 
                           LFC = res_ADvsYoung@listData[["log2FoldChange"]],
                           pvalue = res_ADvsYoung@listData[["pvalue"]],
                           padj = res_ADvsYoung@listData[["padj"]],
                           basemean = res_ADvsYoung@listData[["baseMean"]])
filtered_ADvsYoung_basemean <- ADvsYoung_df[ADvsYoung_df$basemean > 10, ]
filtered_ADvsYoung_basemean <- na.omit(filtered_ADvsYoung_basemean)

## volcano plot no LFC cutoff
volcano <- ggplot(filtered_ADvsYoung_basemean, aes(x = LFC, y = -log10(padj), color = filtered_ADvsYoung_basemean$padj < 0.05 & abs(filtered_ADvsYoung_basemean$LFC) > 1)) +
  geom_point() +
  scale_color_manual(values = c("slategrey", "springgreen")) +  # Color for non-significant and significant genes
  theme_minimal() +
  labs(
    x = "Log Fold Change (LFC)",
    y = "-log10(Adjusted p-value)",
    colour = "Differentially Expressed Genes"
  ) +
  geom_text_repel(label= ifelse(abs(filtered_ADvsYoung_basemean$padj) < 0.05 & abs(filtered_ADvsYoung_basemean$LFC) > 4, filtered_ADvsYoung_basemean$Gene, NA))

sum(abs(filtered_ADvsYoung_basemean$padj) < 0.05 & abs(filtered_ADvsYoung_basemean$LFC) > 4, na.rm = TRUE)

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          # panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
} 
volcano <- edit_plots(volcano)
volcano


##### CLUSTERPROFILER #####
## GENES i put in = all filters applied
## BACKGROUND/UNIVERSE = only base mean filter

library(clusterProfiler)
setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp_update/")
#setwd("D://ABBY.windows_surface_pro/differentialexp/")

### UNIVERSES ###
## AD vs OLD
res_ADvsOld <- readRDS("./DESEQres_pairwise_ADvsOld(3each-bulk).rds")
ADvsOld_df <- data.frame(Gene = res_ADvsOld@rownames, 
                         LFC = res_ADvsOld@listData[["log2FoldChange"]],
                         pvalue = res_ADvsOld@listData[["pvalue"]],
                         padj = res_ADvsOld@listData[["padj"]],
                         basemean = res_ADvsOld@listData[["baseMean"]])
filtered_ADvsOld_basemean <- ADvsOld_df[ADvsOld_df$basemean > 10, ]

## OLD vs YOUNG
res_OldvsYoung <- readRDS("./DESEQres_pairwise_OldvsYoung(3each-bulk).rds")
OldvsYoung_df <- data.frame(Gene = res_OldvsYoung@rownames, 
                            LFC = res_OldvsYoung@listData[["log2FoldChange"]],
                            pvalue = res_OldvsYoung@listData[["pvalue"]],
                            padj = res_OldvsYoung@listData[["padj"]],
                            basemean = res_OldvsYoung@listData[["baseMean"]])
filtered_OldvsYoung_basemean <- OldvsYoung_df[OldvsYoung_df$basemean > 10, ]

## AD vs YOUNG
res_ADvsYoung <- readRDS("./DESEQres_pairwise_ADvsYoung(3each-bulk).rds")
ADvsYoung_df <- data.frame(Gene = res_ADvsYoung@rownames, 
                           LFC = res_ADvsYoung@listData[["log2FoldChange"]],
                           pvalue = res_ADvsYoung@listData[["pvalue"]],
                           padj = res_ADvsYoung@listData[["padj"]],
                           basemean = res_ADvsYoung@listData[["baseMean"]])
filtered_ADvsYoung_basemean <- ADvsYoung_df[ADvsYoung_df$basemean > 10, ]


## INPUT GENES ##

### filtered DESeq2 genes ###
setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp_update/")

ADvsOld <- readRDS("./dataframe_filtered_ADvsOld_LFC1_basemean10_padj0.05.rds")
OldvsYoung <- readRDS("./dataframe_filtered_OldvsYoung_LFC1_basemean10_padj0.05.rds")
ADvsYoung <- readRDS("./dataframe_filtered_ADvsYoung_LFC1_basemean10_padj0.05.rds")


### REMOVING VERSION NUMBERS OF GENES
ADvsOld$Gene = gsub("(ENSG[0-9]+)\\.[0-9]+","\\1",ADvsOld$Gene)
filtered_ADvsOld_basemean$Gene = gsub("(ENSG[0-9]+)\\.[0-9]+","\\1",filtered_ADvsOld_basemean$Gene)
ADvsOld_genes <- ADvsOld$Gene

OldvsYoung$Gene = gsub("(ENSG[0-9]+)\\.[0-9]+","\\1",OldvsYoung$Gene)
filtered_OldvsYoung_basemean$Gene = gsub("(ENSG[0-9]+)\\.[0-9]+","\\1",filtered_OldvsYoung_basemean$Gene)
OldvsYoung_genes <- OldvsYoung$Gene

ADvsYoung$Gene = gsub("(ENSG[0-9]+)\\.[0-9]+","\\1",ADvsYoung$Gene)
filtered_ADvsYoung_basemean$Gene = gsub("(ENSG[0-9]+)\\.[0-9]+","\\1",filtered_ADvsYoung_basemean$Gene)
ADvsYoung_genes <- ADvsYoung$Gene

### to get BPs for AD Old up/down separately
## AD vs OLD - UP- & DOWN-REGULATED
ADvsOld_up <- ADvsOld[ADvsOld$LFC > 0, ]
ADvsOld_down <- ADvsOld[ADvsOld$LFC < 0, ]
## OLD vs YOUNG - UP- & DOWN-REGULATED
OldvsYoung_up <- OldvsYoung[OldvsYoung$LFC > 0, ]
OldvsYoung_down <- OldvsYoung[OldvsYoung$LFC < 0, ]
## AD vs YOUNG - UP- & DOWN-REGULATED
ADvsYoung_up <- ADvsYoung[ADvsYoung$LFC > 0, ]
ADvsYoung_down <- ADvsYoung[ADvsYoung$LFC < 0, ]

## genes
ADOld_genes_up <- ADvsOld_up$Gene
ADOld_genes_down <- ADvsOld_down$Gene

OldYoung_genes_up <- OldvsYoung_up$Gene
OldYoung_genes_down <- OldvsYoung_down$Gene

ADYoung_genes_up <- ADvsYoung_up$Gene
ADYoung_genes_down <- ADvsYoung_down$Gene

unique_ADOld_up <- setdiff(ADOld_genes_up, ADYoung_genes_up)
unique_ADOld_up <- setdiff(unique_ADOld_up, OldYoung_genes_up)

unique_ADOld_down <- setdiff(ADOld_genes_down, ADYoung_genes_down)
unique_ADOld_down <- setdiff(unique_ADOld_down, OldYoung_genes_down)

### AD vs Old ###
## upregulated BPs

## no background
GOenrich_ADOld <- enrichGO(unique_ADOld_up, 'org.Hs.eg.db', keyType = "ENSEMBL", ont = "BP" ,pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
## + background
GOenrich_ADOld2 <- enrichGO(unique_ADOld_up, 'org.Hs.eg.db', keyType = "ENSEMBL", ont = "BP" ,pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05,
                            readable = TRUE,universe = filtered_ADvsOld_basemean$Gene) ##filtered basemean not LFC or p value == what IS expressed in cell

saveRDS(GOenrich_ADOld, "./GOenrich_uniqueADOld_up_nobackground.rds")
saveRDS(GOenrich_ADOld2, "./GOenrich_uniqueADOld_up_+background.rds")

## downregulated BPs

## no background
GOenrich_ADOld3 <- enrichGO(unique_ADOld_down, 'org.Hs.eg.db', keyType = "ENSEMBL", ont = "BP" ,pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
## + background
GOenrich_ADOld4 <- enrichGO(unique_ADOld_down, 'org.Hs.eg.db', keyType = "ENSEMBL", ont = "BP" ,pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05,
                            readable = TRUE,universe = filtered_ADvsOld_basemean$Gene) ##filtered basemean not LFC or p value == what IS expressed in cell

saveRDS(GOenrich_ADOld3, "./GOenrich_uniqueADOld_down_nobackground.rds")
saveRDS(GOenrich_ADOld4, "./GOenrich_uniqueADOld_down_+background.rds")

### AD vs Old altogether  ###

## no background
GOenrich_ADOld <- enrichGO(ADvsOld_genes, 'org.Hs.eg.db', keyType = "ENSEMBL", ont = "BP" ,pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05,readable = TRUE)
## + background
GOenrich_ADOld2 <- enrichGO(ADvsOld_genes, 'org.Hs.eg.db', keyType = "ENSEMBL", ont = "BP" ,pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05,
                               readable = TRUE,universe = filtered_ADvsOld_basemean$Gene) ##filtered basemean not LFC or p value == what IS expressed in cell

saveRDS(GOenrich_ADOld, "./GOenrich_ADOld_all_nobackground.rds")
saveRDS(GOenrich_ADOld2, "./GOenrich_ADOld_all_+background.rds")

### Old vs Young ###

## no background
GOenrich_OldYoung <- enrichGO(OldvsYoung_genes, 'org.Hs.eg.db', keyType = "ENSEMBL", ont = "BP" ,pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05,readable = TRUE)
## + background
GOenrich_OldYoung2 <- enrichGO(OldvsYoung_genes, 'org.Hs.eg.db', keyType = "ENSEMBL", ont = "BP" ,pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05,
                               readable = TRUE,universe = filtered_OldvsYoung_basemean$Gene) ##filtered basemean not LFC or p value == what IS expressed in cell

saveRDS(GOenrich_OldYoung, "./GOenrich_OldYoung_nobackground.rds")
saveRDS(GOenrich_OldYoung2, "./GOenrich_OldYoung_+background.rds")

### AD vs Young ###

## no background
GOenrich_ADYoung <- enrichGO(ADvsYoung_genes, 'org.Hs.eg.db', keyType = "ENSEMBL", ont = "BP" ,pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05,readable = TRUE)
## + background
GOenrich_ADYoung2 <- enrichGO(ADvsYoung_genes, 'org.Hs.eg.db', keyType = "ENSEMBL", ont = "BP" ,pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05,
                              readable = TRUE,universe = filtered_ADvsYoung_basemean$Gene) ##filtered basemean not LFC or p value == what IS expressed in cell

saveRDS(GOenrich_ADYoung, "./GOenrich_ADYoung_nobackground.rds")
saveRDS(GOenrich_ADYoung2, "./GOenrich_ADYoung_+background.rds")


### PLOTTING CLUSTERPROFILER ###
setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp_update/")
enrich_ADOld_up <- readRDS("./GOenrich_uniqueADOld_up_+background.rds")
enrich_ADOld_down <- readRDS("./GOenrich_uniqueADOld_down_+background.rds")

enrich_ADOld <- readRDS("./GOenrich_ADOld_all_+background.rds")
enrich_OldYoung <- readRDS("./GOenrich_OldYoung_+background.rds")
enrich_ADYoung <- readRDS("./GOenrich_ADYoung_+background.rds")

## results dataframes
res_enrich_ADOld_up <- as.data.frame(enrich_ADOld_up@result)
res_enrich_ADOld_down <- as.data.frame(enrich_ADOld_down@result)

res_enrich_ADOld <- as.data.frame(enrich_ADOld@result)
res_enrich_OldYoung <- as.data.frame(enrich_OldYoung@result)
res_enrich_ADYoung <- as.data.frame(enrich_ADYoung@result)

## filtering top 20 BPs

### ordering by increasing p.adjust ###
# Column to order by
column_to_order <- "Count"
# Order the data frame based on the specified column
ordered_res_enrich_ADOld_up <- res_enrich_ADOld_up[order(res_enrich_ADOld_up[[column_to_order]], decreasing = TRUE), ]
ordered_res_enrich_ADOld_down <- res_enrich_ADOld_down[order(res_enrich_ADOld_down[[column_to_order]], decreasing = TRUE), ]

ordered_res_enrich_ADOld <- res_enrich_ADOld[order(res_enrich_ADOld[[column_to_order]], decreasing = TRUE), ]
ordered_res_enrich_OldYoung <- res_enrich_OldYoung[order(res_enrich_OldYoung[[column_to_order]], decreasing = TRUE), ]
ordered_res_enrich_ADYoung <- res_enrich_ADYoung[order(res_enrich_ADYoung[[column_to_order]], decreasing = TRUE), ]


ordered_res_enrich_ADOld_up <- ordered_res_enrich_ADOld_up[1:20, ]
ordered_res_enrich_ADOld_down <- ordered_res_enrich_ADOld_down[1:20, ]

ordered_res_enrich_ADOld <- ordered_res_enrich_ADOld[1:20, ]
ordered_res_enrich_OldYoung <- ordered_res_enrich_OldYoung[1:20, ]
ordered_res_enrich_ADYoung <- ordered_res_enrich_ADYoung[1:20, ]


df_enrich_ADOld_up <- data.frame(BP = ordered_res_enrich_ADOld_up$Description, padj = ordered_res_enrich_ADOld_up$p.adjust, counts = ordered_res_enrich_ADOld_up$Count)
df_enrich_ADOld_down <- data.frame(BP = ordered_res_enrich_ADOld_down$Description, padj = ordered_res_enrich_ADOld_down$p.adjust, counts = ordered_res_enrich_ADOld_down$Count)

df_enrich_ADOld <- data.frame(BP = ordered_res_enrich_ADOld$Description, padj = ordered_res_enrich_ADOld$p.adjust, counts = ordered_res_enrich_ADOld$Count)
df_enrich_OldYoung <- data.frame(BP = ordered_res_enrich_OldYoung$Description, padj = ordered_res_enrich_OldYoung$p.adjust, counts = ordered_res_enrich_OldYoung$Count)
df_enrich_ADYoung <- data.frame(BP = ordered_res_enrich_ADYoung$Description, padj = ordered_res_enrich_ADYoung$p.adjust, counts = ordered_res_enrich_ADYoung$Count)

### manual AD vs OLD
## up
plot_enrich_ADOld_up <- ggplot(df_enrich_ADOld_up, aes(x = counts, y = BP), horiz = TRUE) + # horiz = TRUE horizontal bars
  geom_bar(stat = "identity", fill = "hotpink") +  
  labs(title = "Significant Enriched Terms AD vs Old (UP)", x = "Count", y = "Biological Process")
edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          # panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}

plot_enrich_ADOld_up<- edit_plots(plot_enrich_ADOld_up)
plot_enrich_ADOld_up

## down
plot_enrich_ADOld_down <- ggplot(df_enrich_ADOld_down, aes(x = counts, y = BP), horiz = TRUE) + # horiz = TRUE horizontal bars
  geom_bar(stat = "identity", fill = "hotpink") +  
  labs(title = "Significant Enriched Terms AD vs Old (DOWN)", x = "Count", y = "Biological Process")

plot_enrich_ADOld_down <- edit_plots(plot_enrich_ADOld_down)
plot_enrich_ADOld_down

## altogether
plot_enrich_ADOld <- ggplot(df_enrich_ADOld, aes(x = counts, y = BP, horiz = TRUE)) + # horiz = TRUE horizontal bars
  geom_bar(stat = "identity", fill = "hotpink") +  
  labs(title = "Significant Enriched Terms AD vs Old (all)", x = "Count", y = "Biological Process")

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          # panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}
plot_enrich_ADOld <- edit_plots(plot_enrich_ADOld)
plot_enrich_ADOld

### manual OLD vs YOUNG
#library(stringr)
plot_enrich_OldYoung <- ggplot(df_enrich_OldYoung, aes(x = counts, y = BP, horiz = TRUE)) + # horiz = TRUE horizontal bars
  geom_bar(stat = "identity", fill = "hotpink") +
  labs(title = "Significant Enriched Terms Old vs Young", x = "Count", y = "Biological Process")

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          # panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}

plot_enrich_OldYoung <- edit_plots(plot_enrich_OldYoung)
plot_enrich_OldYoung

# y = str_wrap(BP, width = 50))

## manual AD vs YOUNG
plot_enrich_ADYoung <- ggplot(df_enrich_ADYoung, aes(x = counts, y = BP, horiz = TRUE)) + # horiz = TRUE horizontal bars
  geom_bar(stat = "identity", fill = "hotpink") +
  labs(title = "Significant Enriched Terms AD vs Young", x = "Count", y = "Biological Process")

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          # panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}

plot_enrich_ADYoung <- edit_plots(plot_enrich_ADYoung)
plot_enrich_ADYoung



### plot unique to AD vs Old 
### when i plot in write up, make y axis scale uniform so choose genes with highest exp (base mean)

## get DESeq2 results
tx2gene = readRDS("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/DIFF-EXP/tx2gene.rds")
setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/FASTQ.files/")

files = c('./AD/22-2T-AD.fastq.gz_counts/quant.sf', #3
          './AD/25-5T-AD.fastq.gz_counts/quant.sf', #6
          './AD/28-8T-AD.fastq.gz_counts/quant.sf', #9
          './Old/13-11T-Old.fastq.gz_counts/quant.sf', #4
          './Old/14-7A-Old.fastq.gz_counts/quant.sf', #5
          './Old/15-13T-Old.fastq.gz_counts/quant.sf') #6 

txi.salmon <- tximport(files = files, 
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)

sampleTable <- data.frame(condition = factor(c("AD", "AD","AD","Old","Old","Old")))
names <-  c("AD3","AD6","AD9",
            "Old4","Old5","Old6")
rownames(sampleTable) = names
colnames(txi.salmon$counts) <- names

dds_ADvsOld <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~ condition)
dds_ADvsOld$condition <- relevel(dds_ADvsOld$condition, ref = "Old")
dds_ADvsOld <- DESeq(dds_ADvsOld)
resultsNames(dds_ADvsOld)
res_ADvsOld <- lfcShrink(dds_ADvsOld, coef="condition_AD_vs_Old", type="apeglm")

## genes to filter DESeq results
setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp_update/")
ADvsOld <- readRDS("./dataframe_filtered_ADvsOld_LFC1_basemean10_padj0.05.rds")
ADvsOld_genes <- ADvsOld$Gene

ADvsYoung <- readRDS("./dataframe_filtered_ADvsYoung_LFC1_basemean10_padj0.05.rds")
ADvsYoung_genes <- ADvsYoung$Gene

OldvsYoung <- readRDS("./dataframe_filtered_OldvsYoung_LFC1_basemean10_padj0.05.rds")
OldvsYoung_genes <- OldvsYoung$Gene

uniqueADOld <- setdiff(ADvsOld_genes, ADvsYoung_genes)
uniqueADOld <- setdiff(uniqueADOld, OldvsYoung_genes)

random_gene_ids <- sample(uniqueADOld, 10)
plot_data_list <- list()
for (gene_id in random_gene_ids) {
  plot_data <- plotCounts(dds_ADvsOld, gene=gene_id, intgroup="condition", returnData = TRUE)
  # Store the plot data in the list
  plot_data_list[[gene_id]] <- plot_data
}

for (i in names(plot_data_list)) {
  plot_data_list[[i]]$gene <- i
}

plot_data_df <- bind_rows(plot_data_list) ## combine small data frames in a list into one big dataframe

full_plot <- ggplot(plot_data_df, aes(x = condition, y=count, fill = condition)) +
  geom_boxplot() +
  scale_fill_manual(values = c("plum", "slategray1")) +
  geom_signif(comparisons = list(c("AD","Old")),
              annotations = "*", vjust = 0.7, 
              textsize = 6) +
  facet_wrap(~ gene, scales = "free_y")

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          # panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}

full_plot  <- edit_plots(full_plot )
full_plot 

### plot salmon mapping rate

sample_names <- c("AD3", "AD6", "AD9", "Old4", "Old5", "Old6", "Young1", "Young3", "Young5")
map_perc <- c(52.8, 51.8, 49.3, 54.1, 49.4, 49.7, 68.3, 65.5, 53.3)

# Create a data frame
data <- data.frame(Sample = sample_names, MappedPercentage = map_perc)

# Bar plot

mapping_rate <- ggplot(data, aes(x = Sample, y = MappedPercentage)) +
  geom_bar(stat = "identity", fill = "hotpink") +
  labs(
    x = "Sample", y = "Mapped Percentage (%)") +
  theme_minimal() +
  ylim(0, 100) +  # Set y-axis limits
  theme(axis.text.x = element_text(angle = -90, hjust = 1))  # Rotate x-axis labels

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          # panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}
mapping_rate <- edit_plots(mapping_rate)
mapping_rate

## boxplot
sample_names <- c("AD3", "AD6", "AD9", "Old4", "Old5", "Old6", "Young1", "Young3", "Young5")
map_perc <- c(52.8, 51.8, 49.3, 54.1, 49.4, 49.7, 68.3, 65.5, 53.3)
condition <- c("AD","AD","AD",
               "Old", "Old","Old",
               "Young", "Young","Young")
# Create a data frame
data <- data.frame(Sample = sample_names, MappedPercentage = map_perc, condition)

boxplot_counts <- ggplot(data, aes(x = condition, y = MappedPercentage)) +
  geom_boxplot(fill="slategrey") +
  labs(x = "Sample Group",
       y = "Percent Mapped") +
  scale_fill_discrete(name = "Sample") +
  theme_minimal() +
  geom_signif(comparisons = list(c("AD","Old"),c("Old", "Young"), c("AD", "Young")), 
              map_signif_level=TRUE)

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          # panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}
boxplot_counts <- edit_plots(boxplot_counts)
boxplot_counts



#### genes of interest #### 

### filtered DESeq2 genes ###
setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp_update/")

ADvsOld <- readRDS("./dataframe_filtered_ADvsOld_LFC1_basemean10_padj0.05.rds")
OldvsYoung <- readRDS("./dataframe_filtered_OldvsYoung_LFC1_basemean10_padj0.05.rds")
ADvsYoung <- readRDS("./dataframe_filtered_ADvsYoung_LFC1_basemean10_padj0.05.rds")

## AD vs OLD - UP- & DOWN-REGULATED
ADvsOld_up <- ADvsOld[ADvsOld$LFC > 0, ]
ADvsOld_down <- ADvsOld[ADvsOld$LFC < 0, ]
## OLD vs YOUNG - UP- & DOWN-REGULATED
OldvsYoung_up <- OldvsYoung[OldvsYoung$LFC > 0, ]
OldvsYoung_down <- OldvsYoung[OldvsYoung$LFC < 0, ]
## AD vs YOUNG - UP- & DOWN-REGULATED
ADvsYoung_up <- ADvsYoung[ADvsYoung$LFC > 0, ]
ADvsYoung_down <- ADvsYoung[ADvsYoung$LFC < 0, ]


## genes
ADOld_genes_up <- ADvsOld_up$Gene
ADOld_genes_down <- ADvsOld_down$Gene

OldYoung_genes_up <- OldvsYoung_up$Gene
OldYoung_genes_down <- OldvsYoung_down$Gene

ADYoung_genes_up <- ADvsYoung_up$Gene
ADYoung_genes_down <- ADvsYoung_down$Gene

# unique_ADOld_up <- setdiff(ADOld_genes_up, ADYoung_genes_up)
# unique_ADOld_up <- setdiff(unique_ADOld_up, OldYoung_genes_up)
# 
# unique_ADOld_down <- setdiff(ADOld_genes_down, ADYoung_genes_down)
# unique_ADOld_down <- setdiff(unique_ADOld_down, OldYoung_genes_down)
# 
# saveRDS(unique_ADOld_up, "./unique_genes_ADOld_up.rds")
# saveRDS(unique_ADOld_down, "./unique_genes_ADOld_down.rds")


#### confidence disease-related
AOY_up <- intersect(ADOld_genes_up, ADYoung_genes_up)
AA <- setdiff(AOY_up, OldYoung_genes_up)

saveRDS(AA, "./disease_genes_up.rds")

AOY_down <- intersect(ADOld_genes_down, ADYoung_genes_down)
AA_down <- setdiff(AOY_down, OldYoung_genes_down)

saveRDS(AA_down, "./disease_genes_down.rds")


##### CLUSTERPROFILER #####
### CONFIDENCE AD-REALTED ###
## GENES i put in = all filters applied
## BACKGROUND/UNIVERSE = only base mean filter

library(clusterProfiler)
setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp_update/")
#setwd("D://ABBY.windows_surface_pro/differentialexp/")

### UNIVERSES ###
## AD vs OLD
res_ADvsOld <- readRDS("./DESEQres_pairwise_ADvsOld(3each-bulk).rds")
ADvsOld_df <- data.frame(Gene = res_ADvsOld@rownames, 
                         LFC = res_ADvsOld@listData[["log2FoldChange"]],
                         pvalue = res_ADvsOld@listData[["pvalue"]],
                         padj = res_ADvsOld@listData[["padj"]],
                         basemean = res_ADvsOld@listData[["baseMean"]])
filtered_ADvsOld_basemean <- ADvsOld_df[ADvsOld_df$basemean > 10, ]


## INPUT GENES ##

### filtered DESeq2 genes ###
setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp_update/")

ADvsOld <- readRDS("./dataframe_filtered_ADvsOld_LFC1_basemean10_padj0.05.rds")
OldvsYoung <- readRDS("./dataframe_filtered_OldvsYoung_LFC1_basemean10_padj0.05.rds")
ADvsYoung <- readRDS("./dataframe_filtered_ADvsYoung_LFC1_basemean10_padj0.05.rds")

### REMOVING VERSION NUMBERS OF GENES
ADvsOld$Gene = gsub("(ENSG[0-9]+)\\.[0-9]+","\\1",ADvsOld$Gene)
filtered_ADvsOld_basemean$Gene = gsub("(ENSG[0-9]+)\\.[0-9]+","\\1",filtered_ADvsOld_basemean$Gene)
ADvsOld_genes <- ADvsOld$Gene

OldvsYoung$Gene = gsub("(ENSG[0-9]+)\\.[0-9]+","\\1",OldvsYoung$Gene)
filtered_OldvsYoung_basemean$Gene = gsub("(ENSG[0-9]+)\\.[0-9]+","\\1",filtered_OldvsYoung_basemean$Gene)
OldvsYoung_genes <- OldvsYoung$Gene

ADvsYoung$Gene = gsub("(ENSG[0-9]+)\\.[0-9]+","\\1",ADvsYoung$Gene)
filtered_ADvsYoung_basemean$Gene = gsub("(ENSG[0-9]+)\\.[0-9]+","\\1",filtered_ADvsYoung_basemean$Gene)
ADvsYoung_genes <- ADvsYoung$Gene

### to get BPs for AD Old up/down separately
## AD vs OLD - UP- & DOWN-REGULATED
ADvsOld_up <- ADvsOld[ADvsOld$LFC > 0, ]
ADvsOld_down <- ADvsOld[ADvsOld$LFC < 0, ]
## OLD vs YOUNG - UP- & DOWN-REGULATED
OldvsYoung_up <- OldvsYoung[OldvsYoung$LFC > 0, ]
OldvsYoung_down <- OldvsYoung[OldvsYoung$LFC < 0, ]
## AD vs YOUNG - UP- & DOWN-REGULATED
ADvsYoung_up <- ADvsYoung[ADvsYoung$LFC > 0, ]
ADvsYoung_down <- ADvsYoung[ADvsYoung$LFC < 0, ]

## genes
ADOld_genes_up <- ADvsOld_up$Gene
ADOld_genes_down <- ADvsOld_down$Gene

OldYoung_genes_up <- OldvsYoung_up$Gene
OldYoung_genes_down <- OldvsYoung_down$Gene

ADYoung_genes_up <- ADvsYoung_up$Gene
ADYoung_genes_down <- ADvsYoung_down$Gene

unique_ADOld_up <- intersect(ADOld_genes_up, ADYoung_genes_up)
unique_ADOld_up <- setdiff(unique_ADOld_up, OldYoung_genes_up)

unique_ADOld_down <- intersect(ADOld_genes_down, ADYoung_genes_down)
unique_ADOld_down <- setdiff(unique_ADOld_down, OldYoung_genes_down)

AD_genes <- c(unique_ADOld_up, unique_ADOld_down)
AD_df <- ADvsOld[ADvsOld$Gene %in% AD_genes, ]

gene_list <- AD_df$LFC ## ** see enrichKEGG below
names(gene_list) <- AD_df$Gene

gene_list <-na.omit(gene_list)
gene_list = sort(gene_list, decreasing = TRUE)


### AD vs Old ###
## upregulated BPs

## no background
GOenrich_ADOld <- enrichGO(unique_ADOld_up, 'org.Hs.eg.db', keyType = "ENSEMBL", ont = "BP" ,pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
## + background
GOenrich_ADOld2 <- enrichGO(unique_ADOld_up, 'org.Hs.eg.db', keyType = "ENSEMBL", ont = "BP" ,pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05,
                            readable = TRUE,universe = filtered_ADvsOld_basemean$Gene) ##filtered basemean not LFC or p value == what IS expressed in cell

saveRDS(GOenrich_ADOld, "./GOenrich_diseasegenes_up_nobackground.rds")
saveRDS(GOenrich_ADOld2, "./GOenrich_diseasegenes_up_+background.rds")

## downregulated BPs

## no background
GOenrich_ADOld3 <- enrichGO(unique_ADOld_down, 'org.Hs.eg.db', keyType = "ENSEMBL", ont = "BP" ,pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
## + background
GOenrich_ADOld4 <- enrichGO(unique_ADOld_down, 'org.Hs.eg.db', keyType = "ENSEMBL", ont = "BP" ,pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05,
                            readable = TRUE,universe = filtered_ADvsOld_basemean$Gene) ##filtered basemean not LFC or p value == what IS expressed in cell

saveRDS(GOenrich_ADOld3, "./GOenrich_diseasegenes_down_nobackground.rds")
saveRDS(GOenrich_ADOld4, "./GOenrich_diseasegenes_down_+background.rds")


### KEGG enrichment ###
# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(names(gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb='org.Hs.eg.db')

# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]

# prepare

gene_list_df <- data.frame(ENSEMBL = names(gene_list), LFC = gene_list)
gene_list_full <- inner_join(gene_list_df, ids, by = "ENSEMBL")
kegg_gene_list <- gene_list_full$LFC
names(kegg_gene_list) <- gene_list_full$ENTREZID
kegg_gene_list<-na.omit(kegg_gene_list)

# # Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
# AD_df2 = AD_df[AD_df$Gene %in% dedup_ids$ENSEMBL, ]
# 
# # Create a new column in df2 with the corresponding ENTREZ IDs
# AD_df2$entrezID = dedup_ids$ENTREZID
# 
# # Create a vector of the gene unuiverse
# kegg_gene_list <- AD_df2$LFC
# 
# # Name vector with ENTREZ ids
# names(kegg_gene_list) <- AD_df2$entrezID
# 
# # omit any NA values 
# kegg_gene_list<-na.omit(kegg_gene_list)
# # sort the list in decreasing order (required for clusterProfiler)
# kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
# names(kegg_gene_list) <- as.integer(names(kegg_gene_list))

KEGGenrich <- enrichKEGG(names(kegg_gene_list), 'hsa', keyType = "ncbi-geneid",pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05)


### universe ###
ids2 <-bitr(filtered_ADvsOld_basemean$Gene, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb='org.Hs.eg.db')
dedup_ids2 = ids2[!duplicated(ids2[c("ENSEMBL")]),]
filtered_ADvsOld_basemean <- filtered_ADvsOld_basemean[filtered_ADvsOld_basemean$Gene %in% dedup_ids2$ENSEMBL, ]
filtered_ADvsOld_basemean$entrezID = dedup_ids2$ENTREZID

KEGGenrich2 <- enrichKEGG(names(kegg_gene_list), 'hsa', keyType = "ncbi-geneid",pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05, universe = filtered_ADvsOld_basemean$entrezID)
?enrichKEGG


saveRDS(KEGGenrich, "./KEGGenrich_diseasegenes_nobackground.rds")
saveRDS(KEGGenrich2, "./KEGGenrich_diseasegenes_+background.rds")

### PATHVIEW ###
BiocManager::install("pathview")
library(pathview)

# Produce the native KEGG plot (PNG)
paths <- pathview(gene.data=kegg_gene_list, pathway.id="hsa04080", species = 'hsa', gene.idtype="entrez")
paths2 <- pathview(gene.data=kegg_gene_list, pathway.id="hsa04020", species = 'hsa', gene.idtype="entrez")
paths3 <- pathview(gene.data=kegg_gene_list, pathway.id="hsa05022", species = 'hsa', gene.idtype="entrez")
paths4 <- pathview(gene.data=kegg_gene_list, pathway.id="hsa05010", species = 'hsa', gene.idtype="entrez")
paths5 <- pathview(gene.data=kegg_gene_list, pathway.id="hsa04725", species = 'hsa', gene.idtype="entrez")
paths6 <- pathview(gene.data=kegg_gene_list, pathway.id="hsa04728", species = 'hsa', gene.idtype="entrez")
paths7 <- pathview(gene.data=kegg_gene_list, pathway.id="hsa04911", species = 'hsa', gene.idtype="entrez")

test <- readRDS("./filtered_dataframe_ADOld_down.rds")

knitr::include_graphics("dme04080.pathview.png")

### PLOTTING CLUSTERPROFILER ###
setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp_update/")
setwd("D://ABBY.windows_surface_pro/differentialexp_update/")
enrich_ADOld_up <- readRDS("./GOenrich_diseasegenes_up_+background.rds")
enrich_ADOld_down <- readRDS("./GOenrich_diseasegenes_down_+background.rds")

res_enrich_ADOld_up <- as.data.frame(enrich_ADOld_up@result)
res_enrich_ADOld_down <- as.data.frame(enrich_ADOld_down@result)
## filtering top 20 BPs
### ordering by increasing p.adjust ###
# Column to order by
column_to_order <- "Count"

# Order the data frame based on the specified column
ordered_res_enrich_ADOld_up <- res_enrich_ADOld_up[order(res_enrich_ADOld_up[[column_to_order]], decreasing = TRUE), ]
ordered_res_enrich_ADOld_up <- ordered_res_enrich_ADOld_up[1:20, ]
df_enrich_ADOld_up <- data.frame(BP = ordered_res_enrich_ADOld_up$Description, padj = ordered_res_enrich_ADOld_up$p.adjust, counts = ordered_res_enrich_ADOld_up$Count)

ordered_res_enrich_ADOld_down <- res_enrich_ADOld_down[order(res_enrich_ADOld_down[[column_to_order]], decreasing = TRUE), ]
ordered_res_enrich_ADOld_down <- ordered_res_enrich_ADOld_down[1:20, ]
df_enrich_ADOld_down <- data.frame(BP = ordered_res_enrich_ADOld_down$Description, padj = ordered_res_enrich_ADOld_down$p.adjust, counts = ordered_res_enrich_ADOld_down$Count)

## up
plot_enrich_ADOld_up <- ggplot(df_enrich_ADOld_up, aes(x = counts, y = BP), horiz = TRUE) + # horiz = TRUE horizontal bars
  geom_bar(stat = "identity", fill = "grey50", colour = "springgreen") +  
  labs(title = "AD-specific Significant Enriched Terms (UP)", x = "Count", y = "Biological Process") +
  scale_x_continuous(expand = c(0,0)) 

                       
#Darisia Edit_plot
edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, t = 0.35, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, r = 0.35, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
} 

plot_enrich_ADOld_up<- edit_plots(plot_enrich_ADOld_up)
plot_enrich_ADOld_up

saveRDS(plot_enrich_ADOld_up, "D://ABBY.windows_surface_pro/diffexp_results/clusterprof_disease_up.rds")
                     
                     # plot_enrich_ADOld_up <- ggplot(df_enrich_ADOld_up, aes(x = counts, y = BP), horiz = TRUE) + # horiz = TRUE horizontal bars
                     #   geom_bar(stat = "identity", fill = "s", colour="orchid") +  
                     #   labs(title = "AD-specific Significant Enriched Terms (UP)", x = "Count", y = "Biological Process") +
                     #   scale_x_continuous(expand = c(0,0))                     
saveRDS(plot_enrich_ADOld_up, "./clusterprof_UP_09-09.rds")

## down
plot_enrich_ADOld_down <- ggplot(df_enrich_ADOld_down, aes(x = counts, y = BP), horiz = TRUE) + # horiz = TRUE horizontal bars
  geom_bar(stat = "identity", fill = "grey50", colour = "red") +  
  labs(title = "AD-specific Significant Enriched Terms (DOWN)", x = "Count", y = "Biological Process") +
  scale_x_continuous(expand = c(0,0))

plot_enrich_ADOld_down <- edit_plots(plot_enrich_ADOld_down)
plot_enrich_ADOld_down

saveRDS(plot_enrich_ADOld_down, "D://ABBY.windows_surface_pro/diffexp_results/clusterprof_disease_down.rds")
saveRDS(plot_enrich_ADOld_down, "./clusterprof_DOWN_09-09.rds")

##### disease genes #####

setwd("D://ABBY.windows_surface_pro/differentialexp_update/")

AD_up <- readRDS("./disease_genes_up.rds")
AD_down <- readRDS("./disease_genes_down.rds")


ADOld <- readRDS("./dataframe_filtered_ADvsOld_LFC1_basemean10_padj0.05.rds")
ADYoung <- readRDS("./dataframe_filtered_ADvsYoung_LFC1_basemean10_padj0.05.rds")

filtered_ADOld <-  ADOld[(ADOld$Gene %in% AD_up), ]
filtered_ADYoung <- ADYoung[(ADYoung$Gene %in% AD_up), ]

saveRDS(filtered_ADOld, "./filtered_dataframe_ADOld_up.rds")
saveRDS(filtered_ADYoung, "./filtered_dataframe_ADYoung_up.rds")


filtered_ADOld_down <-  ADOld[(ADOld$Gene %in% AD_down), ]
filtered_ADYoung_down <- ADYoung[(ADYoung$Gene %in% AD_down), ]

saveRDS(filtered_ADOld_down, "./filtered_dataframe_ADOld_down.rds")
saveRDS(filtered_ADYoung_down, "./filtered_dataframe_ADYoung_down.rds")


down_ADO <- filtered_ADOld_down$Gene
downADY <- filtered_ADYoung_down$Gene

downADY %in% down_ADO



#### boxplot for ages ####
samples <- c("AD3", "AD6", "AD9", "Old4", "Old5", "Old6", "Young1", "Young3", "Young5")
ages <- c(71, 61, 79, 68, 70, 77, 42, 57, 59)
condition <- c("AD","AD","AD",
               "Old", "Old","Old",
               "Young", "Young","Young")
donorinfo <- data.frame(Sample = samples, Age = ages, Condition = condition)

df <- data.frame(ages)
rownames(df) <- condition
df


# Convert the result to a data frame
boxplot_counts <- ggplot(donorinfo, aes(x = condition, y = ages)) +
  geom_boxplot(fill="slategrey") +
  labs(x = "Sample Group",
       y = "Age") +
  scale_fill_discrete(name = "Sample") +
  theme_minimal() +
  geom_signif(comparisons = list(c("AD","Old"),c("Old", "Young"), c("AD", "Young")), 
              map_signif_level=TRUE)

#Darisia Edit_plot
edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          # panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
} 
boxplot_counts <- edit_plots(boxplot_counts)
boxplot_counts

## adding individual samples onto plot
boxplot_counts <- boxplot_counts +
  geom_point(data = donorinfo, aes(x = condition, y = ages, colour = condition), ## change shape: shape = condition
             size = 3) +
  geom_text_repel(label= donorinfo$Sample) +
  scale_colour_manual(values = c("AD" = "slategray3", "Old" = "plum", "Young" = "orange"))

boxplot_counts

#### do own stats test and add annotation manually

### plot 10 random genes from DESeq2 results
## genes to filter DESeq results

## Yi's tx2gene

tx2gene = readRDS("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/DIFF-EXP/tx2gene.rds")
setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/FASTQ.files/")

###### AD vs Old ######

files = c('./AD/22-2T-AD.fastq.gz_counts/quant.sf', #3
          './AD/25-5T-AD.fastq.gz_counts/quant.sf', #6
          './AD/28-8T-AD.fastq.gz_counts/quant.sf', #9
          './Old/13-11T-Old.fastq.gz_counts/quant.sf', #4
          './Old/14-7A-Old.fastq.gz_counts/quant.sf', #5
          './Old/15-13T-Old.fastq.gz_counts/quant.sf') #6 

txi.salmon <- tximport(files = files, 
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)

sampleTable <- data.frame(condition = factor(c("AD", "AD","AD","Old","Old","Old")))
names <-  c("AD3","AD6","AD9",
            "Old4","Old5","Old6")
rownames(sampleTable) = names
colnames(txi.salmon$counts) <- names

dds_ADvsOld <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~ condition)
dds_ADvsOld$condition <- relevel(dds_ADvsOld$condition, ref = "Old")
dds_ADvsOld <- DESeq(dds_ADvsOld)
resultsNames(dds_ADvsOld)
res_ADvsOld <- lfcShrink(dds_ADvsOld, coef="condition_AD_vs_Old", type="apeglm")

setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp_update/")
ADvsOld <- readRDS("./dataframe_filtered_ADvsOld_LFC1_basemean10_padj0.05.rds")
ADvsOld_genes <- ADvsOld$Gene

ADvsYoung <- readRDS("./dataframe_filtered_ADvsYoung_LFC1_basemean10_padj0.05.rds")
ADvsYoung_genes <- ADvsYoung$Gene

OldvsYoung <- readRDS("./dataframe_filtered_OldvsYoung_LFC1_basemean10_padj0.05.rds")
OldvsYoung_genes <- OldvsYoung$Gene

uniqueAD <- intersect(ADvsOld_genes, ADvsYoung_genes)
uniqueAD <- setdiff(uniqueAD, OldvsYoung_genes)

### plotting DESeq2 10 random genes
random_gene_ids <- sample(uniqueAD, 10)
plot_data_list <- list()
for (gene_id in random_gene_ids) {
  plot_data <- plotCounts(dds_ADvsOld, gene=gene_id, intgroup="condition", returnData = TRUE)
  # Store the plot data in the list
  plot_data_list[[gene_id]] <- plot_data
}

for (i in names(plot_data_list)) {
  plot_data_list[[i]]$gene <- i
}

plot_data_df <- bind_rows(plot_data_list) ## combine small data frames in a list into one big dataframe

full_plot <- ggplot(plot_data_df, aes(x = condition, y=count, fill = condition)) +
  geom_boxplot() +
  scale_fill_manual(values = c("plum", "slategray1")) +
  geom_signif(comparisons = list(c("AD","Old")),
              annotations = "*", vjust = 0.7, 
              textsize = 6) +
  facet_wrap(~ gene, scales = "free_y")

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          # panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}

full_plot  <- edit_plots(full_plot )
full_plot


### fixing multiqc 
setwd("D://ABBY.windows_surface_pro/differentialexp_update/")
seq_quality <- read_csv("./fastqc_per_sequence_quality_scores_plot.csv")
seq_quality <- as.data.frame(seq_quality)
seq_quality_long <- seq_quality %>%
  pivot_longer(cols = -`Mean Sequence Quality (Phred Score)`, 
               names_to = "Sample", 
               values_to = "Count")

seq_quality_long <- mutate(seq_quality_long, Condition = case_when(grepl("AD", Sample) ~ "AD", grepl("Old", Sample) ~ "Old", grepl("Young", Sample) ~ "Young"))

seq_quality_plot <- ggplot(seq_quality_long, aes(x = `Mean Sequence Quality (Phred Score)`, y = Count, colour = Condition)) +
  geom_line(aes(group=Sample), linewidth = 0.8) +
  labs(x = "Mean Sequence Quality (Phred Score)", y = "Count", color = "") +
  scale_colour_manual(values = c("AD" = "orangered", "Old" = "darkslateblue", "Young" = "gold")) +
  scale_x_continuous(limits = c(0, 36)) +
  theme(legend.position = "top")

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}

seq_quality_plot <- edit_plots(seq_quality_plot) 
seq_quality_plot 

# ggplot(data=seq_quality_long, aes(x=`Mean Sequence Quality (Phred Score)`, y=Count, fill=Condition)) +
#   #geom_line(size=2 * 1.5, aes(group= Sample), colour = "black")+
#   geom_line(size=2, aes(colour= Condition))+
#   scale_fill_brewer(palette = "Set1")

saveRDS(seq_quality_plot, "D://ABBY.windows_surface_pro/diffexp_results/seqqual_plot_all_11sept.rds")

## per seq qual score for selected samples 

samples_keep <- c("22-2T-AD", "25-5T-AD", "28-8T-AD", "13-11T-Old", "14-7A-Old", "15-13T-Old", "2-12A-Young", "4-13A-Young", "6-14A-Young)
")

# Filter the dataframe to include only rows with a specific value in a specific column
selected_seq_quality_long <- seq_quality_long %>%
  filter(Sample %in% samples_keep)

selected_seq_quality_plot <- ggplot(selected_seq_quality_long, aes(x = `Mean Sequence Quality (Phred Score)`, y = Count, color = Condition)) +
  geom_line(aes(group=Sample), linewidth = 0.8) +
  labs(x = "Mean Sequence Quality (Phred Score)", y = "Count", color = "") +
  scale_colour_manual(values = c("AD" = "orangered", "Old" = "darkslateblue", "Young" = "gold")) +
  scale_x_continuous(limits = c(0, 36)) +
  theme(legend.position = "top")

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}

selected_seq_quality_plot <- edit_plots(selected_seq_quality_plot)
selected_seq_quality_plot

saveRDS(selected_seq_quality_plot, "D://ABBY.windows_surface_pro/diffexp_results/seqqual_plot_3each_11sept.rds")

### per base qual score

pb_seq_quality <- read_csv("./fastqc_per_base_sequence_quality_plot.csv")
pb_seq_quality <- as.data.frame(pb_seq_quality)
pb_seq_quality_long <- pb_seq_quality %>%
  pivot_longer(cols = -`Position (bp)`, 
               names_to = "Sample", 
               values_to = "Phred Score")

pb_seq_quality_long <- mutate(pb_seq_quality_long, Condition = case_when(grepl("AD", Sample) ~ "AD", grepl("Old", Sample) ~ "Old", grepl("Young", Sample) ~ "Young"))

pb_seq_quality_plot <- ggplot(pb_seq_quality_long, aes(x = `Position (bp)`, y = `Phred Score`, color = Condition)) +
  geom_line(aes(group=Sample), linewidth = 0.8) +
  labs(x = "Position (bp)", y = "Phred Score", color = "") +
  scale_colour_manual(values = c("AD" = "orangered", "Old" = "darkslateblue", "Young" = "gold")) +
  scale_y_continuous(limits = c(0, 40)) +
  theme(legend.position = "top")

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}

pb_seq_quality_plot <- edit_plots(pb_seq_quality_plot)
pb_seq_quality_plot 

saveRDS(pb_seq_quality_plot, "D://ABBY.windows_surface_pro/diffexp_results/pb_seqqual_plot_all_11sept.rds")

## per base qual score for selected samples 

samples_keep <- c("22-2T-AD", "25-5T-AD", "28-8T-AD", "13-11T-Old", "14-7A-Old", "15-13T-Old", "2-12A-Young", "4-13A-Young", "6-14A-Young)
")

# Filter the dataframe to include only rows with a specific value in a specific column
pb_selected_seq_quality_long <- pb_seq_quality_long %>%
  filter(Sample %in% samples_keep)

pb_selected_seq_quality_plot <- ggplot(pb_selected_seq_quality_long, aes(x = `Position (bp)`, y = `Phred Score`, color = Condition)) +
  geom_line(aes(group=Sample), linewidth = 0.8) +
  labs(x = "Position (bp)", y = "Phred Score", color = "") +
  scale_colour_manual(values = c("AD" = "orangered", "Old" = "darkslateblue", "Young" = "gold")) +
  scale_y_continuous(limits = c(0, 40)) +
  theme(legend.position = "top")

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}

pb_selected_seq_quality_plot <- edit_plots(pb_selected_seq_quality_plot)
pb_selected_seq_quality_plot

saveRDS(pb_selected_seq_quality_plot, "D://ABBY.windows_surface_pro/diffexp_results/pb_seqqual_plot_3each_11sept.rds")

### per base N content

pb_ncont <- read_csv("./fastqc_per_base_n_content_plot.csv")
pb_ncont <- as.data.frame(pb_ncont)
pb_ncont_long <- pb_ncont %>%
  pivot_longer(cols = -`Position in Read (bp)`, 
               names_to = "Sample", 
               values_to = "Percentage N-Count")

pb_ncont_long <- mutate(pb_ncont_long, Condition = case_when(grepl("AD", Sample) ~ "AD", grepl("Old", Sample) ~ "Old", grepl("Young", Sample) ~ "Young"))

pb_ncont_plot <- ggplot(pb_ncont_long, aes(x = `Position in Read (bp)`, y = `Percentage N-Count`, color = Condition)) +
  geom_line(aes(group=Sample), linewidth = 0.8) +
  labs(x = "Position (bp)", y = "Percentage N-Count", color = "") +
  scale_colour_manual(values = c("AD" = "orangered", "Old" = "darkslateblue", "Young" = "gold")) +
  scale_y_continuous(limits = c(0, 6)) +
  theme(legend.position = "top")

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}

pb_ncont_plot <- edit_plots(pb_ncont_plot)
pb_ncont_plot

saveRDS(pb_ncont_plot, "D://ABBY.windows_surface_pro/diffexp_results/pb_Ncontent_plot_all_11sept.rds")

## per base N content for selected samples 

samples_keep <- c("22-2T-AD", "25-5T-AD", "28-8T-AD", "13-11T-Old", "14-7A-Old", "15-13T-Old", "2-12A-Young", "4-13A-Young", "6-14A-Young)
")

# Filter the dataframe to include only rows with a specific value in a specific column
selected_pb_ncont_long <- pb_ncont_long %>%
  filter(Sample %in% samples_keep)

selected_pb_ncont_long_plot <- ggplot(selected_pb_ncont_long, aes(x = `Position in Read (bp)`, y = `Percentage N-Count`, color = Condition)) +
  geom_line(aes(group=Sample), linewidth = 0.8) +
  labs(x = "Position (bp)", y = "Percentage N-Count", color = "") +
  scale_colour_manual(values = c("AD" = "orangered", "Old" = "darkslateblue", "Young" = "gold")) +
  scale_y_continuous(limits = c(0, 6)) +
  theme(legend.position = "top")

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}

selected_pb_ncont_long_plot <- edit_plots(selected_pb_ncont_long_plot)
selected_pb_ncont_long_plot

saveRDS(selected_pb_ncont_long_plot, "D://ABBY.windows_surface_pro/diffexp_results/pb_Ncontent_plot_3each_9sept.rds")


### join 3 QC plots ###

setwd("D://ABBY.windows_surface_pro/differentialexp_update/")

plot1_all <- readRDS("./seqqual_plot_all_9sept.rds")
plot2_all <- readRDS("./pb_seqqual_plot_all_9sept.rds")
plot3_all <- readRDS("./pb_Ncontent_plot_all_9sept.rds")

multi_panel_figure_all <- ggarrange(plot1_all, plot2_all, plot3_all, ncol = 3,
                                    common.legend = TRUE,   # Add a single legend for all plots
                                    legend = "top", # Adjust the legend position as needed
                                    labels = "AUTO",
                                    label.x = 0.15
)


multi_panel_figure_all +
  theme(legend.text = element_text(size = 20),
        legend.key.height  = unit(2, "lines"),
        legend.key.width = unit(2, "lines"))

?ggarrange

# Display the multi-panel figure
multi_panel_figure_all

### selected samples ###
setwd("D://ABBY.windows_surface_pro/differentialexp_update/")

plot1_3each <- readRDS("./seqqual_plot_3each_9sept.rds")
plot2_3each <- readRDS("./pb_seqqual_plot_3each_9sept.rds")
plot3_3each <- readRDS("./pb_Ncontent_plot_3each_9sept.rds")

plot1_3each +
  theme(legend.text = element_text(size=12))

plot2_3each +
  theme(legend.text = element_text(size=12))

plot3_3each +
  theme(legend.text = element_text(size=12))

multi_panel_figure_3each <- ggarrange(plot1_3each, plot2_3each, plot3_3each, ncol = 3,
                                    common.legend = TRUE,   # Add a single legend for all plots
                                    legend = "top", # Adjust the legend position as needed
                                    labels = "AUTO",
                                    label.x = 0.15
)

# Display the multi-panel figure
multi_panel_figure_3each



###### Looking at DESeq2 results ######

setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp_update/")
setwd("D://ABBY.windows_surface_pro/differentialexp_update/")

## AD vs Old ## 
res_ADvsOld <- readRDS("./DESEQres_pairwise_ADvsOld(3each-bulk).rds")
ADvsOld_df <- data.frame(Gene = res_ADvsOld@rownames, 
                         LFC = res_ADvsOld@listData[["log2FoldChange"]],
                         pvalue = res_ADvsOld@listData[["pvalue"]],
                         padj = res_ADvsOld@listData[["padj"]],
                         basemean = res_ADvsOld@listData[["baseMean"]])

filtered_ADvsOld_basemean <- ADvsOld_df[ADvsOld_df$basemean > 10, ]
filtered_ADvsOld_basemean <- na.omit(filtered_ADvsOld_basemean)
filtered_ADvsOld_basemean <- mutate(filtered_ADvsOld_basemean, Direction = case_when(padj < 0.05 & LFC > 1 ~ "UP", padj < 0.05 & LFC < -1 ~ "DOWN", padj >= 0.05 | abs(LFC) <= 1 ~ "notDE"))
saveRDS(filtered_ADvsOld_basemean, "./ADvsOld_dataframe_basemeancutoff_direction.rds")
## volcano plot no LFC cutoff
volcano <- ggplot(filtered_ADvsOld_basemean, aes(x = LFC, y = -log10(padj), colour = Direction)) +
  geom_point(alpha=0.5, size=3) +
  scale_colour_manual(values = c("UP" = "mediumpurple", "DOWN" = "mediumpurple", "notDE" = "grey82")) + # Color for non-significant and significant genes
  #scale_colour_manual(values = c("UP" = "orange", "DOWN" = "springgreen", "notDE" = "grey72")) +
  labs(
    x = "Log Fold Change (LFC)",
    y = "-log10(Adjusted p-value)",
    colour = " ") +
  theme(legend.position = "none")


edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
} 
volcano <- edit_plots(volcano) 
volcano


num_upregulated <- sum(filtered_ADvsOld_basemean$Direction == "UP")
num_downregulated <- sum(filtered_ADvsOld_basemean$Direction == "DOWN")

# Add annotations to the plot
volcano +
  geom_text(aes(label = paste(num_upregulated),
                x = 7,
                y = 10),
            hjust = 0.5, vjust = 0.5, colour = "black", size=6) +
  geom_text(aes(label = paste(num_downregulated),
                x = -6,
                y=10),
            hjust = 0.5, vjust = 0.5, colour = "black", size=6)

saveRDS(volcano, "D://ABBY.windows_surface_pro/diffexp_results/volcano_ADOld.rds")

?geom_label

# ?geom_text 
# geom_label(aes(x = max(filtered_ADvsOld_basemean$LFC),
#                  y = max(-log10(filtered_ADvsOld_basemean$padj)),
#                  label = paste("Up: ", num_upregulated, "\nDown: ", num_downregulated)),
#              hjust = 1, vjust = 1)
# geom_label(aes(label=650), nudge_x = 5, nudge_y = 5)
# ?geom_label
# "\nDown: ", num_downregulated
# max(filtered_ADvsOld_basemean$LFC)
# max(-log10(filtered_ADvsOld_basemean$padj))
# (aes(label = 650,
#      x = 7,
#      y = 10),
#   hjust = 0.5, vjust = 0.5, colour = "black", size=6) +
#   geom_label(aes(label = paste(num_downregulated),
#                  x = -6,
#                  y=10),
#              hjust = 0.5, vjust = 0.5, colour = "black", size=6)

## Old vs Young
res_OldvsYoung <- readRDS("./DESEQres_pairwise_OldvsYoung(3each-bulk).rds")
OldvsYoung_df <- data.frame(Gene = res_OldvsYoung@rownames, 
                            LFC = res_OldvsYoung@listData[["log2FoldChange"]],
                            pvalue = res_OldvsYoung@listData[["pvalue"]],
                            padj = res_OldvsYoung@listData[["padj"]],
                            basemean = res_OldvsYoung@listData[["baseMean"]])
filtered_OldvsYoung_basemean <- OldvsYoung_df[OldvsYoung_df$basemean > 10, ]
filtered_OldvsYoung_basemean <- na.omit(filtered_OldvsYoung_basemean)
filtered_OldvsYoung_basemean <- mutate(filtered_OldvsYoung_basemean, Direction = case_when(padj < 0.05 & LFC > 1 ~ "UP", padj < 0.05 & LFC < -1 ~ "DOWN", padj >= 0.05 | abs(LFC) <= 1 ~ "notDE"))
saveRDS(filtered_OldvsYoung_basemean, "./OldvsYoung_dataframe_basemeancutoff_direction.rds")
## volcano plot no LFC cutoff
volcano <- ggplot(filtered_OldvsYoung_basemean, aes(x = LFC, y = -log10(padj), color = Direction)) +
  geom_point(alpha=0.5, size=3) +
  scale_colour_manual(values = c("UP" = "lightgreen", "DOWN" = "lightgreen", "notDE" = "grey82")) + # Color for non-significant and significant genes
  #scale_colour_manual(values = c("UP" = "orange", "DOWN" = "springgreen", "notDE" = "grey72")) +
  labs(
    x = "Log Fold Change (LFC)",
    y = "-log10(Adjusted p-value)",
    colour = ""
  ) +
  theme(legend.position = "none")

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
} 
volcano <- edit_plots(volcano) 
volcano

num_upregulated <- sum(filtered_OldvsYoung_basemean$Direction == "UP")
num_downregulated <- sum(filtered_OldvsYoung_basemean$Direction == "DOWN")

# Add annotations to the plot
volcano +
  geom_text(aes(label = paste(num_upregulated),
                x = 7,
                y = 10),
            hjust = 0.5, vjust = 0.5, colour = "black", size=6) +
  geom_text(aes(label = paste(num_downregulated),
                x = -6,
                y=10),
            hjust = 0.5, vjust = 0.5, colour = "black", size=6)

saveRDS(volcano, "D://ABBY.windows_surface_pro/diffexp_results/volcano_OldYoung.rds")

# volcano <- ggplot(filtered_OldvsYoung_basemean, aes(x = LFC, y = -log10(padj), color = Direction, fill = Direction)) +
#   geom_point(shape = 21, alpha=0.5, size=3) +
#   scale_fill_manual(values = c("UP" = "pink", "DOWN" = "lightskyblue", "notDE" = "grey82"), labels = c("Down-regulated", "No Change", "Up-regulated")) + # Color for non-significant and significant genes
#   #scale_colour_manual(values = c("UP" = "orange", "DOWN" = "springgreen", "notDE" = "grey72")) +
#   labs(
#     x = "Log Fold Change (LFC)",
#     y = "-log10(Adjusted p-value)",
#     colour = "", fill = " "
#   ) +
#   theme(legend.position = "none")

## AD vs Young
res_ADvsYoung <- readRDS("./DESEQres_pairwise_ADvsYoung(3each-bulk).rds")
ADvsYoung_df <- data.frame(Gene = res_ADvsYoung@rownames, 
                           LFC = res_ADvsYoung@listData[["log2FoldChange"]],
                           pvalue = res_ADvsYoung@listData[["pvalue"]],
                           padj = res_ADvsYoung@listData[["padj"]],
                           basemean = res_ADvsYoung@listData[["baseMean"]])
filtered_ADvsYoung_basemean <- ADvsYoung_df[ADvsYoung_df$basemean > 10, ]
filtered_ADvsYoung_basemean <- na.omit(filtered_ADvsYoung_basemean)
filtered_ADvsYoung_basemean <- mutate(filtered_ADvsYoung_basemean, Direction = case_when(padj < 0.05 & LFC > 1 ~ "UP", padj < 0.05 & LFC < -1 ~ "DOWN", padj >= 0.05 | abs(LFC) <= 1 ~ "notDE"))
saveRDS(filtered_ADvsYoung_basemean, "./ADvsYoung_dataframe_basemeancutoff_direction.rds")

## volcano plot no LFC cutoff
volcano <- ggplot(filtered_ADvsYoung_basemean, aes(x = LFC, y = -log10(padj), color = Direction)) +
  geom_point(alpha=0.5, size=3) +
  scale_colour_manual(values = c("UP" = "orange", "DOWN" = "orange", "notDE" = "grey82")) + # Color for non-significant and significant genes
  #scale_colour_manual(values = c("UP" = "orange", "DOWN" = "springgreen", "notDE" = "grey72")) +
  labs(
    x = "Log Fold Change (LFC)",
    y = "-log10(Adjusted p-value)",
    colour = ""
  ) +
  theme(legend.position = "none")

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
} 
volcano <- edit_plots(volcano) 
volcano

num_upregulated <- sum(filtered_ADvsYoung_basemean$Direction == "UP")
num_downregulated <- sum(filtered_ADvsYoung_basemean$Direction == "DOWN")

# Add annotations to the plot
volcano +
  geom_text(aes(label = paste(num_upregulated),
                x = 7,
                y = 12),
            hjust = 0.5, vjust = 0.5, colour = "black", size=6) +
  geom_text(aes(label = paste(num_downregulated),
                x = -6,
                y=12),
            hjust = 0.5, vjust = 0.5, colour = "black", size=6)

saveRDS(volcano, "D://ABBY.windows_surface_pro/diffexp_results/volcano_ADYoung.rds")


custom_legend <- ggplot() +
  geom_point(aes(x = 1:3, y = 1), fill = c("pink", "lightskyblue", "grey82"), colour=c("orange", "springgreen", "grey72"), shape = 21, size = 4) +
  annotate("text", x = 1.05:3.05, y = 1, label = c("Up-regulated", "Down-regulated", "No Change"), vjust = 0.5, hjust = 0, size = 4, color = "black") +
  theme_void() +
  theme(legend.position = "none")  # Remove default legend

custom_legend


# custom_legend <- ggplot() +
#   geom_point(aes(x = 1, y = 1), shape = 21, size = 4, color = "orange", fill = "pink") +
#   geom_point(aes(x = 2, y = 2), shape = 21, size = 4, color = "springgreen", fill = "lightskyblue") +
#   geom_point(aes(x = 3, y = 3), shape = 21, size = 4, color = "grey72", fill = "grey82") +
#   annotate("text", x = 1, y = 1, label = "Up-regulated", vjust = 0.5, hjust = 0, size = 4, color = "black") +
#   annotate("text", x = 2, y = 2, label = "Down-regulated", vjust = 0.5, hjust = 0, size = 4, color = "black") +
#   annotate("text", x = 3, y = 3, label = "No Change", vjust = 0.5, hjust = 0, size = 4, color = "black") +
#   theme_void()


setwd("D://ABBY.windows_surface_pro/differentialexp_update/")

plot1 <- readRDS("./volcano_ADOld_pinkblue.rds")
plot2 <- readRDS("./volcano_OldYoung_pinkblue.rds")
plot3 <- readRDS("./volcano_ADYoung_pinkblue.rds")


multi <- ggarrange(plot1, plot2, plot3, ncol = 3,
                                                    # Add a single legend for all plots
                                                    # Adjust the legend position as needed
                                    labels = "AUTO",
                                    label.x = 0.0,
                                    align = "hv",
                                    vjust = 1.8)


multi

ggarrange(custom_legend, multi, ncol = 1)


eights = c(1, 5)



### plot salmon mapping rate

sample_names <- c("AD3", "AD6", "AD9", "Old4", "Old5", "Old6", "Young1", "Young3", "Young5")
map_perc <- c(52.8, 51.8, 49.3, 54.1, 49.4, 49.7, 68.3, 65.5, 53.3)
condition <- c("AD", "AD", "AD", "Old", "Old", "Old", "Young", "Young", "Young")

# Create a data frame
data <- data.frame(Sample = sample_names, MappedPercentage = map_perc, Condition = condition)

# Bar plot
mapping_rate <- ggplot(data, aes(x = Sample, y = MappedPercentage, fill = Condition)) +
  geom_bar(stat = "identity") +
  labs(
    x = "Sample", y = "Mapped Percentage (%)") +
  scale_fill_manual(values = c("AD" = "slategray3", "Old" = "plum", "Young" = "orange")) +
  ylim(0, 100) +  # Set y-axis limits
  theme(axis.text.x = element_text(angle = -90, hjust = 1))  # Rotate x-axis labels

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}
mapping_rate <- edit_plots(mapping_rate)
mapping_rate
saveRDS(mapping_rate, "./mapping_plot.rds")
## boxplot
sample_names <- c("AD3", "AD6", "AD9", "Old4", "Old5", "Old6", "Young1", "Young3", "Young5")
map_perc <- c(52.8, 51.8, 49.3, 54.1, 49.4, 49.7, 68.3, 65.5, 53.3)
condition <- c("AD","AD","AD",
               "Old", "Old","Old",
               "Young", "Young","Young")
# Create a data frame
data <- data.frame(Sample = sample_names, MappedPercentage = map_perc, condition)

boxplot_counts <- ggplot(data, aes(x = condition, y = MappedPercentage)) +
  geom_boxplot(fill="grey80") +
  labs(x = "Sample Group",
       y = "Percent Mapped",
       colour = "") +
  theme_minimal() +
  geom_signif(comparisons = list(c("AD","Old"),c("Old", "Young"), c("AD", "Young")), 
              map_signif_level=TRUE)

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
}
boxplot_counts <- edit_plots(boxplot_counts)
boxplot_counts

## adding individual samples onto plot
boxplot_counts <- boxplot_counts +
  geom_point(aes(x = Condition, y = MappedPercentage, colour = Condition), ## change shape: shape = condition
             size = 3) +
  geom_text_repel(label= data$Sample) +
  scale_colour_manual(values = c("AD" = "slategray3", "Old" = "plum", "Young" = "orange")) +
  scale_y_continuous(limits = c(0, 100)) +
  theme(legend.position = "top")
  

boxplot_counts

saveRDS(boxplot_counts, "./mappingrate_boxplot.rds")


###### PLOT TOTAL GENE COUNTS PER SAMPLE (selected samples) ######

#### boxlot ####

tx2gene = readRDS("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/DIFF-EXP/tx2gene.rds")
setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/FASTQ.files/")

files = c('./AD/22-2T-AD.fastq.gz_counts/quant.sf', #3
          './AD/25-5T-AD.fastq.gz_counts/quant.sf', #6
          './AD/28-8T-AD.fastq.gz_counts/quant.sf', #9
          './Old/13-11T-Old.fastq.gz_counts/quant.sf', #4
          './Old/14-7A-Old.fastq.gz_counts/quant.sf', #5
          './Old/15-13T-Old.fastq.gz_counts/quant.sf', #6
          './Young/2-12A-Young.fastq.gz_counts/quant.sf', #1
          './Young/4-13A-Young.fastq.gz_counts/quant.sf', #3
          './Young/6-14A-Young.fastq.gz_counts/quant.sf') #5

txi.salmon <- tximport(files = files, 
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)

sampleTable <- data.frame(condition = factor(c("AD", "AD","AD","Old","Old","Old","Young","Young","Young")))
names <-  c("AD3","AD6","AD9",
            "Old4","Old5","Old6",
            "Young1","Young3","Young5")
rownames(sampleTable) = names
colnames(txi.salmon$counts) <- names
colnames(txi.salmon$abundance) <- names

counts_df <- as.data.frame(txi.salmon$counts)

# Calculate total gene counts per sample
totals <- colSums(counts_df)
samples <- names(totals)
condition <- c("AD","AD","AD",
               "Old", "Old","Old",
               "Young", "Young","Young")

# Convert the result to a data frame
total_counts_df <- data.frame(Sample = samples, TotalCount = totals, condition)

boxplot_counts <- ggplot(total_counts_df, aes(x = condition, y = TotalCount)) +
  geom_boxplot(fill="grey80") +
  labs(x = "Sample Group",
       y = "Total Gene Count",
       colour = "") +
  geom_signif(comparisons = list(c("AD","Old"),c("Old", "Young"), c("AD", "Young")), 
              map_signif_level=TRUE)

#Darisia Edit_plot
edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
} 
boxplot_counts <- edit_plots(boxplot_counts)
boxplot_counts

## adding individual samples onto plot
boxplot_counts <- boxplot_counts +
  geom_point(data = total_counts_df, aes(x = condition, y = TotalCount, colour = condition), ## change shape: shape = condition
             size = 3) +
  geom_text_repel(label= total_counts_df$Sample) +
  scale_colour_manual(values = c("AD" = "slategray3", "Old" = "plum", "Young" = "orange"))

boxplot_counts


setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp_update/")
saveRDS(boxplot_counts, "./totalcounts_boxplot_3each.rds")


### 6 September 2023
###### PLOT TOTAL GENE COUNTS PER SAMPLE ######
## all samples ## 

tx2gene = readRDS("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/DIFF-EXP/tx2gene.rds")
setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/FASTQ.files/")

tx2gene = readRDS("D://ABBY.windows_surface_pro/DIFF-EXP/tx2gene.rds")
setwd("D://ABBY.windows_surface_pro/FASTQ.files/")

files = c('./AD/20-1T-AD.fastq.gz_counts/quant.sf',
          './AD/21-1A-AD.fastq.gz_counts/quant.sf',
          './AD/22-2T-AD.fastq.gz_counts/quant.sf',
          './AD/23-2A-AD.fastq.gz_counts/quant.sf',
          './AD/24-3T-AD.fastq.gz_counts/quant.sf',
          './AD/25-5T-AD.fastq.gz_counts/quant.sf',
          './AD/26-3A-AD.fastq.gz_counts/quant.sf',
          './AD/27-5A-AD.fastq.gz_counts/quant.sf',
          './AD/28-8T-AD.fastq.gz_counts/quant.sf',
          './AD/29-6T-AD.fastq.gz_counts/quant.sf',
          './AD/30-9T-AD.fastq.gz_counts/quant.sf',
          './AD/31-7T-AD.fastq.gz_counts/quant.sf',
          './Old/10-8A-Old.fastq.gz_counts/quant.sf',
          './Old/11-10T-Old.fastq.gz_counts/quant.sf',
          './Old/12-6A-Old.fastq.gz_counts/quant.sf', 
          './Old/13-11T-Old.fastq.gz_counts/quant.sf', 
          './Old/14-7A-Old.fastq.gz_counts/quant.sf', 
          './Old/15-13T-Old.fastq.gz_counts/quant.sf', 
          './Old/16-14T-Old.fastq.gz_counts/quant.sf', 
          './Old/17-9A-Old.fastq.gz_counts/quant.sf',
          './Old/18-10A-Old.fastq.gz_counts/quant.sf',
          './Old/19-11A-Old.fastq.gz_counts/quant.sf',
          './Young/2-12A-Young.fastq.gz_counts/quant.sf',
          './Young/3-17T-Young.fastq.gz_counts/quant.sf',
          './Young/4-13A-Young.fastq.gz_counts/quant.sf',
          './Young/5-18T-Young.fastq.gz_counts/quant.sf',
          './Young/6-14A-Young.fastq.gz_counts/quant.sf',
          './Young/7-19T-Young.fastq.gz_counts/quant.sf',
          './Young/8-15A-Young.fastq.gz_counts/quant.sf',
          './Young/9-16A-Young.fastq.gz_counts/quant.sf')

txi.salmon <- tximport(files = files, 
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)

sampleTable <- data.frame(condition = factor(c("AD", "AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","Old","Old","Old","Old","Old","Old","Old","Old","Old","Old", "Young","Young","Young","Young","Young","Young","Young","Young")))
names <-  c("AD1", "AD2","AD3","AD4","AD5","AD6","AD7","AD8","AD9","AD10","AD11","AD12",
            "Old1","Old2","Old3","Old4","Old5","Old6","Old7","Old8","Old9","Old10",
            "Young1","Young2","Young3","Young4","Young5","Young6","Young7","Young8")
rownames(sampleTable) = names
colnames(txi.salmon$counts) <- names

counts_df <- as.data.frame(txi.salmon$counts)

# Calculate total gene counts per sample
totals <- colSums(counts_df)
samples <- names(totals)
condition <- c("AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD",
               "Old", "Old","Old","Old","Old","Old","Old","Old","Old","Old",
               "Young", "Young","Young","Young","Young","Young","Young","Young")

# Convert the result to a data frame
total_counts_df <- data.frame(Sample = samples, TotalCount = totals, condition)
total_counts_df <- total_counts_df %>%
  mutate(condition = factor(condition, levels = c("AD", "Old", "Young")),
         SampleNumeric = as.numeric(gsub("[^0-9]", "", Sample)),
         SampleLabel = paste0(condition, SampleNumeric))

# Reorder the levels of SampleLabel within each condition
total_counts_df <- total_counts_df %>%
  group_by(condition) %>%
  mutate(SampleLabel = factor(SampleLabel, levels = unique(SampleLabel))) %>%
  ungroup()

# Create a bar plot using ggplot2
plot_count <- ggplot(total_counts_df, aes(x = SampleLabel, y = TotalCount, fill = condition, colour = condition)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample", y = "Cumulative Gene Count", fill = "", colour = "") +
  scale_fill_manual(values = c("AD" = "orangered", "Old" = "darkslateblue", "Young" = "gold")) +
  scale_colour_manual(values = c("AD" = "orangered", "Old" = "darkslateblue", "Young" = "gold")) +
  theme(legend.position = "top") +
  scale_y_continuous(expand = c(0,0)) +
  #theme(axis.title.y = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = -90, hjust = 1)) # Rotate x-axis labels 

#Darisia Edit_plot
edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, t = 0.35, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, r = 0.35, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
} 

plot_count <- edit_plots(plot_count)
plot_count

#setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp_update/")
setwd("D://ABBY.windows_surface_pro/diffexp_results/")
ggsave("./cumulativegenecount_bar_all_11sept.png", plot_count, width = 15, height = 10)
saveRDS(plot_count, "./cumulativegenecount_bar_all_11sept.rds")

### box plot total gene counts all samples ###

boxplot_counts <- ggplot(total_counts_df, aes(x = condition, y = TotalCount)) +
  geom_boxplot(fill="white", alpha=0.15) +
  labs(x = "Sample Group",
       y = "Cumulative Gene Count",
       colour = "") +
  geom_signif(comparisons = list(c("AD","Old"),c("Old", "Young"), c("AD", "Young")),
              annotations = "ns", vjust = -0.5, 
              textsize = 4, colour = "black", y_position = 27900000, fontface="bold")
boxplot_counts <- edit_plots(boxplot_counts)
boxplot_counts

boxplot_counts <- boxplot_counts +
  #geom_point(shape=1) +
  geom_point(aes(x = condition, y = TotalCount, colour = condition), ## change shape: shape = condition
            size = 3) +
  geom_text_repel(label= total_counts_df$SampleLabel, size = 3.5, fontface="bold") +
  scale_colour_manual(values = c("AD" = "orangered", "Old" = "darkslateblue", "Young" = "gold")) +
  theme(legend.position = "top")


boxplot_counts

# boxplot_counts <- boxplot_counts +
#   #geom_point(shape=1) +
#   geom_point(shape=21, aes(x = condition, y = TotalCount, fill = condition), ## change shape: shape = condition
#              colour = "black", size = 3) +
#   geom_text_repel(label= total_counts_df$SampleLabel, size = 3.5, fontface="bold") +
#   scale_fill_manual(values = c("AD" = "black", "Old" = "azure4", "Young" = "wheat3")) +
#   theme(legend.position = "top")

# setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp_update/")
# ggsave("./cumulativegenecounts_all_boxplot_b&w.png", boxplot_counts, width = 15, height = 10)

setwd("D://ABBY.windows_surface_pro/diffexp_results/")
ggsave("./cumulativegenecount_box_all_11sept.png", boxplot_counts, width = 15, height = 10)
saveRDS(boxplot_counts, "./cumulativegenecount_box_all.rds")

###### PLOT TOTAL GENE COUNTS PER SAMPLE (selected samples) ######
#### boxlot ####

# tx2gene = readRDS("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/DIFF-EXP/tx2gene.rds")
# setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/FASTQ.files/")

tx2gene = readRDS("D://ABBY.windows_surface_pro/DIFF-EXP/tx2gene.rds")
setwd("D://ABBY.windows_surface_pro/FASTQ.files/")

files = c('./AD/22-2T-AD.fastq.gz_counts/quant.sf', #3
          './AD/25-5T-AD.fastq.gz_counts/quant.sf', #6
          './AD/28-8T-AD.fastq.gz_counts/quant.sf', #9
          './Old/13-11T-Old.fastq.gz_counts/quant.sf', #4
          './Old/14-7A-Old.fastq.gz_counts/quant.sf', #5
          './Old/15-13T-Old.fastq.gz_counts/quant.sf', #6
          './Young/2-12A-Young.fastq.gz_counts/quant.sf', #1
          './Young/4-13A-Young.fastq.gz_counts/quant.sf', #3
          './Young/6-14A-Young.fastq.gz_counts/quant.sf') #5

txi.salmon <- tximport(files = files, 
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)

sampleTable <- data.frame(condition = factor(c("AD", "AD","AD","Old","Old","Old","Young","Young","Young")))
names <-  c("AD3","AD6","AD9",
            "Old4","Old5","Old6",
            "Young1","Young3","Young5")
rownames(sampleTable) = names
colnames(txi.salmon$counts) <- names
counts_df <- as.data.frame(txi.salmon$counts)

# Calculate total gene counts per sample
totals <- colSums(counts_df)
samples <- names(totals)
condition <- c("AD","AD","AD",
               "Old", "Old","Old",
               "Young", "Young","Young")

# Convert the result to a data frame
total_counts_df <- data.frame(Sample = samples, TotalCount = totals, condition)
total_counts_df <- total_counts_df %>%
  mutate(condition = factor(condition, levels = c("AD", "Old", "Young")),
         SampleNumeric = as.numeric(gsub("[^0-9]", "", Sample)),
         SampleLabel = paste0(condition, SampleNumeric))
# Reorder the levels of SampleLabel within each condition
total_counts_df <- total_counts_df %>%
  group_by(condition) %>%
  mutate(SampleLabel = factor(SampleLabel, levels = unique(SampleLabel))) %>%
  ungroup()

#Darisia Edit_plot
edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
} 

boxplot_counts <- ggplot(total_counts_df, aes(x = condition, y = TotalCount)) +
  geom_boxplot(fill="white", alpha=0.15, colour = "black") +
  labs(x = "Sample Group",
       y = "Cumulative Gene Count",
       fill = "") +
  geom_signif(comparisons = list(c("AD","Old"),c("Old", "Young"), c("AD", "Young")),
              annotations = "ns", vjust = -0.5, 
              textsize = 4, colour = "black", y_position = 27900000, fontface="bold")

boxplot_counts <- edit_plots(boxplot_counts)
boxplot_counts

boxplot_counts <- boxplot_counts +
  #geom_point(shape=1) +
  geom_point(aes(x = condition, y = TotalCount, colour = condition), ## change shape: shape = condition
            size = 3) +
  geom_text_repel(label= total_counts_df$SampleLabel, size = 3.5, fontface="bold") +
  scale_colour_manual(values = c("AD" = "orangered", "Old" = "darkslateblue", "Young" = "gold")) +
  theme(legend.position = "top")


boxplot_counts

# boxplot_counts <- ggplot(total_counts_df, aes(x = condition, y = TotalCount)) +
#   geom_boxplot(fill="grey80") +
#   labs(x = "Sample Group",
#        y = "Cumulative Gene Count",
#        colour = "") +
#   scale_fill_discrete(name = "Sample") +
#   theme_minimal() +
#   geom_signif(comparisons = list(c("AD","Old"),c("Old", "Young"), c("AD", "Young")),
#               annotations = "ns", vjust = -0.5, 
#               textsize = 4, colour = "dimgrey", y_position = 28000000)
# boxplot_counts <- edit_plots(boxplot_counts)
# boxplot_counts
# 
# boxplot_counts <- boxplot_counts +
#   geom_point(aes(x = condition, y = TotalCount, colour = condition), ## change shape: shape = condition
#              size = 2) +
#   geom_text_repel(label= total_counts_df$SampleLabel, size = 3.5) +
#   scale_colour_manual(values = c("AD" = "slategray3", "Old" = "plum", "Young" = "orange")) +
#   theme(legend.position = "top")
# 
# boxplot_counts

setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp_update/")
ggsave("./cumulativegenecounts_3each_boxplot_b&w.png", boxplot_counts, width = 15, height = 10)

setwd("D://ABBY.windows_surface_pro/diffexp_results/")
ggsave("./cumulativegenecount_box_3each_11sept.png", boxplot_counts, width = 15, height = 10)
saveRDS(boxplot_counts, "./cumulativegenecount_box_3each_11.rds")

### plot salmon mapping rate for all samples

sample_names <- c("AD1", "AD2","AD3","AD4","AD5","AD6","AD7","AD8","AD9","AD10","AD11","AD12",
                  "Old1","Old2","Old3","Old4","Old5","Old6","Old7","Old8","Old9","Old10",
                  "Young1","Young2","Young3","Young4","Young5","Young6","Young7","Young8")

map_perc <- c(45.8, 44.9, 52.8, 52.4, 49.4, 51.8, 53, 48, 49.3, 59.1, 49.1, 44.4,
              55.7, 50.1, 54.8, 54.1, 49.4, 49.7, 53.3, 57.3, 56.8, 43.5,
              68.3, 44.4, 65.5, 56.3, 53.3, 52.5, 59.3, 46.2)
condition <- c("AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD",
               "Old", "Old","Old","Old","Old","Old","Old","Old","Old","Old",
               "Young", "Young","Young","Young","Young","Young","Young","Young")

# Create a data frame
data <- data.frame(Sample = sample_names, MappedPercentage = map_perc, condition)
data <- data %>%
  mutate(condition = factor(condition, levels = c("AD", "Old", "Young")),
         SampleNumeric = as.numeric(gsub("[^0-9]", "", Sample)),
         SampleLabel = paste0(condition, SampleNumeric))

# Reorder the levels of SampleLabel within each condition
data <- data %>%
  group_by(condition) %>%
  mutate(SampleLabel = factor(SampleLabel, levels = unique(SampleLabel))) %>%
  ungroup()

# Bar plot
plot_count <- ggplot(data, aes(x = SampleLabel, y = MappedPercentage, fill = condition)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample", y = "Mapped Percentage (%)", fill = "", colour = "") +
  scale_fill_manual(values = c("AD" = "orangered", "Old" = "darkslateblue", "Young" = "gold")) +
#  scale_colour_manual(values = c("AD" = "black", "Old" = "black", "Young" = "black")) +
  theme(legend.position = "top") +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = -90, hjust = 1)) # Rotate x-axis labels 

#Darisia Edit_plot
edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, t = 0.35, unit = 'cm'),
                                      vjust = -5,
                                      size = 12),
          axis.title.y = element_text(margin = margin(l = 1, r = 0.35, unit = 'cm'),
                                      vjust = 5,
                                      size = 12),
          axis.text.x = element_text(size = 10,
                                     colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          legend.text = element_text(size = 10,
                                     colour = 'black'),
          legend.title = element_text(size = 12,
                                      colour = 'black'),
          strip.text = element_text(
            size = 12,
            colour = "black",
            face = "bold"
          ),
          strip.background = element_blank(),
          plot.margin = unit(c(1,1,1,1), "cm"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black")
    ) # Not sure what else I need to add - but for now this if fine
  return(ggplot_object)
} 
plot_count <- edit_plots(plot_count)
plot_count

setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp_update/")
ggsave("./mapping%_all_barplot_b&w.png", plot_count, width = 15, height = 10)

setwd("D://ABBY.windows_surface_pro/diffexp_results/")
ggsave("./maprate_bar_all_11sept.png", plot_count, width = 15, height = 10)
saveRDS(plot_count,"./maprate_bar_all_11sept.rds")
### box plot total gene counts all samples ###

boxplot_counts <- ggplot(data, aes(x = condition, y = MappedPercentage)) +
  geom_boxplot(fill="white", alpha=0.15, colour = "black") +
  labs(x = "Sample Group",
       y = "Mapped Percentage (%)",
       fill = "") +
  geom_signif(comparisons = list(c("AD","Old"),c("Old", "Young"), c("AD", "Young")),
              annotations = "ns", vjust = -0.5, 
              textsize = 4, colour = "black", y_position = 85, fontface="bold")+
  ylim(0, 100)   # Set y-axis limits
boxplot_counts <- edit_plots(boxplot_counts)
boxplot_counts

boxplot_counts <- boxplot_counts +
  #geom_point(shape=1) +
  geom_point(aes(x = condition, y = MappedPercentage, colour = condition), ## change shape: shape = condition
            size = 3) +
  geom_text_repel(label= data$SampleLabel, size = 3.5, fontface="bold") +
  scale_colour_manual(values = c("AD" = "orangered", "Old" = "darkslateblue", "Young" = "gold")) +
  theme(legend.position = "top")


boxplot_counts

# mapping_rate <- ggplot(data, aes(x = SampleLabel, y = MappedPercentage, fill = condition)) +
#   geom_bar(stat = "identity") +
#   labs(
#     x = "Sample", y = "Mapped Percentage (%)", fill = "") +
#   scale_fill_manual(values = c("AD" = "slategray3", "Old" = "plum", "Young" = "orange")) +
#   ylim(0, 100) +  # Set y-axis limits
#   theme(axis.text.x = element_text(angle = -90, hjust = 1)) + # Rotate x-axis labels
#   theme(legend.position = "top") 
# 
# edit_plots <- function(
#     ggplot_object
# ) {
#   ggplot_object <- ggplot_object +
#     theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
#                                       vjust = -5,
#                                       size = 12),
#           axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
#                                       vjust = 5,
#                                       size = 12),
#           axis.text.x = element_text(size = 10,
#                                      colour = 'black'),
#           axis.text.y = element_text(size = 10,
#                                      colour = 'black'),
#           legend.text = element_text(size = 10,
#                                      colour = 'black'),
#           legend.title = element_text(size = 12,
#                                       colour = 'black'),
#           strip.text = element_text(
#             size = 12,
#             colour = "black",
#             face = "bold"
#           ),
#           strip.background = element_blank(),
#           plot.margin = unit(c(1,1,1,1), "cm"),
#           panel.background = element_blank(),
#           panel.grid = element_blank(),
#           axis.line = element_line(colour = "black")
#     ) # Not sure what else I need to add - but for now this if fine
#   return(ggplot_object)
# }
# mapping_rate <- edit_plots(mapping_rate)
# mapping_rate
# 
# saveRDS(mapping_rate, "./mapping_plot.rds")
# 
# ### boxplot mapping rate all samples ### 
# boxplot_counts <- ggplot(data, aes(x = condition, y = MappedPercentage)) +
#   geom_boxplot(fill="grey80") +
#   labs(x = "Sample Group",
#        y = "Mapped Percentage (%)",
#        colour = "") +
#   ylim(0, 100) +  # Set y-axis limits
#   geom_signif(comparisons = list(c("AD","Old"),c("Old", "Young"), c("AD", "Young")),
#               annotations = "ns", vjust = -0.5, 
#               textsize = 4, colour = "dimgrey", y_position = 85)
# boxplot_counts <- edit_plots(boxplot_counts)
# boxplot_counts
# 
# boxplot_counts <- boxplot_counts +
#   geom_point(aes(x = condition, y = MappedPercentage, colour = condition), ## change shape: shape = condition
#              size = 2) +
#   geom_text_repel(label= data$SampleLabel, size = 3.5) +
#   scale_colour_manual(values = c("AD" = "slategray3", "Old" = "plum", "Young" = "orange")) +
#   theme(legend.position = "top")
# 
# 
# boxplot_counts

setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp_update/")
ggsave("./mappingrate_all_boxplot_b&w.png", boxplot_counts, width = 15, height = 10)

setwd("D://ABBY.windows_surface_pro/diffexp_results/")
ggsave("./maprate_box_all_11sept.png", boxplot_counts, width = 15, height = 10)
saveRDS(boxplot_counts,"./maprate_box_all_11sept.rds")

# Bar plot
# mapping_rate <- ggplot(data, aes(x = Sample, y = MappedPercentage, fill = Condition)) +
#   geom_bar(stat = "identity") +
#   labs(
#     x = "Sample", y = "Mapped Percentage (%)") +
#   scale_fill_manual(values = c("AD" = "slategray3", "Old" = "plum", "Young" = "orange")) +
#   ylim(0, 100) +  # Set y-axis limits
#   theme(axis.text.x = element_text(angle = -90, hjust = 1))  # Rotate x-axis labels


### plot salmon mapping rate selected samples ###

sample_names <- c("AD3", "AD6", "AD9", "Old4", "Old5", "Old6", "Young1", "Young3", "Young5")
map_perc <- c(52.8, 51.8, 49.3, 54.1, 49.4, 49.7, 68.3, 65.5, 53.3)
condition<- c("AD", "AD", "AD", "Old", "Old", "Old", "Young", "Young", "Young")

# Create a data frame
data <- data.frame(Sample = sample_names, MappedPercentage = map_perc, condition)
data <- data %>%
  mutate(condition = factor(condition, levels = c("AD", "Old", "Young")),
         SampleNumeric = as.numeric(gsub("[^0-9]", "", Sample)),
         SampleLabel = paste0(condition, SampleNumeric))

# Reorder the levels of SampleLabel within each condition
data <- data %>%
  group_by(condition) %>%
  mutate(SampleLabel = factor(SampleLabel, levels = unique(SampleLabel))) %>%
  ungroup()

### boxplot mapping rate ### 

boxplot_counts <- ggplot(data, aes(x = condition, y = MappedPercentage)) +
  geom_boxplot(fill="white", alpha=0.15, colour = "black") +
  labs(x = "Sample Group",
       y = "Mapped Percentage (%)",
       fill = "") +
  geom_signif(comparisons = list(c("AD","Old"),c("Old", "Young"), c("AD", "Young")),
              annotations = "ns", vjust = -0.5, 
              textsize = 4, colour = "black", y_position = 85, fontface="bold")+
  ylim(0, 100)   # Set y-axis limits
boxplot_counts <- edit_plots(boxplot_counts)
boxplot_counts

boxplot_counts <- boxplot_counts +
  #geom_point(shape=1) +
  geom_point(aes(x = condition, y = MappedPercentage, colour = condition), ## change shape: shape = condition
            size = 3) +
  geom_text_repel(label= data$SampleLabel, size = 3.5, fontface="bold") +
  scale_colour_manual(values = c("AD" = "orangered", "Old" = "darkslateblue", "Young" = "gold")) +
  theme(legend.position = "top")


boxplot_counts

# boxplot_counts <- ggplot(data, aes(x = condition, y = MappedPercentage)) +
#   geom_boxplot(fill="grey80") +
#   labs(x = "Sample Group",
#        y = "Mapped Percentage (%)",
#        colour = "") +
#   ylim(0, 100) +  # Set y-axis limits
#   geom_signif(comparisons = list(c("AD","Old"),c("Old", "Young"), c("AD", "Young")),
#               annotations = "ns", vjust = -0.5, 
#               textsize = 4, colour = "dimgrey", y_position = 85)
# boxplot_counts <- edit_plots(boxplot_counts)
# boxplot_counts
# 
# boxplot_counts <- boxplot_counts +
#   geom_point(aes(x = condition, y = MappedPercentage, colour = condition), ## change shape: shape = condition
#              size = 2) +
#   geom_text_repel(label= data$SampleLabel, size = 3.5) +
#   scale_colour_manual(values = c("AD" = "slategray3", "Old" = "plum", "Young" = "orange")) +
#   theme(legend.position = "top")
# 
# 
# boxplot_counts

setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp_update/")
ggsave("./mappingrate_3each_boxplot_b&w.png", boxplot_counts, width = 15, height = 10)

setwd("D://ABBY.windows_surface_pro/diffexp_results/")
ggsave("./maprate_box_3each_11sept.png", boxplot_counts, width = 15, height = 10)
saveRDS(boxplot_counts,"./maprate_box_3each_11sept.rds")
#### boxplot for ages ####

samples <- c("AD3", "AD6", "AD9", "Old4", "Old5", "Old6", "Young1", "Young3", "Young5")
ages <- c(71, 61, 79, 68, 70, 77, 42, 57, 59)
condition <- c("AD","AD","AD",
               "Old", "Old","Old",
               "Young", "Young","Young")
donorinfo <- data.frame(Sample = samples, Age = ages, Condition = condition)


# # Create a scatter plot
# scatter_plot <- ggplot(donorinfo, aes(x = Condition, y = Age)) +
#   geom_point(aes(colour = Condition)) +
#   labs(x = "Condition", y = "Age", colour = "") +
#   geom_text_repel(label= donorinfo$Sample, size = 3.5) +
#   scale_colour_manual(values = c("AD" = "black", "Old" = "plum", "Young" = "orange")) +
#   theme(legend.position = "top") +
#   ylim(30, 80)  # Set y-axis limits
# 
# edit_plots <- function(
#     ggplot_object
#   ) {
#     ggplot_object <- ggplot_object +
#       theme(axis.title.x = element_text(margin = margin(b = 1, unit = 'cm'),
#                                         vjust = -5,
#                                         size = 12),
#             axis.title.y = element_text(margin = margin(l = 1, unit = 'cm'),
#                                         vjust = 5,
#                                         size = 12),
#             axis.text.x = element_text(size = 10,
#                                        colour = 'black'),
#             axis.text.y = element_text(size = 10,
#                                        colour = 'black'),
#             legend.text = element_text(size = 10,
#                                        colour = 'black'),
#             legend.title = element_text(size = 12,
#                                         colour = 'black'),
#             strip.text = element_text(
#               size = 12,
#               colour = "black",
#               face = "bold"
#             ),
#             strip.background = element_blank(),
#             plot.margin = unit(c(1,1,1,1), "cm"),
#             panel.background = element_blank(),
#             panel.grid = element_blank(),
#             axis.line = element_line(colour = "black")
#       ) # Not sure what else I need to add - but for now this if fine
#     return(ggplot_object)
#   }
# scatter_plot <- edit_plots(scatter_plot)
# scatter_plot 

boxplot_ages <- ggplot(donorinfo, aes(x = Condition, y = Age)) +
  geom_boxplot(fill="white", alpha=0.15, colour = "black") +
  labs(x = "Sample Group",
       y = "Age",
       colour = "") +
  ylim(30, 80)   # Set y-axis limits
boxplot_ages <- edit_plots(boxplot_ages)
boxplot_ages

## adding individual samples onto plot
# boxplot_ages <- boxplot_ages +
#   geom_point(data = donorinfo, aes(x = condition, y = ages, colour = condition), ## change shape: shape = condition
#              size = 3) +
#   geom_text_repel(label= donorinfo$Sample) +
#   scale_colour_manual(values = c("AD" = "slategray3", "Old" = "plum", "Young" = "orange"))
# 
# boxplot_ages <- edit_plots(boxplot_ages)
# boxplot_ages

boxplot_ages <- boxplot_ages +
  #geom_point(shape=1) +
  geom_point(aes(x = condition, y = ages, colour = condition), ## change shape: shape = condition
            size = 3) +
  geom_text_repel(label= donorinfo$Sample, size = 3.5, fontface="bold") +
  scale_colour_manual(values = c("AD" = "orangered", "Old" = "darkslateblue", "Young" = "gold")) +
  theme(legend.position = "top")


boxplot_ages


# boxplot_ages <- ggplot(donorinfo, aes(x = Condition, y = Age)) +
#   geom_boxplot(fill="grey80") +
#   labs(x = "Sample Group",
#        y = "Age", colour = "") +
#   ylim(30, 80)  # Set y-axis limits
# 
# ## adding individual samples onto plot
# boxplot_ages <- boxplot_ages +
#   geom_point(data = donorinfo, aes(x = condition, y = ages, colour = condition), ## change shape: shape = condition
#              size = 3) +
#   geom_text_repel(label= donorinfo$Sample) +
#   scale_colour_manual(values = c("AD" = "slategray3", "Old" = "plum", "Young" = "orange"))
# 
# boxplot_ages <- edit_plots(boxplot_ages)
# boxplot_ages

setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp_update/")
ggsave("./ages_3each_boxplot_b&w.png", boxplot_ages, width = 15, height = 10)

setwd("D://ABBY.windows_surface_pro/diffexp_results/")
ggsave("./ages_box_3each_11sept.png", boxplot_ages, width = 15, height = 10)
saveRDS(boxplot_ages, "./ages_box_3each.rds")

setwd("D://ABBY.windows_surface_pro/differentialexp_update/")

filtered_ADOld_down <- readRDS("./filtered_dataframe_ADOld_down.rds")
filtered_ADYoung_down <- readRDS("./filtered_dataframe_ADYoung_down.rds")
filtered_ADOld_up <- readRDS("./filtered_dataframe_ADOld_up.rds")
filtered_ADYoung_up <- readRDS ("./filtered_dataframe_ADYoung_up.rds")

filtered_ADOld_down$Gene = gsub("(ENSG[0-9]+)\\.[0-9]+","\\1",filtered_ADOld_down$Gene)
filtered_ADYoung_down$Gene = gsub("(ENSG[0-9]+)\\.[0-9]+","\\1",filtered_ADYoung_down$Gene)
filtered_ADOld_up$Gene = gsub("(ENSG[0-9]+)\\.[0-9]+","\\1",filtered_ADOld_up$Gene)
filtered_ADYoung_up$Gene = gsub("(ENSG[0-9]+)\\.[0-9]+","\\1",filtered_ADYoung_up$Gene)

filtered_ADOld_down$symbol <- mapIds(org.Hs.eg.db, keys = filtered_ADOld_down$Gene, keytype = "ENSEMBL", column = "SYMBOL")
filtered_ADYoung_down$symbol <- mapIds(org.Hs.eg.db, keys = filtered_ADYoung_down$Gene, keytype = "ENSEMBL", column = "SYMBOL")
filtered_ADOld_up$symbol <- mapIds(org.Hs.eg.db, keys = filtered_ADOld_up$Gene, keytype = "ENSEMBL", column = "SYMBOL")
filtered_ADYoung_up$symbol <- mapIds(org.Hs.eg.db, keys = filtered_ADYoung_up$Gene, keytype = "ENSEMBL", column = "SYMBOL")

saveRDS(filtered_ADOld_down, "./disease_genes_down_ADOld_+symbol.rds")
saveRDS(filtered_ADOld_up, "./disease_genes_up_ADOld_+symbol.rds")
saveRDS(filtered_ADYoung_down, "./disease_genes_down_ADYoung_+symbol.rds")
saveRDS(filtered_ADYoung_up, "./disease_genes_up_ADYoung_+symbol.rds")

# Column to order by
column_to_order <- "LFC"
# Order the data frame based on the specified column
ordered_filtered_ADOld_down <- filtered_ADOld_down[order(filtered_ADOld_down[[column_to_order]]), ]
ordered_filtered_ADYoung_down <- filtered_ADYoung_down[order(filtered_ADYoung_down[[column_to_order]]), ]

AO_down <- ordered_filtered_ADOld_down[1:50, ]
AY_down <- ordered_filtered_ADYoung_down[1:50, ]

geneslook <- AO_down$Gene %in% AY_down$Gene

AO_down[geneslook, ]



ordered_filtered_ADOld_up <- filtered_ADOld_up[order(filtered_ADOld_up[[column_to_order]], decreasing = TRUE), ]
ordered_filtered_ADYoung_up <- filtered_ADYoung_up[order(filtered_ADYoung_up[[column_to_order]], decreasing = TRUE), ]

AO_up <- ordered_filtered_ADOld_up[1:50, ]
AY_up <- ordered_filtered_ADYoung_up[1:50, ]

geneslook <- AO_up$Gene %in% AY_up$Gene

AO_up[geneslook, ]



setwd("D://ABBY.windows_surface_pro/differentialexp_update/")
disease_AO_down <- readRDS("./disease_genes_down_ADOld_+symbol.rds")
disease_AO_up <- readRDS("./disease_genes_up_ADOld_+symbol.rds")
disease_AY_down <- readRDS("./disease_genes_down_ADYoung_+symbol.rds")
disease_AY_up <- readRDS("./disease_genes_up_ADYoung_+symbol.rds")

enrich_ADOld <- readRDS("./GOenrich_ADOld_all_+background.rds")
enrich_OldYoung <- readRDS("./GOenrich_OldYoung_+background.rds")
enrich_ADYoung <- readRDS("./GOenrich_ADYoung_+background.rds")


disease_up <- readRDS("./GOenrich_diseasegenes_up_+background.rds")
disease_down <- readRDS("./GOenrich_diseasegenes_down_+background.rds")


res_enrich_OldYoung <- as.data.frame(enrich_OldYoung@result)
res_disease_up <- as.data.frame(disease_up@result)
res_disease_down <- as.data.frame(disease_down@result)

# Column to order by
column_to_order <- "Count"
# Order the data frame based on the specified column
ordered_res_enrich_OldYoung <- res_enrich_OldYoung[order(res_enrich_OldYoung[[column_to_order]], decreasing = TRUE), ]
ordered_disease_up <- res_disease_up[order(res_disease_up[[column_to_order]], decreasing = TRUE), ]
ordered_disease_down <- res_disease_down[order(res_disease_down[[column_to_order]], decreasing = TRUE), ]




saveRDS(ordered_disease_down, "./ordered_diseasegenes_down.rds")
saveRDS(ordered_disease_up, "./ordered_diseasegenes_up.rds")

































?ggarrange
# c("AD" = "darkolivegreen3", "Old" = "mediumvioletred", "Young" = "midnightblue")
# files = c(
#           './AD/22-2T-AD.fastq.gz_counts/quant.sf', 3
#           './AD/25-5T-AD.fastq.gz_counts/quant.sf', 6
#           './AD/26-3A-AD.fastq.gz_counts/quant.sf', 7
#           './AD/27-5A-AD.fastq.gz_counts/quant.sf', 8
#           './AD/28-8T-AD.fastq.gz_counts/quant.sf', 9
#           './Old/10-8A-Old.fastq.gz_counts/quant.sf', 1
#           './Old/11-10T-Old.fastq.gz_counts/quant.sf',2
#           './Old/12-6A-Old.fastq.gz_counts/quant.sf', 3
#           './Old/13-11T-Old.fastq.gz_counts/quant.sf', 4
#           './Old/14-7A-Old.fastq.gz_counts/quant.sf', 5
#           './Old/15-13T-Old.fastq.gz_counts/quant.sf', 6
#           './Old/16-14T-Old.fastq.gz_counts/quant.sf', 7
#           './Young/2-12A-Young.fastq.gz_counts/quant.sf', 1
#           './Young/3-17T-Young.fastq.gz_counts/quant.sf', 2
#           './Young/4-13A-Young.fastq.gz_counts/quant.sf', 3
#           './Young/6-14A-Young.fastq.gz_counts/quant.sf', 5






# total_counts_df$SampleLabel <- fct_reorder(total_counts_df$SampleLabel, total_counts_df$SampleNumeric)

# total_counts_df <- total_counts_df %>%
#   group_by(condition) %>%
#   mutate(row_number = as.numeric(gsub("[^0-9]", "", Sample))) %>%
#   ungroup() %>%
#   mutate(SampleLabel = paste0(condition, row_number))

