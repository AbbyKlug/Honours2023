library(ggplot2)
library(DESeq2)
library(tximport)
library(tidyverse)
library(ggpubr)
library(ggsignif)
library(ggrepel)

setwd("D://ABBY.windows_surface_pro/diffexp_results/")

## MultiQC ALL

plot1_all <- readRDS("./seqqual_biglegend_all.rds")
plot2_all <- readRDS("./pbseqqual_biglegend_all.rds")
plot3_all <- readRDS("./ncont_biglegend_all.rds")

seqqual <- plot1_all +
theme(legend.key.size = unit(3, "lines"),
                               legend.text = element_text(size = 15))
# saved pdf 10 x 10
saveRDS(seqqual, "./seqqual_biglegend_all.rds")

pbseqqual <- plot2_all +
  theme(legend.key.size = unit(3, "lines"),
        legend.text = element_text(size = 15))
pbseqqual

# saved pdf 10 x 10
saveRDS(pbseqqual, "./pbseqqual_biglegend_all.rds")

ncont <- plot3_all +
  theme(legend.key.size = unit(3, "lines"),
        legend.text = element_text(size = 15))
ncont

# saved pdf 10 x 10
saveRDS(ncont, "./ncont_biglegend_all.rds")


## MultiQC plots

multi_panel_figure_all <- ggarrange(plot1_all, plot2_all, plot3_all, ncol = 3,
                                    common.legend = TRUE,   # Add a single legend for all plots
                                    legend = "top", # Adjust the legend position as needed
                                    labels = "AUTO",
                                    label.x = 0.17,
                                    align = "hv",
                                    vjust = 1.8)
multi_panel_figure_all

## MultiQC 3 each

plot1_all <- readRDS("./seqqual_plot_3each_11sept.rds")
plot2_all <- readRDS("./pb_seqqual_plot_3each_11sept.rds")
plot3_all <- readRDS("./pb_Ncontent_plot_3each_9sept.rds")

seqqual <- plot1_all +
  theme(legend.key.size = unit(3, "lines"),
        legend.text = element_text(size = 15))
seqqual
# saved pdf 10 x 10
saveRDS(seqqual, "./seqqual_biglegend_3each.rds")

pbseqqual <- plot2_all +
  theme(legend.key.size = unit(3, "lines"),
        legend.text = element_text(size = 15))
pbseqqual

# saved pdf 10 x 10
saveRDS(pbseqqual, "./pbseqqual_biglegend_3each.rds")

ncont <- plot3_all +
  theme(legend.key.size = unit(3, "lines"),
        legend.text = element_text(size = 15))
ncont

# saved pdf 10 x 10
saveRDS(ncont, "./ncont_biglegend_3each.rds")


## MultiQC plots

plot1 <- readRDS("./seqqual_biglegend_3each.rds")
plot2 <- readRDS("./pbseqqual_biglegend_3each.rds")
plot3 <- readRDS("./ncont_biglegend_3each.rds")

multi_panel_figure_3each <- ggarrange(plot1, plot2, plot3, ncol = 3,
                                    common.legend = TRUE,   # Add a single legend for all plots
                                    legend = "top", # Adjust the legend position as needed
                                    labels = "AUTO",
                                    label.x = 0.17,
                                    align = "hv",
                                    vjust = 1.8)
multi_panel_figure_3each


## mapping rate

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

boxplot_counts <- ggplot(data, aes(x = condition, y = MappedPercentage)) +
  geom_boxplot(fill="white", alpha=0.15, colour = "black") +
  labs(x = "Sample Group",
       y = "Mapped Percentage (%)",
       fill = "", colour=" ") +
  geom_signif(comparisons = list(c("AD","Old"),c("Old", "Young"), c("AD", "Young")),
              annotations = "ns", vjust = -0.5, 
              textsize = 6, colour = "black", y_position = 85, fontface="bold")+
  ylim(0, 100)   # Set y-axis limits

boxplot_counts <- boxplot_counts +
  geom_point(aes(x = condition, y = MappedPercentage, colour = condition), ## change shape: shape = condition
             size = 3) +
  geom_text_repel(label= data$SampleLabel, size = 3.5, fontface="bold") +
  scale_colour_manual(values = c("AD" = "orangered", "Old" = "darkslateblue", "Young" = "gold")) +
  theme(legend.position = "top")


boxplot_counts

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, t = 0.35, unit = 'cm'),
                                      vjust = -5,
                                      size = 14, face="bold"),
          axis.title.y = element_text(margin = margin(l = 1, r = 0.35, unit = 'cm'),
                                      vjust = 5,
                                      size = 14, face="bold"),
          axis.text.x = element_text(size = 13,
                                     colour = 'black'),
          axis.text.y = element_text(size = 13,
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
boxplot_counts <- boxplot_counts + theme(legend.key.size = unit(3, "lines"),
                       legend.text = element_text(size = 15))

boxplot_counts
saveRDS(boxplot_counts, "./maprate_box_all_new.rds")


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
              textsize = 6, colour = "black", y_position = 85, fontface="bold")+
  ylim(0, 100)   # Set y-axis limits


boxplot_counts <- boxplot_counts +
  geom_point(aes(x = condition, y = MappedPercentage, colour = condition), ## change shape: shape = condition
             size = 3) +
  geom_text_repel(label= data$SampleLabel, size = 3.5, fontface="bold") +
  scale_colour_manual(values = c("AD" = "orangered", "Old" = "darkslateblue", "Young" = "gold")) +
  theme(legend.position = "top")
edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, t = 0.35, unit = 'cm'),
                                      vjust = -5,
                                      size = 14, face="bold"),
          axis.title.y = element_text(margin = margin(l = 1, r = 0.35, unit = 'cm'),
                                      vjust = 5,
                                      size = 14, face="bold"),
          axis.text.x = element_text(size = 13,
                                     colour = 'black'),
          axis.text.y = element_text(size = 13,
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
boxplot_counts <- boxplot_counts + theme(legend.key.size = unit(3, "lines"),
                                         legend.text = element_text(size = 15))

boxplot_counts
saveRDS(boxplot_counts, "./maprate_box_3each_new.rds")




#### cumulative gene counts
## all samples

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

boxplot_counts <- ggplot(total_counts_df, aes(x = condition, y = TotalCount)) +
  geom_boxplot(fill="white", alpha=0.15) +
  labs(x = "Sample Group",
       y = "Cumulative Gene Count",
       colour = "") +
  geom_signif(comparisons = list(c("AD","Old"),c("Old", "Young"), c("AD", "Young")),
              annotations = "ns", vjust = -0.5, 
              textsize = 6, colour = "black", y_position = 27900000, fontface="bold")

boxplot_counts <- boxplot_counts +
  geom_point(aes(x = condition, y = TotalCount, colour = condition), ## change shape: shape = condition
             size = 3) +
  geom_text_repel(label= total_counts_df$SampleLabel, size = 3.5, fontface="bold") +
  scale_colour_manual(values = c("AD" = "orangered", "Old" = "darkslateblue", "Young" = "gold")) +
  theme(legend.position = "top")

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, t = 0.35, unit = 'cm'),
                                      vjust = -5,
                                      size = 14, face="bold"),
          axis.title.y = element_text(margin = margin(l = 1, r = 0.35, unit = 'cm'),
                                      vjust = 5,
                                      size = 14, face="bold"),
          axis.text.x = element_text(size = 13,
                                     colour = 'black'),
          axis.text.y = element_text(size = 13,
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
boxplot_counts <- boxplot_counts + theme(legend.key.size = unit(3, "lines"),
                                         legend.text = element_text(size = 15))

boxplot_counts
saveRDS(boxplot_counts, "./cumulativegenecount_box_all_new.rds")


## selected samples

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
boxplot_counts <- ggplot(total_counts_df, aes(x = condition, y = TotalCount)) +
  geom_boxplot(fill="white", alpha=0.15) +
  labs(x = "Sample Group",
       y = "Cumulative Gene Count",
       colour = "") +
  geom_signif(comparisons = list(c("AD","Old"),c("Old", "Young"), c("AD", "Young")),
              annotations = "ns", vjust = -0.5, 
              textsize = 6, colour = "black", y_position = 27900000, fontface="bold")

boxplot_counts <- boxplot_counts +
  geom_point(aes(x = condition, y = TotalCount, colour = condition), ## change shape: shape = condition
             size = 3) +
  geom_text_repel(label= total_counts_df$SampleLabel, size = 3.5, fontface="bold") +
  scale_colour_manual(values = c("AD" = "orangered", "Old" = "darkslateblue", "Young" = "gold")) +
  theme(legend.position = "top")

edit_plots <- function(
    ggplot_object
) {
  ggplot_object <- ggplot_object +
    theme(axis.title.x = element_text(margin = margin(b = 1, t = 0.35, unit = 'cm'),
                                      vjust = -5,
                                      size = 14, face="bold"),
          axis.title.y = element_text(margin = margin(l = 1, r = 0.35, unit = 'cm'),
                                      vjust = 5,
                                      size = 14, face="bold"),
          axis.text.x = element_text(size = 13,
                                     colour = 'black'),
          axis.text.y = element_text(size = 13,
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
boxplot_counts <- boxplot_counts + theme(legend.key.size = unit(3, "lines"),
                                         legend.text = element_text(size = 15))

boxplot_counts
setwd("D://ABBY.windows_surface_pro/diffexp_results/")
saveRDS(boxplot_counts, "./cumulativegenecount_box_3each_new.rds")
