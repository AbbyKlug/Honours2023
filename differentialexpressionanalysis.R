library(ggplot2)
library(DESeq2)
library(tximport)
library(tidyverse)

## Yi's tx2gene
tx2gene = readRDS("D://ABBY.windows_surface_pro/DIFF-EXP/tx2gene.rds")
setwd("D://ABBY.windows_surface_pro/FASTQ.files/")
setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/FASTQ.files/")


###### AD vs Old ######

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
          './Old/19-11A-Old.fastq.gz_counts/quant.sf')

txi.salmon <- tximport(files = files, 
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)

sampleTable <- data.frame(condition = factor(c("AD", "AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","Old","Old","Old","Old","Old","Old","Old","Old","Old","Old")))
names <-  (c("AD1", "AD2","AD3","AD4","AD5","AD6","AD7","AD8","AD9","AD10","AD11","AD12",
             "Old1","Old2","Old3","Old4","Old5","Old6","Old7","Old8","Old9","Old10"))

rownames(sampleTable) = names
colnames(txi.salmon$counts) <- names

dds_ADvsOld <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~ condition)
dds_ADvsOld$condition <- relevel(dds_ADvsOld$condition, ref = "Old")
dds_ADvsOld <- DESeq(dds_ADvsOld)
resultsNames(dds_ADvsOld)

## PCA ##

vsd <- vst(dds_ADvsOld, blind = FALSE) # or varianceStabilizingTransformation instead of vst

plotPCA(vsd, intgroup="condition") +
  geom_label(aes(label = name))

## results ##

res_ADvsOld <- lfcShrink(dds_ADvsOld, coef="condition_AD_vs_Old", type="apeglm")
summary(res_ADvsOld)
sum(res_ADvsOld$padj < 0.05, na.rm = TRUE)

## saving ##

setwd("D://ABBY.windows_surface_pro/differentialexp/")
saveRDS(res_ADvsOld, file = "./DESEQres_pairwise_ADvsOld(all-bulk).rds")


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

sampleTable <- data.frame(condition = factor(c("Old","Old","Old","Old","Old","Old","Old","Old","Old","Old", "Young","Young","Young","Young","Young","Young","Young","Young")))
names <-  (c("Old1","Old2","Old3","Old4","Old5","Old6","Old7","Old8","Old9","Old10",
             "Young1","Young2","Young3","Young4","Young5","Young6","Young7","Young8"))

rownames(sampleTable) = names
colnames(txi.salmon$counts) <- names

dds_OldvsYoung <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~ condition)
dds_OldvsYoung$condition <- relevel(dds_OldvsYoung$condition, ref = "Young")
dds_OldvsYoung <- DESeq(dds_OldvsYoung)
resultsNames(dds_OldvsYoung)

## PCA ##

vsd <- vst(dds_OldvsYoung, blind = FALSE) # or varianceStabilizingTransformation instead of vst

plotPCA(vsd, intgroup="condition") +
  geom_label(aes(label = name))

## results ##

res_OldvsYoung <- lfcShrink(dds_OldvsYoung, coef="condition_Old_vs_Young", type="apeglm")
summary(res_OldvsYoung)
sum(res_OldvsYoung$padj < 0.05, na.rm = TRUE)

## saving ##

setwd("D://ABBY.windows_surface_pro/differentialexp/")
saveRDS(res_OldvsYoung, file = "./DESEQres_pairwise_OldvsYoung(all-bulk).rds")


###### AD vs Young ######

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

sampleTable <- data.frame(condition = factor(c("AD", "AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD", "Young","Young","Young","Young","Young","Young","Young","Young")))
names <-  (c("AD1", "AD2","AD3","AD4","AD5","AD6","AD7","AD8","AD9","AD10","AD11","AD12",
             "Young1","Young2","Young3","Young4","Young5","Young6","Young7","Young8"))


rownames(sampleTable) = names
colnames(txi.salmon$counts) <- names

dds_ADvsYoung <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~ condition)
dds_ADvsYoung$condition <- relevel(dds_ADvsYoung$condition, ref = "Young")
dds_ADvsYoung <- DESeq(dds_ADvsYoung)
resultsNames(dds_ADvsYoung)

## PCA ##

vsd <- vst(dds_ADvsYoung, blind = FALSE) # or varianceStabilizingTransformation instead of vst

plotPCA(vsd, intgroup="condition") +
  geom_label(aes(label = name))

## results ##

res_ADvsYoung <- lfcShrink(dds_ADvsYoung, coef="condition_AD_vs_Young", type="apeglm")
summary(res_ADvsYoung)
sum(res_ADvsYoung$padj < 0.05, na.rm = TRUE)

## saving ##

setwd("D://ABBY.windows_surface_pro/differentialexp/")
saveRDS(res_ADvsYoung, file = "./DESEQres_pairwise_ADvsYoung(all-bulk).rds")



####### DGE SELECTING 4 SAMPLES PER CONDITION #######
# AD 4,5,7,9 # Old 7,8,9,10 # Young 1,3,5,8

## Yi's tx2gene
tx2gene = readRDS("D://ABBY.windows_surface_pro/DIFF-EXP/tx2gene.rds")
setwd("D://ABBY.windows_surface_pro/FASTQ.files/")


# ## files = c('./AD/23-2A-AD.fastq.gz_counts/quant.sf',
#           './AD/24-3T-AD.fastq.gz_counts/quant.sf',
#           './AD/26-3A-AD.fastq.gz_counts/quant.sf',
#           './AD/28-8T-AD.fastq.gz_counts/quant.sf',
#           './Old/16-14T-Old.fastq.gz_counts/quant.sf',
#           './Old/17-9A-Old.fastq.gz_counts/quant.sf',
#           './Old/18-10A-Old.fastq.gz_counts/quant.sf',
#           './Old/19-11A-Old.fastq.gz_counts/quant.sf',
#           './Young/2-12A-Young.fastq.gz_counts/quant.sf',
#           './Young/4-13A-Young.fastq.gz_counts/quant.sf',
#           './Young/6-14A-Young.fastq.gz_counts/quant.sf',
#           './Young/9-16A-Young.fastq.gz_counts/quant.sf') 


###### AD vs Old ######

files = c('./AD/23-2A-AD.fastq.gz_counts/quant.sf',
          './AD/24-3T-AD.fastq.gz_counts/quant.sf',
          './AD/26-3A-AD.fastq.gz_counts/quant.sf',
          './AD/28-8T-AD.fastq.gz_counts/quant.sf',
          './Old/16-14T-Old.fastq.gz_counts/quant.sf',
          './Old/17-9A-Old.fastq.gz_counts/quant.sf',
          './Old/18-10A-Old.fastq.gz_counts/quant.sf',
          './Old/19-11A-Old.fastq.gz_counts/quant.sf')

txi.salmon <- tximport(files = files, 
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)

sampleTable <- data.frame(condition = factor(c("AD", "AD","AD","AD","Old","Old","Old","Old")))
names <-  (c("AD4","AD5","AD7","AD9",
             "Old7","Old8","Old9","Old10"))
rownames(sampleTable) = names
colnames(txi.salmon$counts) <- names

dds_ADvsOld <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~ condition)
dds_ADvsOld$condition <- relevel(dds_ADvsOld$condition, ref = "Old")
dds_ADvsOld <- DESeq(dds_ADvsOld)
resultsNames(dds_ADvsOld)

## PCA ##

vsd <- vst(dds_ADvsOld, blind = FALSE) # or varianceStabilizingTransformation instead of vst

plotPCA(vsd, intgroup="condition") +
  geom_label(aes(label = name))

## results ##

res_ADvsOld <- lfcShrink(dds_ADvsOld, coef="condition_AD_vs_Old", type="apeglm")
summary(res_ADvsOld)
sum(res_ADvsOld$padj < 0.05, na.rm = TRUE)

## saving ##

setwd("D://ABBY.windows_surface_pro/differentialexp/")
saveRDS(res_ADvsOld, file = "./DESEQres_pairwise_ADvsOld(4each-bulk).rds")



###### Old vs Young ######

files = c('./Old/16-14T-Old.fastq.gz_counts/quant.sf',
          './Old/17-9A-Old.fastq.gz_counts/quant.sf',
          './Old/18-10A-Old.fastq.gz_counts/quant.sf',
          './Old/19-11A-Old.fastq.gz_counts/quant.sf',
          './Young/2-12A-Young.fastq.gz_counts/quant.sf',
          './Young/4-13A-Young.fastq.gz_counts/quant.sf',
          './Young/6-14A-Young.fastq.gz_counts/quant.sf',
          './Young/9-16A-Young.fastq.gz_counts/quant.sf')

txi.salmon <- tximport(files = files, 
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)

sampleTable <- data.frame(condition = factor(c("Old","Old","Old","Old", "Young","Young","Young","Young")))
names <-  (c("Old7","Old8","Old9","Old10", 
             "Young1","Young3","Young5","Young8"))
rownames(sampleTable) = names
colnames(txi.salmon$counts) <- names

dds_OldvsYoung <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~ condition)
dds_OldvsYoung$condition <- relevel(dds_OldvsYoung$condition, ref = "Young")
dds_OldvsYoung <- DESeq(dds_OldvsYoung)
resultsNames(dds_OldvsYoung)

## PCA ##

vsd <- vst(dds_OldvsYoung, blind = FALSE) # or varianceStabilizingTransformation instead of vst

plotPCA(vsd, intgroup="condition") +
  geom_label(aes(label = name))

## results ##

res_OldvsYoung <- lfcShrink(dds_OldvsYoung, coef="condition_Old_vs_Young", type="apeglm")
summary(res_OldvsYoung)
sum(res_OldvsYoung$padj < 0.05, na.rm = TRUE)

## saving ##

setwd("D://ABBY.windows_surface_pro/differentialexp/")
saveRDS(res_OldvsYoung, file = "./DESEQres_pairwise_OldvsYoung(4each-bulk).rds")


###### AD vs Young ######

files = c('./AD/23-2A-AD.fastq.gz_counts/quant.sf',
          './AD/24-3T-AD.fastq.gz_counts/quant.sf',
          './AD/26-3A-AD.fastq.gz_counts/quant.sf',
          './AD/28-8T-AD.fastq.gz_counts/quant.sf',
          './Young/2-12A-Young.fastq.gz_counts/quant.sf',
          './Young/4-13A-Young.fastq.gz_counts/quant.sf',
          './Young/6-14A-Young.fastq.gz_counts/quant.sf',
          './Young/9-16A-Young.fastq.gz_counts/quant.sf')

txi.salmon <- tximport(files = files, 
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)


sampleTable <- data.frame(condition = factor(c("AD","AD","AD","AD", "Young","Young","Young","Young")))
names <-  (c("AD4","AD5","AD7","AD9", 
             "Young1","Young3","Young5","Young8"))
rownames(sampleTable) = names
colnames(txi.salmon$counts) <- names

dds_ADvsYoung <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~ condition)
dds_ADvsYoung$condition <- relevel(dds_ADvsYoung$condition, ref = "Young")
dds_ADvsYoung <- DESeq(dds_ADvsYoung)
resultsNames(dds_ADvsYoung)

## PCA ##

vsd <- vst(dds_ADvsYoung, blind = FALSE) # or varianceStabilizingTransformation instead of vst

plotPCA(vsd, intgroup="condition") +
  geom_label(aes(label = name))

## results ##

res_ADvsYoung <- lfcShrink(dds_ADvsYoung, coef="condition_AD_vs_Young", type="apeglm")
summary(res_ADvsYoung)
sum(res_ADvsYoung$padj < 0.05, na.rm = TRUE)

## saving ##

setwd("D://ABBY.windows_surface_pro/differentialexp/")
saveRDS(res_ADvsYoung, file = "./DESEQres_pairwise_ADvsYoung(4each-bulk).rds")


###### PLOT TOTAL GENE COUNTS PER SAMPLE ######

setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp/")
tx2gene = readRDS("./tx2gene.rds")
setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/FASTQ.files/")
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
colnames(txi.salmon$abundance) <- names

counts_df <- as.data.frame(txi.salmon$counts)

# Calculate total gene counts per sample
totals <- colSums(counts_df)
samples <- names(totals)
condition <- c("AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD",
               "Old", "Old","Old","Old","Old","Old","Old","Old","Old","Old",
               "Young", "Young","Young","Young","Young","Young","Young","Young")

# Convert the result to a data frame
total_counts_df <- data.frame(Sample = samples, TotalCount = totals, condition)

# Create a bar plot using ggplot2
plot_counts <- ggplot(total_counts_df, aes(x = Sample, y = TotalCount, fill = condition)) +
  geom_bar(stat = "identity") +
  labs(title = "Total Gene Counts per Sample", x = "Sample", y = "Total Gene Count") +
  scale_fill_manual(values = c("AD" = "slategrey", "Old" = "slategray2", "Young" = "lightblue3"))

plot_counts
setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp/")
ggsave("./plot_totalgenecounts_persample.png", plot_counts, width = 15, height = 10)


# Calculate total abundances per sample
#### The tximport package has a single function for importing transcript-level estimates. The type argument is used to specify what software was used for estimation. 
## A simple list with matrices, "abundance", "counts", and "length", is returned, where the transcript level information is summarized to the gene-level. 
## Typically, abundance is provided by the quantification tools as TPM (transcripts-per-million)
## the counts are estimated counts (possibly fractional)
## the "length" matrix contains the effective gene lengths. The "length" matrix can be used to generate an offset matrix for downstream gene-level differential analysis of count matrices

abundances_df <- as.data.frame(txi.salmon$abundance)
total_ab <- colSums(abundances_df)
samples_ab <- names(total_ab)
condition <- c("AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD","AD",
               "Old", "Old","Old","Old","Old","Old","Old","Old","Old","Old",
               "Young", "Young","Young","Young","Young","Young","Young","Young")

# Convert the result to a data frame
total_abundances_df <- data.frame(Sample = samples_ab, TotalAbundance = total_ab, condition)

# Create a bar plot using ggplot2
plot_abundances <- ggplot(total_abundances_df, aes(x = Sample, y = TotalAbundance, fill = condition)) +
  geom_bar(stat = "identity") +
  labs(title = "Total Gene Abundance per Sample", x = "Sample", y = "Total Gene Abundance") +
  scale_fill_manual(values = c("AD" = "mistyrose4", "Old" = "mistyrose3", "Young" = "mistyrose2"))

plot_abundances
setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp/")
ggsave("./plot_totalgeneabundances(normalised)_persample.png", plot_abundances, width = 15, height = 10)


### checking ###
colSums(txi.salmon[["abundance"]])


###### PLOT TOTAL GENE COUNTS PER SAMPLE (selected samples) ######

setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp/")
tx2gene = readRDS("./tx2gene.rds")
setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/FASTQ.files/")
files = c('./AD/23-2A-AD.fastq.gz_counts/quant.sf',
          './AD/24-3T-AD.fastq.gz_counts/quant.sf',
          './AD/26-3A-AD.fastq.gz_counts/quant.sf',
          './AD/28-8T-AD.fastq.gz_counts/quant.sf',
          './Old/16-14T-Old.fastq.gz_counts/quant.sf',
          './Old/17-9A-Old.fastq.gz_counts/quant.sf',
          './Old/18-10A-Old.fastq.gz_counts/quant.sf',
          './Old/19-11A-Old.fastq.gz_counts/quant.sf',
          './Young/2-12A-Young.fastq.gz_counts/quant.sf',
          './Young/4-13A-Young.fastq.gz_counts/quant.sf',
          './Young/6-14A-Young.fastq.gz_counts/quant.sf',
          './Young/9-16A-Young.fastq.gz_counts/quant.sf')

txi.salmon <- tximport(files = files, 
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)

sampleTable <- data.frame(condition = factor(c("AD", "AD","AD","AD","Old","Old","Old","Old", "Young","Young","Young","Young")))
names <-  c("AD4","AD5","AD7","AD9",
            "Old7","Old8","Old9","Old10",
            "Young1","Young3","Young5","Young8")
rownames(sampleTable) = names
colnames(txi.salmon$counts) <- names
colnames(txi.salmon$abundance) <- names

counts_df <- as.data.frame(txi.salmon$counts)

# Calculate total gene counts per sample
totals <- colSums(counts_df)
samples <- names(totals)
condition <- c("AD","AD","AD","AD",
               "Old", "Old","Old","Old",
               "Young", "Young","Young","Young")

# Convert the result to a data frame
total_counts_df <- data.frame(Sample = samples, TotalCount = totals, condition)

# Create a bar plot using ggplot2
plot_counts <- ggplot(total_counts_df, aes(x = Sample, y = TotalCount, fill = condition)) +
  geom_bar(stat = "identity") +
  labs(title = "Total Gene Counts per Sample", x = "Sample", y = "Total Gene Count") +
  scale_fill_manual(values = c("AD" = "slategrey", "Old" = "slategray2", "Young" = "lightblue3")) +
  theme_minimal()

setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp/")
ggsave("./plot_totalgenecounts_persample_selectedsamples.png", plot_counts, width = 15, height = 10)

#### boxlot ####

library(ggsignif)
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
# Calculate total abundances per sample
#### The tximport package has a single function for importing transcript-level estimates. The type argument is used to specify what software was used for estimation. 
## A simple list with matrices, "abundance", "counts", and "length", is returned, where the transcript level information is summarized to the gene-level. 
## Typically, abundance is provided by the quantification tools as TPM (transcripts-per-million)
## the counts are estimated counts (possibly fractional)
## the "length" matrix contains the effective gene lengths. The "length" matrix can be used to generate an offset matrix for downstream gene-level differential analysis of count matrices

abundances_df <- as.data.frame(txi.salmon$abundance)
total_ab <- colSums(abundances_df)
samples_ab <- names(total_ab)
condition <- c("AD","AD","AD","AD",
               "Old", "Old","Old","Old",
               "Young", "Young","Young","Young")

# Convert the result to a data frame
total_abundances_df <- data.frame(Sample = samples_ab, TotalAbundance = total_ab, condition)

# Create a bar plot using ggplot2
plot_abundances <- ggplot(total_abundances_df, aes(x = Sample, y = TotalAbundance, fill = condition)) +
  geom_bar(stat = "identity") +
  labs(title = "Total Gene Abundance per Sample", x = "Sample", y = "Total Gene Abundance") +
  scale_fill_manual(values = c("AD" = "mistyrose4", "Old" = "mistyrose3", "Young" = "mistyrose2")) +
  theme_minimal()

plot_abundances
setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp/")
ggsave("./plot_totalgeneabundances(normalised)_persample_selectedsamples.png", plot_abundances, width = 15, height = 10)


###### Looking at DESeq2 results ######
setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp/")
setwd("D://ABBY.windows_surface_pro/differentialexp/")
## AD vs Old ## 

res_ADvsOld <- readRDS("./DESEQres_pairwise_ADvsOld(4each-bulk).rds")
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

filtered_ADvsOld_LFC1 <- ADvsOld_df[ADvsOld_df$basemean > 10, ]

## volcano plot no LFC cutoff
ggplot(filtered_ADvsOld_LFC1, aes(x = LFC, y = -log10(padj), color = filtered_ADvsOld_LFC1$padj < 0.05)) +
  geom_point() +
  scale_color_manual(values = c("slategrey", "springgreen")) +  # Color for non-significant and significant genes
  theme_minimal() +
  labs(
    x = "Log Fold Change (LFC)",
    y = "-log10(Adjusted p-value)"
  )

## volcano plot with LFC cutoff
filtered_ADvsOld_LFC1 <- filtered_ADvsOld_LFC1[filtered_ADvsOld_LFC1$LFC > 1 | filtered_ADvsOld_LFC1$LFC < -1, ]
ggplot(filtered_ADvsOld_LFC1, aes(x = LFC, y = -log10(padj), color = filtered_ADvsOld_LFC1$padj < 0.05)) +
  geom_point() +
  scale_color_manual(values = c("slategrey", "springgreen")) +  # Color for non-significant and significant genes
  theme_minimal() +
  labs(
    x = "Log Fold Change (LFC)",
    y = "-log10(Adjusted p-value)"
  )


filtered_ADvsOld_LFC1<- filtered_ADvsOld_LFC1[filtered_ADvsOld_LFC1$padj < 0.05, ]
filtered_ADvsOld_LFC1 <- filtered_ADvsOld_LFC1[filtered_ADvsOld_LFC1$LFC > 1 | filtered_ADvsOld_LFC1$LFC < -1, ]
filtered_ADvsOld_LFC1 <- na.omit(filtered_ADvsOld_LFC1)
saveRDS(filtered_ADvsOld_LFC1, "./dataframe_filtered_ADvsOld_LFC1_basemean10_padj0.05.rds")

## Old vs Young
res_OldvsYoung <- readRDS("./DESEQres_pairwise_OldvsYoung(4each-bulk).rds")
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

## volcano plot no LFC cutoff
filtered_OldvsYoung_LFC1 <- OldvsYoung_df[OldvsYoung_df$basemean > 10, ]
ggplot(filtered_OldvsYoung_LFC1, aes(x = LFC, y = -log10(padj), color = filtered_OldvsYoung_LFC1$padj < 0.05)) +
  geom_point() +
  scale_color_manual(values = c("slategrey", "springgreen")) +  # Color for non-significant and significant genes
  theme_minimal() +
  labs(
    x = "Log Fold Change (LFC)",
    y = "-log10(Adjusted p-value)"
  )

## volcano plot with LFC cutoff
filtered_OldvsYoung_LFC1 <- filtered_OldvsYoung_LFC1[filtered_OldvsYoung_LFC1$LFC > 1 | filtered_OldvsYoung_LFC1$LFC < -1, ]
ggplot(filtered_OldvsYoung_LFC1, aes(x = LFC, y = -log10(padj), color = filtered_OldvsYoung_LFC1$padj < 0.05)) +
  geom_point() +
  scale_color_manual(values = c("slategrey", "springgreen")) +  # Color for non-significant and significant genes
  theme_minimal() +
  labs(
    x = "Log Fold Change (LFC)",
    y = "-log10(Adjusted p-value)"
  )

filtered_OldvsYoung_LFC1<- filtered_OldvsYoung_LFC1[filtered_OldvsYoung_LFC1$padj < 0.05, ]
filtered_OldvsYoung_LFC1 <- filtered_OldvsYoung_LFC1[filtered_OldvsYoung_LFC1$LFC > 1 | filtered_OldvsYoung_LFC1$LFC < -1, ]
filtered_OldvsYoung_LFC1 <- na.omit(filtered_OldvsYoung_LFC1)
saveRDS(filtered_OldvsYoung_LFC1, "./dataframe_filtered_OldvsYoung_LFC1_basemean10_padj0.05.rds")


## AD vs YOUNG
res_ADvsYoung <- readRDS("./DESEQres_pairwise_ADvsYoung(4each-bulk).rds")
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

## volcano plot no LFC cutoff
filtered_ADvsYoung_LFC1 <- ADvsYoung_df[ADvsYoung_df$basemean > 10, ]
ggplot(filtered_ADvsYoung_LFC1, aes(x = LFC, y = -log10(padj), color = filtered_ADvsYoung_LFC1$padj < 0.05)) +
  geom_point() +
  scale_color_manual(values = c("slategrey", "springgreen")) +  # Color for non-significant and significant genes
  theme_minimal() +
  labs(
    x = "Log Fold Change (LFC)",
    y = "-log10(Adjusted p-value)"
  )

## volcano plot with LFC cutoff
filtered_ADvsYoung_LFC1 <- filtered_ADvsYoung_LFC1[filtered_ADvsYoung_LFC1$LFC > 1 | filtered_ADvsYoung_LFC1$LFC < -1, ]
ggplot(filtered_ADvsYoung_LFC1, aes(x = LFC, y = -log10(padj), color = filtered_ADvsYoung_LFC1$padj < 0.05)) +
  geom_point() +
  scale_color_manual(values = c("slategrey", "springgreen")) +  # Color for non-significant and significant genes
  theme_minimal() +
  labs(
    x = "Log Fold Change (LFC)",
    y = "-log10(Adjusted p-value)"
  )


filtered_ADvsYoung_LFC1<- filtered_ADvsYoung_LFC1[filtered_ADvsYoung_LFC1$padj < 0.05, ]
filtered_ADvsYoung_LFC1 <- filtered_ADvsYoung_LFC1[filtered_ADvsYoung_LFC1$LFC > 1 | filtered_ADvsYoung_LFC1$LFC < -1, ]
filtered_ADvsYoung_LFC1 <- na.omit(filtered_ADvsYoung_LFC1)
saveRDS(filtered_ADvsYoung_LFC1, "./dataframe_filtered_ADvsYoung_LFC1_basemean10_padj0.05.rds")


##### CLUSTERPROFILER #####
## GENES i put in = all filters applied
## BACKGROUND/UNIVERSE = only base mean filter

library(clusterProfiler)
setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp/")
#setwd("D://ABBY.windows_surface_pro/differentialexp/")

### UNIVERSES ###
## AD vs OLD
res_ADvsOld <- readRDS("./DESEQres_pairwise_ADvsOld(4each-bulk).rds")
ADvsOld_df <- data.frame(Gene = res_ADvsOld@rownames, 
                         LFC = res_ADvsOld@listData[["log2FoldChange"]],
                         pvalue = res_ADvsOld@listData[["pvalue"]],
                         padj = res_ADvsOld@listData[["padj"]],
                         basemean = res_ADvsOld@listData[["baseMean"]])
filtered_ADvsOld_basemean <- ADvsOld_df[ADvsOld_df$basemean > 10, ]

## OLD vs YOUNG
res_OldvsYoung <- readRDS("./DESEQres_pairwise_OldvsYoung(4each-bulk).rds")
OldvsYoung_df <- data.frame(Gene = res_OldvsYoung@rownames, 
                            LFC = res_OldvsYoung@listData[["log2FoldChange"]],
                            pvalue = res_OldvsYoung@listData[["pvalue"]],
                            padj = res_OldvsYoung@listData[["padj"]],
                            basemean = res_OldvsYoung@listData[["baseMean"]])
filtered_OldvsYoung_basemean <- OldvsYoung_df[OldvsYoung_df$basemean > 10, ]

## AD vs YOUNG
res_ADvsYoung <- readRDS("./DESEQres_pairwise_ADvsYoung(4each-bulk).rds")
ADvsYoung_df <- data.frame(Gene = res_ADvsYoung@rownames, 
                           LFC = res_ADvsYoung@listData[["log2FoldChange"]],
                           pvalue = res_ADvsYoung@listData[["pvalue"]],
                           padj = res_ADvsYoung@listData[["padj"]],
                           basemean = res_ADvsYoung@listData[["baseMean"]])
filtered_ADvsYoung_basemean <- ADvsYoung_df[ADvsYoung_df$basemean > 10, ]


## INPUT GENES ##
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


### AD vs Old ###

## no background
GOenrich_ADOld <- enrichGO(ADvsOld_genes, 'org.Hs.eg.db', keyType = "ENSEMBL", ont = "BP" ,pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
## + background
GOenrich_ADOld2 <- enrichGO(ADvsOld_genes, 'org.Hs.eg.db', keyType = "ENSEMBL", ont = "BP" ,pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05,
                            readable = TRUE,universe = filtered_ADvsOld_basemean$Gene) ##filtered basemean not LFC or p value == what IS expressed in cell

saveRDS(GOenrich_ADOld, "./recent_GOenrich_ADOld_nobackground.rds")
saveRDS(GOenrich_ADOld2, "./recent_GOenrich_ADOld_+background.rds")

### Old vs Young ###

## no background
GOenrich_OldYoung <- enrichGO(OldvsYoung_genes, 'org.Hs.eg.db', keyType = "ENSEMBL", ont = "BP" ,pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05,readable = TRUE)
## + background
GOenrich_OldYoung2 <- enrichGO(OldvsYoung_genes, 'org.Hs.eg.db', keyType = "ENSEMBL", ont = "BP" ,pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05,
                               readable = TRUE,universe = filtered_OldvsYoung_basemean$Gene) ##filtered basemean not LFC or p value == what IS expressed in cell

saveRDS(GOenrich_OldYoung, "./recent_GOenrich_OldYoung_nobackground.rds")
saveRDS(GOenrich_OldYoung2, "./recent_GOenrich_OldYoung_+background.rds")

### AD vs Young ###

## no background
GOenrich_ADYoung <- enrichGO(ADvsYoung_genes, 'org.Hs.eg.db', keyType = "ENSEMBL", ont = "BP" ,pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05,readable = TRUE)
## + background
GOenrich_ADYoung2 <- enrichGO(ADvsYoung_genes, 'org.Hs.eg.db', keyType = "ENSEMBL", ont = "BP" ,pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05,
                              readable = TRUE,universe = filtered_ADvsYoung_basemean$Gene) ##filtered basemean not LFC or p value == what IS expressed in cell

saveRDS(GOenrich_ADYoung, "./recent_GOenrich_ADYoung_nobackground.rds")
saveRDS(GOenrich_ADYoung2, "./recent_GOenrich_ADYoung_+background.rds")



### filtered DESeq2 genes ###
setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp/")

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
BiocManager::install("VennDiagram")
library(VennDiagram)

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
  cat.dist = 0.03,  # Distance of labels from the sets
  cat.cex = 0.8,   # Label text size
  ext.text = TRUE,
  label.col = "black",
  fill = c("olivedrab1", "limegreen", "palegreen"),  # Colors for the sets
  alpha = 0.5
)

grid.newpage()
grid.draw(venn_result_up)

OldYoung_genes_up %in% ADOld_genes_up

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
  cat.pos = c(0, 0, 0),  # Place labels at default position (inside the sets)
  cat.dist = 0.03,  # Distance of labels from the sets
  cat.cex = 0.8,   # Label text size
  ext.text = TRUE,
  label.col = "black",
  fill = c("red", "indianred3", "salmon"),  # Colors for the sets
  alpha = 0.5
)

grid.newpage()
grid.draw(venn_result_down)

# COMMON UPREGULATED
common_ADOld_ADYoung_up <- intersect(ADOld_genes_up, ADYoung_genes_up)
common_OldYoung_ADYoung_up <- intersect(OldYoung_genes_up, ADYoung_genes_up)
common_ADOld_OldYoung_up <- intersect(ADOld_genes_up, OldYoung_genes_up)

saveRDS(common_ADOld_ADYoung_up, "./common_ADOld_ADYoung_up.rds")
saveRDS(common_OldYoung_ADYoung_up, "./common_OldYoung_ADYoung_up.rds")
saveRDS(common_ADOld_OldYoung_up, "./common_ADOld_OldYoung_up.rds")


## DEGs IN OLD/YOUNG+AD/YOUNG TO EXCLUDE FROM DEGs IN AD/OLD+AD/YOUNG UPregulated
genes_exclude_up <- intersect(common_ADOld_ADYoung_up, common_OldYoung_ADYoung_up)
## NONE
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
## NONE
common_OldYoung_ADYoung_down %in% common_ADOld_ADYoung_down

# not necessary:
# saveRDS(genes_exclude_up, "./genes_exclude_up.rds")
# saveRDS(genes_exclude_down, "./genes_exclude_down.rds")





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


ggplot(aes(x = respondent_wall_type)) +
  geom_bar(aes(fill = village), position = "dodge")

### PLOTTING UP- AND DOWN-REGULATED DEGs PER PAIRWISE COMPARISON ###
### filtered DESeq2 genes ###
setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp/")

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

saveRDS(total_counts_df, "./dataframe_UP-DOWN_DEGs_percomparison_forplotting.rds")


### PLOTTING CLUSTERPROFILER ###
setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp/")
enrich_ADOld <- readRDS("./recent_GOenrich_ADOld_+background.rds")
enrich_OldYoung <- readRDS("./recent_GOenrich_OldYoung_+background.rds")
enrich_ADYoung <- readRDS("./recent_GOenrich_ADYoung_+background.rds")

library(enrichplot)
## barplot() displays most significant or selected enriched terms
barplot(enrich_ADOld, showCategory=20, label_format = 50)
barplot(enrich_OldYoung, showCategory=20, label_format = 50)
barplot(enrich_ADYoung, showCategory=20, label_format = 50)

?par
??enrichplot

## results dataframes
res_enrich_ADOld <- as.data.frame(enrich_ADOld@result)
res_enrich_OldYoung <- as.data.frame(enrich_OldYoung@result)
res_enrich_ADYoung <- as.data.frame(enrich_ADYoung@result)

## filtering top 20 BPs

### ordering by increasing p.adjust ###
# Column to order by
column_to_order <- "Count"
# Order the data frame based on the specified column
ordered_res_enrich_ADOld <- res_enrich_ADOld[order(res_enrich_ADOld[[column_to_order]], decreasing = TRUE), ]
ordered_res_enrich_OldYoung <- res_enrich_OldYoung[order(res_enrich_OldYoung[[column_to_order]], decreasing = TRUE), ]
ordered_res_enrich_ADYoung <- res_enrich_ADYoung[order(res_enrich_ADYoung[[column_to_order]], decreasing = TRUE), ]


ordered_res_enrich_ADOld <- ordered_res_enrich_ADOld[1:20, ]
ordered_res_enrich_OldYoung <- ordered_res_enrich_OldYoung[1:20, ]
ordered_res_enrich_ADYoung <- ordered_res_enrich_ADYoung[1:20, ]


df_enrich_ADOld <- data.frame(BP = ordered_res_enrich_ADOld$Description, padj = ordered_res_enrich_ADOld$p.adjust, counts = ordered_res_enrich_ADOld$Count)
df_enrich_OldYoung <- data.frame(BP = ordered_res_enrich_OldYoung$Description, padj = ordered_res_enrich_OldYoung$p.adjust, counts = ordered_res_enrich_OldYoung$Count)
df_enrich_ADYoung <- data.frame(BP = ordered_res_enrich_ADYoung$Description, padj = ordered_res_enrich_ADYoung$p.adjust, counts = ordered_res_enrich_ADYoung$Count)

### manual AD vs OLD
plot_enrich_ADOld <- ggplot(df_enrich_ADOld, aes(x = counts, y = BP), horiz = TRUE) + # horiz = TRUE horizontal bars
  geom_bar(stat = "identity", fill = "hotpink") +  
  labs(title = "Significant Enriched Terms AD vs Old", x = "Count", y = "Biological Process")

plot_enrich_ADOld

### manual OLD vs YOUNG
library(stringr)
plot_enrich_OldYoung <- ggplot(df_enrich_OldYoung, aes(x = counts, y = str_wrap(BP, width = 50)), horiz = TRUE) + # horiz = TRUE horizontal bars
  geom_bar(stat = "identity", fill = "hotpink") +  
  labs(title = "Significant Enriched Terms Old vs Young", x = "Count", y = "Biological Process")

plot_enrich_OldYoung

## manual AD vs YOUNG
plot_enrich_ADYoung <- ggplot(df_enrich_ADYoung, aes(x = counts, y = BP), horiz = TRUE) + # horiz = TRUE horizontal bars
  geom_bar(stat = "identity", fill = "hotpink") +  
  labs(title = "Significant Enriched Terms AD vs Young", x = "Count", y = "Biological Process")

plot_enrich_ADYoung


#### research/looking @ DEGs ####
### filtered DESeq2 genes ###
setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp/")

ADvsOld <- readRDS("./dataframe_filtered_ADvsOld_LFC1_basemean10_padj0.05.rds")
OldvsYoung <- readRDS("./dataframe_filtered_OldvsYoung_LFC1_basemean10_padj0.05.rds")
ADvsYoung <- readRDS("./dataframe_filtered_ADvsYoung_LFC1_basemean10_padj0.05.rds")

### ordering by increasing padj ###
# Column to order by
column_to_order <- "padj" ## change to LFC
# Order the data frame based on the specified column
ordered_ADvsOld <- ADvsOld[order(ADvsOld[[column_to_order]]), ]
ordered_OldvsYoung <- OldvsYoung[order(OldvsYoung[[column_to_order]]), ]
ordered_ADvsYoung <- ADvsYoung[order(ADvsYoung[[column_to_order]]), ]

saveRDS(ordered_ADvsOld, "./orderedbyPADJ_dataframe_filtered_ADvsOld_LFC1_basemean10_padj0.05.rds")
saveRDS(ordered_OldvsYoung, "./orderedbyPADJ_dataframe_filtered_OldvsYoung_LFC1_basemean10_padj0.05.rds")
saveRDS(ordered_ADvsYoung, "./orderedbyPADJ_dataframe_filtered_ADvsYoung_LFC1_basemean10_padj0.05.rds")





###### Looking at DESeq2 results ######
### label volcano plots with genes LFC > 4 or <-4 ###
setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp/")

## AD vs Old ## 
res_ADvsOld <- readRDS("./DESEQres_pairwise_ADvsOld(4each-bulk).rds")
ADvsOld_df <- data.frame(Gene = res_ADvsOld@rownames, 
                         LFC = res_ADvsOld@listData[["log2FoldChange"]],
                         pvalue = res_ADvsOld@listData[["pvalue"]],
                         padj = res_ADvsOld@listData[["padj"]],
                         basemean = res_ADvsOld@listData[["baseMean"]])

filtered_ADvsOld_basemean <- ADvsOld_df[ADvsOld_df$basemean > 10, ]
filtered_ADvsOld_basemean <- na.omit(filtered_ADvsOld_basemean)

## volcano plot no LFC cutoff
ggplot(filtered_ADvsOld_basemean, aes(x = LFC, y = -log10(padj), color = filtered_ADvsOld_basemean$padj < 0.05 & abs(filtered_ADvsOld_basemean$LFC) > 1)) +
  geom_point() +
  scale_color_manual(values = c("slategrey", "springgreen")) +  # Color for non-significant and significant genes
  theme_minimal() +
  labs(
    x = "Log Fold Change (LFC)",
    y = "-log10(Adjusted p-value)",
    colour = "Differentially Expressed Genes"
  ) +
  geom_text_repel(label= ifelse(abs(filtered_ADvsOld_basemean$LFC) > 4, filtered_ADvsOld_basemean$Gene, NA))

sum(abs(filtered_ADvsOld_basemean$LFC) > 4, na.rm = TRUE)


## Old vs Young
res_OldvsYoung <- readRDS("./DESEQres_pairwise_OldvsYoung(4each-bulk).rds")
OldvsYoung_df <- data.frame(Gene = res_OldvsYoung@rownames, 
                            LFC = res_OldvsYoung@listData[["log2FoldChange"]],
                            pvalue = res_OldvsYoung@listData[["pvalue"]],
                            padj = res_OldvsYoung@listData[["padj"]],
                            basemean = res_OldvsYoung@listData[["baseMean"]])
filtered_OldvsYoung_basemean <- OldvsYoung_df[OldvsYoung_df$basemean > 10, ]
filtered_OldvsYoung_basemean <- na.omit(filtered_OldvsYoung_basemean)

## volcano plot no LFC cutoff
ggplot(filtered_OldvsYoung_basemean, aes(x = LFC, y = -log10(padj), color = filtered_OldvsYoung_basemean$padj < 0.05 & abs(filtered_OldvsYoung_basemean$LFC) > 1)) +
  geom_point() +
  scale_color_manual(values = c("slategrey", "springgreen")) +  # Color for non-significant and significant genes
  theme_minimal() +
  labs(
    x = "Log Fold Change (LFC)",
    y = "-log10(Adjusted p-value)",
    colour = "Differentially Expressed Genes"
  ) +
  geom_text_repel(label= ifelse(abs(filtered_OldvsYoung_basemean$LFC) > 4, filtered_OldvsYoung_basemean$Gene, NA))

sum(abs(filtered_OldvsYoung_basemean$LFC) > 4, na.rm = TRUE)

## AD vs Young
res_ADvsYoung <- readRDS("./DESEQres_pairwise_ADvsYoung(4each-bulk).rds")
ADvsYoung_df <- data.frame(Gene = res_ADvsYoung@rownames, 
                           LFC = res_ADvsYoung@listData[["log2FoldChange"]],
                           pvalue = res_ADvsYoung@listData[["pvalue"]],
                           padj = res_ADvsYoung@listData[["padj"]],
                           basemean = res_ADvsYoung@listData[["baseMean"]])
filtered_ADvsYoung_basemean <- ADvsYoung_df[ADvsYoung_df$basemean > 10, ]
filtered_ADvsYoung_basemean <- na.omit(filtered_ADvsYoung_basemean)

## volcano plot no LFC cutoff
ggplot(filtered_ADvsYoung_basemean, aes(x = LFC, y = -log10(padj), color = filtered_ADvsYoung_basemean$padj < 0.05 & abs(filtered_ADvsYoung_basemean$LFC) > 1)) +
  geom_point() +
  scale_color_manual(values = c("slategrey", "springgreen")) +  # Color for non-significant and significant genes
  theme_minimal() +
  labs(
    x = "Log Fold Change (LFC)",
    y = "-log10(Adjusted p-value)",
    colour = "Differentially Expressed Genes"
  ) +
  geom_text_repel(label= ifelse(abs(filtered_ADvsYoung_basemean$LFC) > 4, filtered_ADvsYoung_basemean$Gene, NA))

sum(abs(filtered_ADvsYoung_basemean$LFC) > 4, na.rm = TRUE)


genes2lookat <- filtered_ADvsYoung_basemean[(abs(filtered_ADvsYoung_basemean$LFC) > 4), ]
saveRDS(genes2lookat, "./genes2lookat_ADvsYoung.rds")



###### home ######
setwd("D:/ABBY.windows_surface_pro/differentialexp/")
genes_interest <- readRDS("./genes2lookat_ADvsYoung.rds")
genes_interest$Gene






































































# files = c('./AD/20-1T-AD.fastq.gz_counts/quant.sf', 1
#           './AD/21-1A-AD.fastq.gz_counts/quant.sf', 2
#           './AD/22-2T-AD.fastq.gz_counts/quant.sf', 3
#           './AD/23-2A-AD.fastq.gz_counts/quant.sf', 4
#           './AD/24-3T-AD.fastq.gz_counts/quant.sf', 5
#           './AD/25-5T-AD.fastq.gz_counts/quant.sf', 6
#           './AD/26-3A-AD.fastq.gz_counts/quant.sf', 7
#           './AD/27-5A-AD.fastq.gz_counts/quant.sf', 8
#           './AD/28-8T-AD.fastq.gz_counts/quant.sf', 9
#           './AD/29-6T-AD.fastq.gz_counts/quant.sf', 10
#           './AD/30-9T-AD.fastq.gz_counts/quant.sf', 11
#           './AD/31-7T-AD.fastq.gz_counts/quant.sf', 12
#           './Old/10-8A-Old.fastq.gz_counts/quant.sf', 1
#           './Old/11-10T-Old.fastq.gz_counts/quant.sf',2
#           './Old/12-6A-Old.fastq.gz_counts/quant.sf', 3
#           './Old/13-11T-Old.fastq.gz_counts/quant.sf', 4
#           './Old/14-7A-Old.fastq.gz_counts/quant.sf', 5
#           './Old/15-13T-Old.fastq.gz_counts/quant.sf', 6
#           './Old/16-14T-Old.fastq.gz_counts/quant.sf', 7
#           './Old/17-9A-Old.fastq.gz_counts/quant.sf', 8
#           './Old/18-10A-Old.fastq.gz_counts/quant.sf', 9
#           './Old/19-11A-Old.fastq.gz_counts/quant.sf', 10
#           './Young/2-12A-Young.fastq.gz_counts/quant.sf', 1
#           './Young/3-17T-Young.fastq.gz_counts/quant.sf', 2
#           './Young/4-13A-Young.fastq.gz_counts/quant.sf', 3
#           './Young/5-18T-Young.fastq.gz_counts/quant.sf', 4
#           './Young/6-14A-Young.fastq.gz_counts/quant.sf', 5
#           './Young/7-19T-Young.fastq.gz_counts/quant.sf', 6
#           './Young/8-15A-Young.fastq.gz_counts/quant.sf', 7
#           './Young/9-16A-Young.fastq.gz_counts/quant.sf') 8



# # Calculate total gene counts per sample
# total_counts <- colSums(counts_df)
# 
# # Convert the result to a data frame
# total_counts_df <- data.frame(Sample = names(total_counts), TotalCount = total_counts)
# 
# # Create a bar plot using ggplot2
# plot_total_counts <- ggplot(total_counts_df, aes(x = Sample, y = TotalCount)) +
#   geom_bar(stat = "identity") +
#   labs(title = "Total Gene Counts per Sample", x = "Sample", y = "Total Count")
# 
# ggsave("./plot_total_counts.png", plot_total_counts, width = 15, height = 10)
# 
# colSums(counts_df)
# names(colSums(counts_df))
# 
# 
# 
# ### OR ###
# count_long <- counts_df %>%
#   rownames_to_column("Gene") %>%
#   pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Count")
# 
# 
# plot_counts <- ggplot(count_long, aes(x = Sample, y = Count, color = Sample)) +
#   geom_point() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   labs(x = "Sample", y = "Count", title = "Count Matrix Scatter Plot") +
#   
#   
#   ggsave("./plot_counts_txi.png", plot_total_counts, width = 15, height = 10)
# 
# 
# 
# 
# samples[13]


# Calculate average gene counts per group
# avg_gene_counts <- total_counts_df %>%
#   group_by(condition) %>%
#   summarise(AveGeneCount = mean(TotalCount))
#             
# # Create a box plot
# ggplot(avg_gene_counts, aes(x = condition, y = AveGeneCount)) +
#   geom_boxplot() +
#   labs(title = "Average Total Gene Counts per Sample",
#        x = "Sample Group",
#        y = "Average Total Gene Count") +
#   theme_minimal()
