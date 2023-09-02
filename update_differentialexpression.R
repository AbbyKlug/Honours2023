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

## Yi's tx2gene
tx2gene = readRDS("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/DIFF-EXP/tx2gene.rds")
setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/FASTQ.files/")

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

plotPCA(vsd, intgroup="condition") +
  geom_label(aes(label = name))

setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp_update/")

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

plotPCA(vsd, intgroup="condition") +
  geom_label(aes(label = name))


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

plotPCA(vsd, intgroup="condition") +
  geom_label(aes(label = name))


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
tx2gene = readRDS("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/DIFF-EXP/tx2gene.rds")
setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/FASTQ.files/")
###### AD vs Old ######

files = c('./AD/22-2T-AD.fastq.gz_counts/quant.sf', #3
          './AD/25-5T-AD.fastq.gz_counts/quant.sf', #6
          './AD/26-3A-AD.fastq.gz_counts/quant.sf', #7
          './AD/27-5A-AD.fastq.gz_counts/quant.sf', #8
          './AD/28-8T-AD.fastq.gz_counts/quant.sf', #9
          './Old/10-8A-Old.fastq.gz_counts/quant.sf', #1
          './Old/11-10T-Old.fastq.gz_counts/quant.sf',#2
          './Old/12-6A-Old.fastq.gz_counts/quant.sf', #3
          './Old/13-11T-Old.fastq.gz_counts/quant.sf', #4
          './Old/14-7A-Old.fastq.gz_counts/quant.sf', #5
          './Old/15-13T-Old.fastq.gz_counts/quant.sf', #6
          './Old/16-14T-Old.fastq.gz_counts/quant.sf') #7

txi.salmon <- tximport(files = files, 
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)

sampleTable <- data.frame(condition = factor(c("AD", "AD","AD","AD","AD","Old","Old","Old","Old","Old","Old","Old")))
names <-  c("AD3","AD6","AD7", "AD8","AD9",
             "Old1","Old2","Old3","Old4","Old5","Old6","Old7")
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


###### Old vs Young ######
files = c('./Old/10-8A-Old.fastq.gz_counts/quant.sf', #1
          './Old/11-10T-Old.fastq.gz_counts/quant.sf',#2
          './Old/12-6A-Old.fastq.gz_counts/quant.sf', #3
          './Old/13-11T-Old.fastq.gz_counts/quant.sf', #4
          './Old/14-7A-Old.fastq.gz_counts/quant.sf', #5
          './Old/15-13T-Old.fastq.gz_counts/quant.sf', #6
          './Old/16-14T-Old.fastq.gz_counts/quant.sf', #7
          './Young/2-12A-Young.fastq.gz_counts/quant.sf', #1
          './Young/3-17T-Young.fastq.gz_counts/quant.sf', #2
          './Young/4-13A-Young.fastq.gz_counts/quant.sf', #3
          './Young/6-14A-Young.fastq.gz_counts/quant.sf') #5
txi.salmon <- tximport(files = files, 
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)

sampleTable <- data.frame(condition = factor(c("Old","Old","Old","Old","Old","Old","Old", "Young","Young","Young","Young")))
names <-  c("Old1","Old2","Old3","Old4","Old5","Old6","Old7", 
             "Young1","Young2","Young3","Young5")
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


###### AD vs Young ######

files = c('./AD/22-2T-AD.fastq.gz_counts/quant.sf', #3
          './AD/25-5T-AD.fastq.gz_counts/quant.sf', #6
          './AD/26-3A-AD.fastq.gz_counts/quant.sf', #7
          './AD/27-5A-AD.fastq.gz_counts/quant.sf', #8
          './AD/28-8T-AD.fastq.gz_counts/quant.sf', #9
          './Young/2-12A-Young.fastq.gz_counts/quant.sf', #1
          './Young/3-17T-Young.fastq.gz_counts/quant.sf', #2
          './Young/4-13A-Young.fastq.gz_counts/quant.sf', #3
          './Young/6-14A-Young.fastq.gz_counts/quant.sf') #5

txi.salmon <- tximport(files = files, 
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)


sampleTable <- data.frame(condition = factor(c("AD", "AD","AD","AD","AD", "Young","Young","Young","Young")))
names <-  c("AD3","AD6","AD7", "AD8","AD9", 
            "Young1","Young2","Young3","Young5")
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



####### 2nd round selection: DGE SELECTING SAMPLES PER CONDITION #######

# files = c('./AD/22-2T-AD.fastq.gz_counts/quant.sf', #3
#           './AD/25-5T-AD.fastq.gz_counts/quant.sf', #6
#           './AD/26-3A-AD.fastq.gz_counts/quant.sf', #7
#           './AD/28-8T-AD.fastq.gz_counts/quant.sf', #9
#           './Old/13-11T-Old.fastq.gz_counts/quant.sf', #4
#           './Old/14-7A-Old.fastq.gz_counts/quant.sf', #5
#           './Old/15-13T-Old.fastq.gz_counts/quant.sf', #6
#           './Old/16-14T-Old.fastq.gz_counts/quant.sf', #7
#           './Young/2-12A-Young.fastq.gz_counts/quant.sf', #1
#           './Young/3-17T-Young.fastq.gz_counts/quant.sf', #2
#           './Young/4-13A-Young.fastq.gz_counts/quant.sf', #3
#           './Young/6-14A-Young.fastq.gz_counts/quant.sf') #5

## Yi's tx2gene

tx2gene = readRDS("D://ABBY.windows_surface_pro/DIFF-EXP/tx2gene.rds")
setwd("D://ABBY.windows_surface_pro/FASTQ.files/")

###### AD vs Old ######

files = c('./AD/22-2T-AD.fastq.gz_counts/quant.sf', #3
          './AD/25-5T-AD.fastq.gz_counts/quant.sf', #6
          './AD/26-3A-AD.fastq.gz_counts/quant.sf', #7
          './AD/28-8T-AD.fastq.gz_counts/quant.sf', #9
          './Old/13-11T-Old.fastq.gz_counts/quant.sf', #4
          './Old/14-7A-Old.fastq.gz_counts/quant.sf', #5
          './Old/15-13T-Old.fastq.gz_counts/quant.sf', #6
          './Old/16-14T-Old.fastq.gz_counts/quant.sf') #7

txi.salmon <- tximport(files = files, 
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)

sampleTable <- data.frame(condition = factor(c("AD", "AD","AD","AD","Old","Old","Old","Old")))
names <-  c("AD3","AD6","AD7","AD9",
            "Old4","Old5","Old6","Old7")
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


###### Old vs Young ######
files = c('./Old/13-11T-Old.fastq.gz_counts/quant.sf', #4
          './Old/14-7A-Old.fastq.gz_counts/quant.sf', #5
          './Old/15-13T-Old.fastq.gz_counts/quant.sf', #6
          './Old/16-14T-Old.fastq.gz_counts/quant.sf', #7
          './Young/2-12A-Young.fastq.gz_counts/quant.sf', #1
          './Young/3-17T-Young.fastq.gz_counts/quant.sf', #2
          './Young/4-13A-Young.fastq.gz_counts/quant.sf', #3
          './Young/6-14A-Young.fastq.gz_counts/quant.sf') #5
txi.salmon <- tximport(files = files, 
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)

sampleTable <- data.frame(condition = factor(c("Old","Old","Old","Old", "Young","Young","Young","Young")))
names <-  c("Old4","Old5","Old6","Old7", 
            "Young1","Young2","Young3","Young5")
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


###### AD vs Young ######

files = c('./AD/22-2T-AD.fastq.gz_counts/quant.sf', #3
          './AD/25-5T-AD.fastq.gz_counts/quant.sf', #6
          './AD/26-3A-AD.fastq.gz_counts/quant.sf', #7
          './AD/28-8T-AD.fastq.gz_counts/quant.sf', #9
          './Young/2-12A-Young.fastq.gz_counts/quant.sf', #1
          './Young/3-17T-Young.fastq.gz_counts/quant.sf', #2
          './Young/4-13A-Young.fastq.gz_counts/quant.sf', #3
          './Young/6-14A-Young.fastq.gz_counts/quant.sf') #5

txi.salmon <- tximport(files = files, 
                       type = "salmon",
                       txOut = FALSE,
                       tx2gene = tx2gene,
                       ignoreAfterBar = T)


sampleTable <- data.frame(condition = factor(c("AD","AD","AD","AD", "Young","Young","Young","Young")))
names <-  c("AD3","AD6","AD7", "AD9", 
            "Young1","Young2","Young3","Young5")
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



#### 3 per condition ####


## Yi's tx2gene

tx2gene = readRDS("D://ABBY.windows_surface_pro/DIFF-EXP/tx2gene.rds")
setwd("D://ABBY.windows_surface_pro/FASTQ.files/")

###### AD vs Old ######

files = c('./AD/22-2T-AD.fastq.gz_counts/quant.sf', #3
          './AD/26-3A-AD.fastq.gz_counts/quant.sf', #7
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
names <-  c("AD3","AD7","AD9",
            "Old4","Old5","Old6")
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

## PCA ##

vsd <- vst(dds_OldvsYoung, blind = FALSE) # or varianceStabilizingTransformation instead of vst

plotPCA(vsd, intgroup="condition") +
  geom_label(aes(label = name))


###### AD vs Young ######

files = c('./AD/22-2T-AD.fastq.gz_counts/quant.sf', #3
          './AD/26-3A-AD.fastq.gz_counts/quant.sf', #7
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
names <-  c("AD3","AD7", "AD9", 
            "Young1","Young3","Young5")
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


#### AD 3,6,9 

#### 3 per condition ####


## Yi's tx2gene

tx2gene = readRDS("D://ABBY.windows_surface_pro/DIFF-EXP/tx2gene.rds")
setwd("D://ABBY.windows_surface_pro/FASTQ.files/")

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

plotPCA(vsd, intgroup="condition") +
  geom_label(aes(label = name))


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

plotPCA(vsd, intgroup="condition") +
  geom_label(aes(label = name))


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

plotPCA(vsd, intgroup="condition") +
  geom_label(aes(label = name))



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
#filtered_OldvsYoung_basemean$Gene = gsub("(ENSG[0-9]+)\\.[0-9]+","\\1",filtered_OldvsYoung_basemean$Gene)
OldvsYoung_genes <- OldvsYoung$Gene

ADvsYoung$Gene = gsub("(ENSG[0-9]+)\\.[0-9]+","\\1",ADvsYoung$Gene)
#filtered_ADvsYoung_basemean$Gene = gsub("(ENSG[0-9]+)\\.[0-9]+","\\1",filtered_ADvsYoung_basemean$Gene)
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




### PLOTTING CLUSTERPROFILER ###
setwd("~/../../media/intel6700/Seagate Expansion Drive/ABBY.windows_surface_pro/differentialexp_update/")
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
  geom_bar(stat = "identity", fill = "hotpink") +  
  labs(title = "AD-specific Significant Enriched Terms (UP)", x = "Count", y = "Biological Process")
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
  labs(title = "AD-specific Significant Enriched Terms (DOWN)", x = "Count", y = "Biological Process")

plot_enrich_ADOld_down <- edit_plots(plot_enrich_ADOld_down)
plot_enrich_ADOld_down











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
