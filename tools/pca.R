library("ggplot2")
library(DESeq2)
library("RColorBrewer")
library("gplots")
library(reshape2)

args <- commandArgs(TRUE)
csvFilePath = args[1]
pcaOutputFpath = args[2]
geneCombinedCounts = args[3]

meta <- read.csv(csvFilePath, stringsAsFactors=F, header=T)
rownames(meta) <- meta$description
if(!is.element("condition", colnames(meta))){{
  meta <- data.frame(meta, condition1=rep("group1", nrow(meta)), condition2=meta$description)
}} else {{
  meta <- data.frame(meta, condition1=meta$condition, condition2=meta$condition)
}}

metax <- data.frame(name=as.character(meta$description), condition=as.character(meta$condition2))
rownames(metax) <-  metax$name

countsx <- read.delim(geneCombinedCounts, header=F, stringsAsFactors = F)
rownames(countsx) <- countsx[,1]
colnames(countsx) <- countsx[1,]
countsx <- countsx[-1,-1]
countsx <- countsx[, rownames(metax)]
counts <- apply(countsx, 2, as.numeric)
rownames(counts) <- rownames(countsx)

dds <- DESeqDataSetFromMatrix(countData = counts, colData=metax, design= ~ condition)
rld <- rlog(dds)

data <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE, ntop=1000)
percentVar <- round(100 * attr(data, "percentVar"))
write.table(data, file=pcaOutputFpath, sep="\t")
line = paste0("#Variance:", percentVar[1], ",", percentVar[2])
write(line, file=pcaOutputFpath, append=TRUE)