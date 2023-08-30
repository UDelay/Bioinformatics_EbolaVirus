### RNA-seq analysis ###
require(DESeq2)
require(ggplot2)

ensSymbMapFile = "ensembl.symbol.map.txt"
ensSymbMap = read.table(ensSymbMapFile, header=T, stringsAsFactors = F)

rcFile = "rcMat.RData" ##integer로바꿔줌 
rcMetadataFile  =  "rcMetadata.RData"

load(rcFile)
load(rcMetadataFile)

rcMat = round(rcMat)


dds.rcMat <- DESeqDataSetFromMatrix(countData = rcMat,
                                     colData = rcMetadata,
                                     design = ~ Condition) #-1넣는 경우도 있음 
dds.rcMat <- estimateSizeFactors(dds.rcMat)

dds.rcMat.normalized = DESeq2::counts(dds.rcMat, normalized = T)
write.table(dds.rcMat.normalized, "dds.rcMat.normalized.coun.txt", sep="\t")

dds.rcMat <- DESeq(dds.rcMat)
resultsNames(dds.rcMat)

dds.rcMat.lfc <- lfcShrink(dds.rcMat, coef="Condition_Treated_vs_Control", type="apeglm")
dds.rcMat.lfc <- as.data.frame(dds.rcMat.lfc)
dds.rcMat.lfc$gene = rownames(dds.rcMat.lfc)
dds.rcMat.lfc$symbol = ensSymbMap$Symbol[match(dds.rcMat.lfc$gene, ensSymbMap$Ensembl)]
dds.rcMat.lfc = dds.rcMat.lfc[dds.rcMat.lfc$baseMean != 0, ]
dds.rcMat.lfc = dds.rcMat.lfc[!is.na(dds.rcMat.lfc$padj), ]
dds.rcMat.lfc$class = "lightgray"
dds.rcMat.lfc$class[dds.rcMat.lfc$padj < 0.01] = "gray"
dds.rcMat.lfc$class[dds.rcMat.lfc$log2FoldChange > 1 & dds.rcMat.lfc$padj < 0.01] = "maroon"
dds.rcMat.lfc$class[dds.rcMat.lfc$log2FoldChange < -1 & dds.rcMat.lfc$padj < 0.01] = "skyblue"

write.table(dds.rcMat.lfc, "dds.rcMat.lfc.stats.txt", sep="\t")


ggplot(dds.rcMat.lfc, aes(x=log2FoldChange, y=-log10(padj), colour = class)) + 
  geom_point() + xlim(c(-5,5)) +geom_vline(xintercept = c(-1,1)) + geom_hline(yintercept = 2) + 
  scale_color_manual(values = c("lightgray", "gray", "maroon", "skyblue")) +
  xlab("log2 fold change") + ylab("-log adj. p-value")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),  
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.5),
        legend.position = "none")



