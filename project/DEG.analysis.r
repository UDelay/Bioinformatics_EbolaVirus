#project_code_mine
### RNA-seq analysis ###
require(DESeq2)
require(ggplot2)

#after install -> 3 lines #
install.packages("ggplot2")
require(BiocManager)
BiocManager::install("DESeq2")

ensSymbMapFile = "ensembl.symbol.map.txt"
ensSymbMap = read.table(ensSymbMapFile, header=T, stringsAsFactors = F)

#rcFile = "rcMat.RData"
#rcMetadataFile  =  "rcMetadata.RData"

#load(rcFile)
#load(rcMetadataFile)

restv_s1 = read.table("2Drestv_s1.tsv", sep="\t")
restv_s2 = read.table("2Drestv_s2.tsv", sep="\t")
restv_s3 = read.table("2Drestv_s3.tsv", sep="\t")
mock_s1 = read.table("2Dmock_s1.tsv", sep="\t")
mock_s2 = read.table("2Dmock_s2.tsv", sep="\t")
mock_s3 = read.table("2Dmock_s3.tsv", sep="\t")
ebov_s1 = read.table("2Debov_s1.tsv", sep="\t")
ebov_s2 = read.table("2Debov_s2.tsv", sep="\t")
ebov_s3 = read.table("2Debov_s3.tsv", sep="\t")


rcMatR_ = cbind(mock_s1[2:dim(mock_s1)[1],4], mock_s2[2:dim(mock_s1)[1],4])
rcMatR_ = cbind(rcMatR_, mock_s3[2:dim(mock_s1)[1],4])
rcMatR_ = cbind(rcMatR_, restv_s1[2:dim(mock_s1)[1],4])
rcMatR_ = cbind(rcMatR_, restv_s2[2:dim(mock_s1)[1],4])
rcMatR_ = cbind(rcMatR_, restv_s3[2:dim(mock_s1)[1],4])

rcMatE_ = cbind(mock_s1[2:dim(mock_s1)[1],4], mock_s2[2:dim(mock_s1)[1],4])
rcMatE_ = cbind(rcMatE_, mock_s3[2:dim(mock_s1)[1],4])
rcMatE_ = cbind(rcMatE_, ebov_s1[2:dim(mock_s1)[1],4])
rcMatE_ = cbind(rcMatE_, ebov_s2[2:dim(mock_s1)[1],4])
rcMatE_ = cbind(rcMatE_, ebov_s3[2:dim(mock_s1)[1],4])


rcMatR0 = matrix(0,dim(rcMatR_)[1],dim(rcMatR_)[2])
for(i in 1:dim(rcMatR0)[1]){
  rcMatR0[i,] = round(as.numeric(rcMatR_[i,]))
}
rcMatE0 = matrix(0,dim(rcMatE_)[1],dim(rcMatE_)[2])
for(i in 1:dim(rcMatE0)[1]){
  rcMatE0[i,] = round(as.numeric(rcMatE_[i,]))
}
# rcMat0 may contain a row [0, 0, ... , 0]


rcMatR0_m = rcMatR0[,1:3]
rcMatR0_v = rcMatR0[,4:6]
row_sumR_m = apply(rcMatR0_m, 1, sum)
row_sumR_v = apply(rcMatR0_v, 1, sum)

rcMatE0_m = rcMatE0[,1:3]
rcMatE0_v = rcMatE0[,4:6]
row_sumE_m = apply(rcMatE0_m, 1, sum)
row_sumE_v = apply(rcMatE0_v, 1, sum)
  
rcMatR = c()
rcMatR_rownames = c()
count = 0
for(i in 1:dim(rcMatR_)[1]){
  if ((rcMatR0_m[i]>0) & (rcMatR0_v[i]>0)){
    rcMatR = rbind(rcMatR, rcMatR0[i,])
    count = count + 1
    rcMatR_rownames[count] = restv_s1[1+i,1]
  }
}

rcMatE = c()
rcMatE_rownames = c()
count = 0
for(i in 1:dim(rcMatE_)[1]){
  if ((rcMatE0_m[i]>0) & (rcMatE0_v[i]>0)){
    rcMatE = rbind(rcMatE, rcMatE0[i,])
    count = count + 1
    rcMatE_rownames[count] = restv_s1[1+i,1]
  }
}

rownames(rcMatR) = rcMatR_rownames
colnames(rcMatR) = c("x001","x002","x003","x004","x005","x006")
rownames(rcMatE) = rcMatE_rownames
colnames(rcMatE) = c("x001","x002","x003","x004","x005","x006")
rcMetadata = data.frame(Sample = c('x001','x002',"x003","x004","x005","x006"),Condition = c('Control','Control','Control','Treated','Treated','Treated'))




dds.rcMatR <- DESeqDataSetFromMatrix(countData = rcMatR,
                                     colData = rcMetadata,
                                     design = ~ Condition)
dds.rcMatR <- estimateSizeFactors(dds.rcMatR)

dds.rcMatR.normalized = DESeq2::counts(dds.rcMatR, normalized = T)
write.table(dds.rcMatR.normalized, "dds.rcMatR.normalized.coun.txt", sep="\t")

dds.rcMatR <- DESeq(dds.rcMatR)
resultsNames(dds.rcMatR)

dds.rcMatR.lfc <- lfcShrink(dds.rcMatR, coef="Condition_Treated_vs_Control", type="apeglm")
dds.rcMatR.lfc <- as.data.frame(dds.rcMatR.lfc)
dds.rcMatR.lfc$gene = rownames(dds.rcMatR.lfc)
dds.rcMatR.lfc$symbol = ensSymbMap$Symbol[match(dds.rcMatR.lfc$gene, ensSymbMap$Ensembl)]
dds.rcMatR.lfc = dds.rcMatR.lfc[dds.rcMatR.lfc$baseMean != 0, ]
dds.rcMatR.lfc = dds.rcMatR.lfc[!is.na(dds.rcMatR.lfc$padj), ]
dds.rcMatR.lfc$class = "lightgray"
dds.rcMatR.lfc$class[dds.rcMatR.lfc$padj < 0.01] = "gray"
dds.rcMatR.lfc$class[dds.rcMatR.lfc$log2FoldChange > 1 & dds.rcMatR.lfc$padj < 0.01] = "maroon"
dds.rcMatR.lfc$class[dds.rcMatR.lfc$log2FoldChange < -1 & dds.rcMatR.lfc$padj < 0.01] = "skyblue"

write.table(dds.rcMatR.lfc, "dds.rcMatR.lfc.stats.txt", sep="\t")




dds.rcMatE <- DESeqDataSetFromMatrix(countData = rcMatE,
                                     colData = rcMetadata,
                                     design = ~ Condition)
dds.rcMatE <- estimateSizeFactors(dds.rcMatE)

dds.rcMatE.normalized = DESeq2::counts(dds.rcMatE, normalized = T)
write.table(dds.rcMatE.normalized, "dds.rcMatE.normalized.coun.txt", sep="\t")

dds.rcMatE <- DESeq(dds.rcMatE)
resultsNames(dds.rcMatE)

dds.rcMatE.lfc <- lfcShrink(dds.rcMatE, coef="Condition_Treated_vs_Control", type="apeglm")
dds.rcMatE.lfc <- as.data.frame(dds.rcMatE.lfc)
dds.rcMatE.lfc$gene = rownames(dds.rcMatE.lfc)
dds.rcMatE.lfc$symbol = ensSymbMap$Symbol[match(dds.rcMatE.lfc$gene, ensSymbMap$Ensembl)]
dds.rcMatE.lfc = dds.rcMatE.lfc[dds.rcMatE.lfc$baseMean != 0, ]
dds.rcMatE.lfc = dds.rcMatE.lfc[!is.na(dds.rcMatE.lfc$padj), ]
dds.rcMatE.lfc$class = "lightgray"
dds.rcMatE.lfc$class[dds.rcMatE.lfc$padj < 0.01] = "gray"
dds.rcMatE.lfc$class[dds.rcMatE.lfc$log2FoldChange > 1 & dds.rcMatE.lfc$padj < 0.01] = "maroon"
dds.rcMatE.lfc$class[dds.rcMatE.lfc$log2FoldChange < -1 & dds.rcMatE.lfc$padj < 0.01] = "skyblue"

write.table(dds.rcMatE.lfc, "dds.rcMatE.lfc.stats.txt", sep="\t")


# ggplot(dds.rcMatR.lfc, aes(x=log2FoldChange, y=-log10(padj), colour = class)) +
#   geom_point() + xlim(c(-10,10)) + ylim(c(0,70)) + geom_vline(xintercept = c(-1,1)) + geom_hline(yintercept = 2) +
#   scale_color_manual(values = c("lightgray", "gray", "maroon", "skyblue")) +
#   xlab("log2 fold change") + ylab("-log adj. p-value")+
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         panel.border = element_rect(colour = "black", fill=NA, size=.5),
#         legend.position = "none")

# ggplot(dds.rcMatE.lfc, aes(x=log2FoldChange, y=-log10(padj), colour = class)) +
#   geom_point() + xlim(c(-10,10)) + ylim(c(0,70)) + geom_vline(xintercept = c(-1,1)) + geom_hline(yintercept = 2) +
#   scale_color_manual(values = c("lightgray", "gray", "maroon", "skyblue")) +
#   xlab("log2 fold change") + ylab("-log adj. p-value")+
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         panel.border = element_rect(colour = "black", fill=NA, size=.5),
#         legend.position = "none")



