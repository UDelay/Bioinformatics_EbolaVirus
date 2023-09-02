require(ggplot2)

#load("DEG_ebov_restv.Rdata")

draw_plot <- function(RorE, xlim, ylim, dotsize) {
  # xlim : default 1
  # ylim : default 2
  
  dds.rcMatR.lfc$class = "lightgray"
  dds.rcMatR.lfc$class[dds.rcMatE.lfc$padj < 10^(-ylim)] = "gray"
  dds.rcMatR.lfc$class[dds.rcMatR.lfc$log2FoldChange > xlim & dds.rcMatR.lfc$padj < 10^(-ylim)] = "maroon"
  dds.rcMatR.lfc$class[dds.rcMatR.lfc$log2FoldChange < -xlim & dds.rcMatR.lfc$padj < 10^(-ylim)] = "skyblue"
  
  dds.rcMatE.lfc$class = "lightgray"
  dds.rcMatE.lfc$class[dds.rcMatE.lfc$padj < 10^(-ylim)] = "gray"
  dds.rcMatE.lfc$class[dds.rcMatE.lfc$log2FoldChange > xlim & dds.rcMatE.lfc$padj < 10^(-ylim)] = "maroon"
  dds.rcMatE.lfc$class[dds.rcMatE.lfc$log2FoldChange < -xlim & dds.rcMatE.lfc$padj < 10^(-ylim)] = "skyblue"
  
  if (RorE == 'R' | RorE == 'r') {
    ggplot(dds.rcMatR.lfc, aes(x=log2FoldChange, y=-log10(padj), colour = class)) +
      geom_point(size = dotsize) + xlim(c(-10,10)) + ylim(c(0,70)) +
      geom_vline(xintercept = c(-xlim,xlim),linetype = 'dashed') + geom_hline(yintercept = ylim,linetype = 'dashed') +
      scale_color_manual(values = c("lightgray", "gray", "maroon", "skyblue")) +
      xlab("log2 fold change") + ylab("-log adj. p-value")+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=.5),
            legend.position = "none")
  } else if (RorE == 'E' | RorE == 'e') {
    ggplot(dds.rcMatE.lfc, aes(x=log2FoldChange, y=-log10(padj), colour = class)) +
      geom_point(size = dotsize) + xlim(c(-10,10)) + ylim(c(0,70)) + 
      geom_vline(xintercept = c(-xlim,xlim),linetype = 'dashed') + geom_hline(yintercept = ylim,linetype = 'dashed') +
      scale_color_manual(values = c("lightgray", "gray", "maroon", "skyblue")) +
      xlab("log2 fold change") + ylab("-log adj. p-value")+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=.5),
            legend.position = "none")
  } else {
    print("Wrong input.")
  }
}




