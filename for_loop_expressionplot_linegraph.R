for (gene in gene_list$V1) {
gene_expr <- getGenes(cuff, gene)
p1 <- expressionPlot(gene_expr, showErrorbars = F) + scale_x_discrete(breaks=c("E11_5", "E15_5", "PN"), labels=c("E11.5", "E15.5", "PN")) + xlab("Developmental Stage") + ylab("Expression value (FPKM)") + expand_limits(y=0)
file <- sprintf("%s.linegraph.pdf",gene)
pdf(file)
print(p1)
dev.off()
}
