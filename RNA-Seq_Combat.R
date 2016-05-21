require(edgeR)
require(sva)
source('code/annotate_edgeR.R')
files = data.frame(files=c('counts.control0', 'counts.control1',
'counts.control2', 'counts.control3', 'counts.treatment0',
'counts.treatment1', 'counts.treatment2', 'counts.treatment3'),
                   group=c('control', 'control', 'control', 'control',
'treatment', 'treatment', 'treatment', 'treatment'),
                   day=rep(0:3,2)
)
labels <- paste0(files$group, files$day)
dge <- readDGE(files=files, path='data/HTSeq/', labels=labels)
rownames(dge$counts) <- paste0('Gene', 1:nrow(dge$counts)) #Change gene
names to anonymize data
################################
# save(dge, file='objs/dge.Rdata')
#     SEE ATTACHED FILE     #
###############################

##  filter out the no_feature etc. rows
dge <- dge[1:(nrow(dge)-5), ]
##  This mitochondrial rRNA gene takes up a massive portion of my libraries
dge <- dge[!rownames(dge)%in%'Gene13515', ]
##  filter out lowly expressed genes
keep <- rowSums(cpm(dge) > 1) >= 3 ## gene has at least 3 columns where cpm
is > 1
dge <- dge[keep, ]
##  Recompute library sizes
dge$samples$lib.size <- colSums(dge$counts)
##  Normalize for lib size
dge <- calcNormFactors(dge)

## ComBat
mod <- model.matrix(~as.factor(group), data=dge$sample)
mod0 <- model.matrix(~1, data=dge$sample)
batch <- dge$sample$day

combat <- ComBat(dat=cpm(dge), batch=batch, mod=mod)
pval_combat = f.pvalue(combat, mod, mod0)
padj_combat = p.adjust(pval_combat, method="BH")
mean_control <- rowMeans(combat[, 1:4])
mean_treatment <- rowMeans(combat[, 5:8])
logFC <- log2(mean_treatment/mean_control)

res <- data.frame(mean_control, mean_treatment, logFC, pval=pval_combat,
padj=padj_combat)
res <- res[order(res$padj), ]
