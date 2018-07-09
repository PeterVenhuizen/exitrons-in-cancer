library(tximport)
library(DESeq2)

# http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#salmon-sailfish-with-inferential-replicates
# https://master.bioconductor.org/help/course-materials/2016/CSAMA/lab-3-rnaseq/rnaseq_gene_CSAMA2016.html#pca-plot

sampleTable <- read.csv("sample_table.csv", row.names=1)
files <- paste0(sampleTable$samplePath, "quant.sf")
names(files) <- row.names(sampleTable)

if (all(file.exists(files))) {
    tx2gene <- read.csv("/home/venhuip8/EI_in_cancer/hg38/tx2gene.csv", sep="\t")
    txi <- tximport(files, type="salmon", tx2gene=tx2gene)

    sampleTable <- data.frame(genotype = factor(sampleTable$genotype), pairs = factor(sampleTable$pairs))
    rownames(sampleTable) <- colnames(txi$counts)
    
    dds <- DESeqDataSetFromTximport(txi, sampleTable, ~genotype)
    vsd <- vst(dds)
    
    pcaData <- plotPCA(vsd, intgroup=c("genotype", "pairs"), returnData=TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    ggplot(pcaData, aes(PC1, PC2, color=genotype)) +
    #ggplot(pcaData, aes(PC1, PC2, color=genotype, label=pairs)) +
        scale_shape_manual(values=sprintf("%s", seq(1:length(sampleTable$pairs)/2))) +
        geom_point(size=3) +
        #geom_text() +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        coord_fixed()
    
}
