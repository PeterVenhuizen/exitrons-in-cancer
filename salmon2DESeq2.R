#!/usr/bin/env Rscript
# http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#salmon-sailfish-with-inferential-replicates
# https://master.bioconductor.org/help/course-materials/2016/CSAMA/lab-3-rnaseq/rnaseq_gene_CSAMA2016.html#pca-plot

suppressMessages(library("optparse"))
suppressMessages(library("tximport"))
suppressMessages(library("readr"))
suppressMessages(library("DESeq2"))

option_list = list(
    make_option(c("-s", "--sampleTable"), type="character", default=NULL, help="Sample table csv file."),
    make_option("--tx2gene", type="character", default=NULL, help="Transcript gene association list."),
    make_option(c("-r", "--reference"), type="character", default=NULL, help="Name of reference condition"),
    make_option(c("-o", "--output"), type="character", default=NULL, help="Output file name.")
)
opt_parser = OptionParser(option_list=option_list)
args = parse_args(opt_parser)

sampleTable <- read.csv(args$sampleTable, row.names=1)
files <- paste0(sampleTable$samplePath, "quant.sf")
names(files) <- row.names(sampleTable)

if (all(file.exists(files))) {
    tx2gene <- read.csv("/home/venhuip8/EI_in_cancer/hg38/tx2gene.csv", sep="\t")
    txi <- tximport(files, type="salmon", tx2gene=tx2gene)
    
    samples <- data.frame(genotype = factor(sampleTable$genotype), pairs = factor(sampleTable$pairs))
    rownames(sampleTable) <- colnames(txi$counts)
    
    dds <- DESeqDataSetFromTximport(txi, samples, ~genotype)
    design(dds) <- formula(~ pairs + genotype)
    dds <- DESeq(dds)

    res <- results(dds)
    resOrdered <- res[order(res$padj),]
    #res10 <- subset(resOrdered, padj < 0.1) # Everything with an adjusted p value < 0.1 (default)
    #write.table(res10, file=args$output, quote=FALSE, sep="\t")
    write.table(resOrdered, file=args$output, quote=FALSE, sep="\t")
    
}
