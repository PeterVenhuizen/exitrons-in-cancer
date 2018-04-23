#!/usr/bin/env Rscript
# https://bioconductor.org/packages/3.7/bioc/vignettes/tximport/inst/doc/tximport.html

suppressMessages(library("optparse"))
suppressMessages(library("tximport"))
suppressMessages(library("readr"))
suppressMessages(library("edgeR"))

option_list = list(
    make_option(c("-s", "--sampleTable"), type="character", default=NULL, help="Sample table csv file."),
    make_option("--tx2gene", type="character", default=NULL, help="Transcript gene association list."),
    make_option(c("-o", "--output"), type="character", default=NULL, help="Output file name."),
    make_option(c("-r", "--reference", type="character", default=NULL, help="Reference condition."))
)
opt_parser = OptionParser(option_list=option_list)
args = parse_args(opt_parser)

sampleTable <- read.csv(args$sampleTable, row.names=1)

# tximport Salmon data
files <- paste0(sampleTable$samplePath, "quant.sf")
names(files) <- row.names(sampleTable)
tx2gene <- read.csv(args$tx2gene, sep="\t")
txi.salmon <- tximport(files, type="salmon", tx2gene = tx2gene)

# edgeR with Salmon
cts <- txi.salmon$counts
normMat <- txi.salmon$length
normMat <- normMat/exp(rowMeans(log(normMat)))
o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))

y <- DGEList(cts, group=sampleTable$genotype)
y <- y[rowSums(1e+6 * y$counts/expandAsMatrix(y$samples$lib.size, dim(y)) > 1) >= 3, ]
y$offset <- t(t(log(normMat)) + o)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)

#et <- exactTest(y, pair=levels(sampleTable$condition))
et <- exactTest(y, pair=c(args$ref, setdiff(sampleTable$genotype, args$ref)))
results <- topTags(et, n=Inf, p.value=0.05)
write.table(results, file=args$output, quote=FALSE, sep="\t")