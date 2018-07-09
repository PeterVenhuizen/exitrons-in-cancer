library("ggpubr")

setwd("/home/venhuip8/EI_in_cancer/TCGA-KICH/features")

data <- read.table("exitron_length_gc.txt", header=T, row.names=1, sep="\t")
data$ORIGIN <- factor(data$ORIGIN, levels=c("Normal", "Tumour", "Normal & Tumour", "gencode"))

my_comparisons <- list( c("Normal", "Tumour"), c("Normal", "Normal & Tumour"), c("Tumour", "Normal & Tumour") )

### LENGTH ###
len <- subset(data, ORIGIN != "gencode")
ggboxplot(len, x = "ORIGIN", y = "LENGTH", title = "Exitron length (nt)") +
    scale_y_continuous(trans="sqrt") +
    stat_compare_means(comparisons = my_comparisons) + 
    xlab("") +
    ylab("Length (nt)")

### GC ###
gc <- subset(data, ORIGIN != "gencode")
ggboxplot(gc, x = "ORIGIN", y = "GC", title = "Exitron GC content") + 
    stat_compare_means(comparisons = my_comparisons) +
    xlab("")
