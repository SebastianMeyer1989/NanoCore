library(pheatmap, warn.conflicts = FALSE)

args = commandArgs(trailingOnly=TRUE)

min_cov <- data.matrix(read.table(args[1], header=TRUE, row.names=1, sep="\t",check.names = FALSE))
hete <- data.matrix(read.table(args[2], header=TRUE, row.names=1, sep="\t",check.names = FALSE))
cov <- data.matrix(read.table(args[3], header=TRUE, row.names=1, sep="\t",check.names = FALSE))
mapq <- data.matrix(read.table(args[4], header=TRUE, row.names=1, sep="\t",check.names = FALSE))
stats <- read.table(args[5], header=TRUE, sep="\t",check.names = FALSE)
reads <- read.table(args[6], header=TRUE, sep="\t",check.names = FALSE)


### Boxplot

pdf(args[7])
par(mar=c(5,7.5,1,1))
boxplot(t(min_cov),
        horizontal = TRUE,
        col="dodgerblue3",
        xlab = "Ratio of gene positions above chosen covereage",
        las = 2,
        )
title(ylab = "Samples", line = 6)
dev.off()


### Heatmaps

pdf(args[8], h=50)
pheatmap(t(min_cov), fontsize_row=2, cluster_cols=F, cluster_rows=F)
dev.off()

pdf(args[9], h=20)
pheatmap(t(hete[,colSums(hete!=0)>0]), fontsize_row=2, treeheight_row=0, treeheight_col=0, cluster_cols=F)
dev.off()

pdf(args[10], h=20)
pheatmap(t(cov[,colSums(cov!=0)>0]), fontsize_row=2, treeheight_row=0, treeheight_col=0, cluster_cols=F)
dev.off()

pdf(args[11], h=20)
pheatmap(t(mapq[,colSums(mapq!=0)>0]), fontsize_row=2, treeheight_row=0, treeheight_col=0, cluster_cols=F)
dev.off()

### Basic Stats

pdf(args[12])
par(mfrow=c(2,2))
par(mar=c(9,6,1,1))

barplot(stats$Total_Bases, names.arg=stats$ID, las=2, cex.names=0.6)
title(xlab = "Samples", line = 6)
title(ylab = "Total bases", line = 5)

barplot(stats$Mapped_Bases, names.arg=stats$ID, las=2, cex.names=0.6)
title(xlab = "Samples", line = 6)
title(ylab = "Mapped bases", line = 4)

barplot(stats$Total_Reads, names.arg=stats$ID, las=2, cex.names=0.6)
title(xlab = "Samples", line = 6)
title(ylab = "Total reads", line = 4)

barplot(stats$Mapped_Reads, names.arg=stats$ID, las=2, cex.names=0.6)
title(xlab = "Samples", line = 6)
title(ylab = "Mapped reads", line = 4)

dev.off()

pdf(args[13])
par(mfrow=c(2,2))
par(mar=c(9,6,1,1))

barplot(stats$Average_Coverage, names.arg=stats$ID, las=2, cex.names=0.6)
title(xlab = "Samples", line = 6)
title(ylab = "Average coverage", line = 3)

hist(stats$Average_Coverage, xlab="Average coverage (hist)", main="", breaks=max(stats$Average_Coverage)/30)

hist(reads$Length_of_all_reads, xlab="Read length (hist)", main="")

hist(reads$Length_of_all_reads, xlab="Read length (hist, cutoff)", main="", ylim=c(0,20000))

dev.off()
