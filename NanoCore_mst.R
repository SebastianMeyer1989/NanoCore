library(ape, warn.conflicts = FALSE)
library(pegas, warn.conflicts = FALSE)
library(sna, warn.conflicts = FALSE)

args = commandArgs(trailingOnly=TRUE)
alleles <- read.table(args[1], header=TRUE, row.names=1, sep="\t",check.names = FALSE)

distanceMatrix <- as.matrix(alleles, matrix.type = c("adjacency"))

plot_indices <- c(1:length(rownames(distanceMatrix)))
myMST <- mst(distanceMatrix)
myMST_matrix <- matrix(c(0), ncol = length(plot_indices), nrow = length(plot_indices))
rownames(myMST_matrix) <- rownames(distanceMatrix)
colnames(myMST_matrix) <- rownames(distanceMatrix)
labels_myMST <- attr(myMST, "labels")
for(linkI in 1:dim(myMST)[[1]])
{
  label1 <- labels_myMST[[myMST[linkI,][[1]]]]
  label2 <- labels_myMST[[myMST[linkI,][[2]]]]
  myMST_matrix[c(label1),c(label2)] <- 1
  myMST_matrix[c(label2),c(label1)] <- 1
}

plotcord <- data.frame(gplot.layout.fruchtermanreingold(myMST_matrix, NULL))
colnames(plotcord) = c("X1","X2")
rownames(plotcord) = rownames(distanceMatrix)

pdf(args[2])

plot(plotcord, bg = "white", col = "white", cex = 0, axes = F, xlab = "", ylab = "", main = "")
coord <- par("usr")
  
eList <- NULL
for ( i in 1:nrow(myMST_matrix) ){
  for ( j in 1:ncol(myMST_matrix)) {
    if (!is.na(myMST_matrix[i,j])){
      if (myMST_matrix[i,j]>0){
        eList <- rbind(eList,c( rownames(myMST_matrix)[i], colnames(myMST_matrix)[j], myMST_matrix[i,j]))
      }
    }
  }
}

sum_distances <- 0
edges <- data.frame(plotcord[eList[,1],1:2], plotcord[eList[,2],1:2], eList[,1], eList[,2])
for(eI in 1:(dim(edges)[[1]]))
{
    n1 <- as.character(edges[eI, 5])
    n2 <- as.character(edges[eI, 6])
    distance <- distanceMatrix[[n1, n2]]
    
    lty <- "solid"
    if(distance == 0)
    {
      lty <- "dashed"
    }
    lines(c(edges[eI, 1], edges[eI, 3]), c(edges[eI, 2], edges[eI, 4]), lty = lty, lwd = 0.8) # edges
    
    print(c(n1, n2, distance))
    midPoint_x <- mean(c(edges[eI, 1], edges[eI, 3]))
    midPoint_y <- mean(c(edges[eI, 2], edges[eI, 4]))
    if(distance >= 1)
    {
      text(midPoint_x, midPoint_y, labels = distance, adj = 0.5, cex = 0.8) # edge labels
    }
    sum_distances <- sum_distances + distance
}
  
points(plotcord, pch = 21, bg = "gray", cex = 3, xlab = "", ylab = "", main = "")	# nodes


for(i in 1:length(rownames(plotcord)))	
{
  colForText <- "black"
  text(plotcord[[1]][[i]], plotcord[[2]][[i]], labels = rownames(plotcord)[[i]], adj = 0.5, cex = 0.6, col = colForText, xpd = TRUE) # node labels
}

dev.off()
