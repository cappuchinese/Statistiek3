##############################################################################
#
# This file contains several functions for clustering and dendrograms.
#
# Emile Apol, feb-april 2013 
#
# ----------------------------------------------------------------------------
#
# 1.  makeColDendroK          Color labels of k subclusters of a dendrogram
# 2.  makeKMeansDendrogram    Makes a "dendrogram" for the output of kmeans clustering
# 3.  traditionalScreeplot    Plots a "traditional" screeplot, 
#                             based on pca object with class "prcomp" or class "princomp".
# 4.  traditionalPCPlot       Plots a "traditional" PC1-PC2-plot, i.e., individual cases 
#                             projected in the plane of the first 2 eigenvectors PC1 and PC2, 
#                             based on pca object with class "prcomp" or class "princomp".
#
##############################################################################



#######################################################################
#
# 1. function makeColDendroK
#
#
# Emile Apol, based on Gaston Sanchez, 2012
# 14 march 2013
#
# Color labels of k subclusters of a dendrogram 
# x = hclust object or a dendrogram object
#
# Returns: dendrogram object
#
#######################################################################

makeColDendroK <- function(x, k){
  
  # vector of colors
  # labelColors <- rainbow(k)
  labelColors <- palette("default")
  # cut dendrogram in k clusters
  clusMember <- cutree(x, k)
  # if necessary, change class of x into dendrogram
  xd <- as.dendrogram(x)
  # function to get color labels
  colLab <- function(n) {
    if (is.leaf(n)) {
      a <- attributes(n)
      labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
      attr(n, "nodePar") <- c(lab.col = labCol, pch="")
    }
    n
  }
  # using dendrapply
  clusDendro <- dendrapply(xd, colLab)
  return(clusDendro)
}

#######################################################################
# 2. function makeKMeansDendrogram
#
# Emile Apol
# 22 march 2013
#
# Makes a "dendrogram" for the output of kmeans clustering.
# This dendrogram object can be visualized via plot().
#
# Usage:  makeKMeansDendrogram(kcl)
#
# kcl       = cluster object from kmeans function
#
#######################################################################

makeKMeansDendrogram <- function(kcl, dbg=F){
  clusters <- kcl$cluster
  size.clusters <- kcl$size
  n.genes <- length(clusters) # so there will be n.genes-1 hierarchical clusters!
  gene.names <- names(clusters)
  n.clusters <- length(size.clusters)
  if(dbg) cat("Gene names are: ",gene.names,"\n")
  if(dbg) cat("size.clusters = ",size.clusters,"\n")
  NONSINGLE.g <- which(size.clusters > 1)  # which cluster nrs are non-singleton?
  SINGLE.g <- which(size.clusters == 1)    # which cluster nrs are singleton?
  if(dbg) cat("Non-singleton cluster nrs are ",NONSINGLE.g,"\n")
  if(dbg) cat("    Singleton cluster nrs are ",SINGLE.g,"\n")
  n.clusters.NONSINGLE <- length(NONSINGLE.g)
  n.clusters.SINGLE <- length(SINGLE.g)
  if(dbg) cat("Nr of non-singleton clusters = ",n.clusters.NONSINGLE,"\n")
  if(dbg) cat("Nr of     singleton clusters = ",n.clusters.SINGLE,"\n")
  ORDER <- order(clusters)
  if(dbg) cat("ORDER = ",ORDER,"\n")
  gene.names.ord <- gene.names[ORDER]
  
  # Initialisation
  a <- list()         # make empty return list
  mergeMat <- matrix(rep(0, 2*(n.genes-1)), ncol=2) # make merge matrix
  height <- c()       # make empty height vector
  leaf.labels <- c()  # keep track of the labels per new leaf
  singletons.g <- c() # keep track of group index g of the singleton clusters
  singletons.i <- c() # keep track of array index i of the singleton clusters
  singletons.n <- c() # keep track of the names of singleton clusters
  stemRow <- c()      # keep track of the row in mergeMat that contains the stem of the group
  leaf.height <- 1.0  # default height of each leaf in dendrogram
  stem.height <- 2.0  # default height to connect the groups in dendrogram
  rowIndex <- 1
  leafIndex <- 1
  itemIndex <- 1
  # calculate the row nrs in mergeMat to connect the stems of the non-singleton clusters
  rowStem <- c()
  accumulate <- 0
  for(i in 1 : n.clusters.NONSINGLE){
    thisSize <- size.clusters[NONSINGLE.g[i]]
    # rowStem[i] <- accumulate + round(thisSize/2) # this does not seem to work...
    rowStem[i] <- accumulate + thisSize - 1 # this is the old way...
    accumulate <- accumulate + (thisSize - 1)
  }
  if(dbg) cat("rowStem = ",rowStem,"\n")
  
  # STEP I. Make the leaf clusters
  if(dbg) print("Step I mergeMat")
  for (g in 1 : n.clusters){
    newCluster <- T
    if(dbg) cat("g (1) = ",g,"\n")
    if(dbg) cat("j from 1 to ",size.clusters[g]-1,"\n")
    
    if(size.clusters[g] == 1){ 
      # this is a singleton cluster!
      singletons.g <- c(singletons.g, g)
      singletons.i <- c(singletons.i, itemIndex)
      singletons.n <- c(singletons.n, gene.names.ord[itemIndex])
      itemIndex <- itemIndex + 1
      
    } else { 
      # this is not a singleton cluster!
      for(j in 1 : (size.clusters[g]-1)){
        if(dbg) cat("j (1) = ",j,"\n")
        if(newCluster){
          if(dbg) cat("This is a NEW cluster (",g,")!\n")
          mergeMat[rowIndex, 1] <- -(leafIndex) # add new leaf
          leaf.labels <- c(leaf.labels, gene.names.ord[itemIndex])
          leafIndex <- leafIndex + 1
          itemIndex <- itemIndex + 1
          mergeMat[rowIndex, 2] <- -(leafIndex) # add new leaf
          leaf.labels <- c(leaf.labels, gene.names.ord[itemIndex])
          leafIndex <- leafIndex + 1
          itemIndex <- itemIndex + 1
          height[rowIndex] <- leaf.height
          newCluster <- F
        } else{
          if(dbg) cat("This is an old cluster (",g,")!\n")
          mergeMat[rowIndex, 1] <- rowIndex - 1 # previous row nr
          mergeMat[rowIndex, 2] <- -(leafIndex) # add new leaf
          leaf.labels <- c(leaf.labels, gene.names.ord[itemIndex])
          leafIndex <- leafIndex + 1
          itemIndex <- itemIndex + 1
          height[rowIndex] <- leaf.height
        }
        rowIndex <- rowIndex + 1
      }
    }
  }
  if(dbg) print(mergeMat)
  if(dbg) cat("Singleton names = ",singletons.n,"\n")
  if(dbg) cat("Singleton nrs   = ",singletons.g,"\n")
  if(dbg) cat("Singleton index = ",singletons.i,"\n")
  
  # STEP II. Make the stem connections between non-singleton clusters
  if(dbg) print("Step II mergeMat")
  if(n.clusters.NONSINGLE > 1){
    newCluster <- T
    for (i in 1 : (n.clusters.NONSINGLE-1)){ # i runs over all non-singleton clusters
      if(newCluster){
        mergeMat[rowIndex, 1] <- rowStem[i] # add new stem
        mergeMat[rowIndex, 2] <- rowStem[i+1] # add new stem
        height[rowIndex] <- stem.height
        newCluster <- F
      } else {
        mergeMat[rowIndex, 1] <- rowIndex - 1 # previous row nr
        mergeMat[rowIndex, 2] <- rowStem[i+1] # add new stem
        height[rowIndex] <- stem.height
      }  
      rowIndex <- rowIndex + 1
    }
  }
  if(dbg) print(mergeMat)
  
  # STEP 3. Connect the singleton clusters as single leaves
  if(dbg) print("Step III mergeMat")
  if(dbg) cat("n.clusters.SINGLE = ",n.clusters.SINGLE,"\n")
  if(n.clusters.SINGLE > 0){
    for(i in 1 : (n.clusters.SINGLE)){ # i runs over all singleton clusters)
      mergeMat[rowIndex, 1] <- rowIndex - 1 # previous row nr
      mergeMat[rowIndex, 2] <- -(leafIndex) # add new stem    
      height[rowIndex] <- stem.height
      leaf.labels <- c(leaf.labels, singletons.n[i])
      leafIndex <- leafIndex + 1
      rowIndex <- rowIndex + 1
    }
  }
  if(dbg) print(mergeMat)
  a$merge <- mergeMat
  
  # make height vector
  a$height <- height
  
  # make order vector
  # I. Indexing the non-singleton clusters:
  accumulate <- 0
  leafOrder <- c()
  for(i in 1 : n.clusters.NONSINGLE){
    thisSize <- size.clusters[NONSINGLE.g[i]]
    ndxStem <- round(thisSize/2)
    leafOrder <- c(leafOrder, accumulate + (1 : ndxStem))
    leafOrder <- c(leafOrder, accumulate + (thisSize : (ndxStem+1)))
    accumulate <- accumulate + thisSize
  }
  # II. Indexing singleton clusters:
  if(n.clusters.SINGLE > 0){
    leafOrder <- c(leafOrder, (accumulate+1) : n.genes)
  }
  if(dbg) cat("order = ",leafOrder,"\n")
  a$order <- leafOrder
  
  # make and reshuffle label vector
  a$labels <- leaf.labels[leafOrder]
  
  # make class
  class(a) <- "hclust"
  # return clustering object a
  return(a)
}


################################################################
#
# 3. function traditionalScreeplot
#
# Emile Apol, 24 march 2013
#
# Plots a "traditional" screeplot, 
# based on pca object with class "prcomp" or class "princomp".
#
# xNames = F: only PC index is used as x-label
#        = T: names of PC's are being used as x-label
#
################################################################

traditionalScreeplot <- function(x, xNames=F){
  nX <- length(x$sdev); xplot <- (1:length(x$sdev))
  plot(x$sdev^2 ~ xplot, type="b", pch=20, 
       axes=F, xlab="Component", 
       ylab="Eigenvalue", main="Screeplot")
  if(xNames==T){
    axis(1, at=(1 : nX), labels=colnames(x$rotation), las=2)
  } else{
    axis(1, at=(1 : nX), labels=xplot)
  }
  axis(2)
  box()
}

################################################################
#
# 4. function traditionalPCPlot
#
# Emile Apol, 24 march 2013
#
# Plots a "traditional" PC1-PC2-plot, i.e., individual cases 
# projected in the plane of the first 2 eigenvectors PC1 and PC2, 
# based on pca object with class "prcomp" or class "princomp".
#
# Text = T (F) logical: print labels at points in PC Plot?
# Select =   (optional) logical vector with selection of datapoints to be 
#            plotted
#
# col = (optional) vector of colors of the labels per point.
# main = (optional) title of plot
#
################################################################

traditionalPCPlot <- function(x, Text=T, Select=T, dbg=F, ...){
  if(class(x)=="prcomp"){
    X <- x$x[,1]
    Y <- x$x[,2]
    theNames <- rownames(x$x)
  } else if(class(x)=="princomp"){
    X <- x$scores[,1]
    Y <- x$scores[,2]
    theNames <- rownames(x$scores)
  }
  X <- X[Select]
  Y <- Y[Select]
  theNames <- theNames[Select]
  plotMax <- max(abs(X), abs(Y))
  nNames <- length(theNames)
  cexNames <- 1.0/(nNames^0.10)
  if(dbg) print(plotMax)
  plot(Y ~ X, type="n", pch=20, asp=1, 
       xlab="PC 1", ylab="PC 2",
       xlim=c(-plotMax, plotMax),
       ylim=c(-plotMax, plotMax), ...)
  points(Y ~ X, type="p", pch=20, ...)
  if(Text) text(Y ~ X, pos=3, labels=theNames, xpd=T, cex=cexNames)
  abline(h=0, lty="dashed", col="black")
  abline(v=0, lty="dashed", col="black")
}


#######################################################################
#
# 5. function makeColDendroG
#
#
# Emile Apol, based on Gaston Sanchez, 2012
# 5 april 2014
#
# Color labels a dendrogram according to a grouping vector g
# x = hclust object or a dendrogram object
#
# Returns: dendrogram object
#
#######################################################################

makeColDendroG <- function(x, g){
  
  # vector of colors
  # labelColors <- rainbow(k)
  labelColors <- palette("default")
  # make names to the grouping vector
  g <- as.factor(g)
  names(g) <- x$labels
  clusMember <- g
  # cut dendrogram in k clusters
  # clusMember <- cutree(x, k)
  # if necessary, change class of x into dendrogram
  xd <- as.dendrogram(x)
  # function to get color labels
  colLab <- function(n) {
    if (is.leaf(n)) {
      a <- attributes(n)
      labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
      attr(n, "nodePar") <- c(lab.col = labCol, pch="")
    }
    n
  }
  # using dendrapply
  clusDendro <- dendrapply(xd, colLab)
  return(clusDendro)
}


#######################################################################
#
# 6. function plotColDendroG
#
#
# Emile Apol, based on Gaston Sanchez, 2012
# 6 april 2014
#
# Color labels a dendrogram according to a grouping vector g
# x = hclust object or a dendrogram object
# ... = other graphical parameters (like xlab, ylab etc.)
#
# Returns: NULL (plot)
#
#######################################################################

plotColDendroG <- function(x, g, ...){
  
  # number of labels
  n.lab <- length(labels(x))
  n.grp <- length(g)
  # set cex (size of labels)
  cexLab <- 1.0/n.lab^0.10
  par(cex=cexLab)
  # vector of colors
  # labelColors <- rainbow(k)
  labelColors <- palette("default")
  # make names to the grouping vector
  g <- as.factor(g)
  names(g) <- x$labels
  clusMember <- g
  # cut dendrogram in k clusters
  # clusMember <- cutree(x, k)
  # if necessary, change class of x into dendrogram
  xd <- as.dendrogram(x)
  # function to get color labels
  colLab <- function(n) {
    if (is.leaf(n)) {
      a <- attributes(n)
      labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
      attr(n, "nodePar") <- c(lab.col = labCol, pch="")
    }
    n
  }
  # using dendrapply
  clusDendro <- dendrapply(xd, colLab)
  plot(clusDendro, ...)
  par(cex=1.0)
  return()
}





cat("Sourced: Dendrograms_v8.r\n")

test=F
if(test!=F){
  
######################################################################
#
# Test data sets for clusters / dendrograms
#
######################################################################

# Make test dataset
#
# cluster 1: A, B, D, E       size = 4    center = 0
# cluster 2: C, F, G, K, I    size = 5    center = 1
# cluster 3: H, J             size = 2    center = 2
# cluster 4: L                size = 1    center = 3
# cluster 5: M                size = 1    center = 4
kclust <- list()
clusters <- c(1, 1, 2, 1, 1, 2, 2, 3, 2, 3, 2, 4, 5)
names(clusters) <- paste("gene", LETTERS[1:length(clusters)])
( kclust$cluster <- clusters )
rm(clusters)
kclust$size <- c(4, 5, 2, 1, 1)
class(kclust) <- "kmeans"

clust <- makeKMeansDendrogram(kclust)
clust <- makeKMeansDendrogram(kclust, dbg=T)
plot(clust, hang=-1, ann=F, axes=F)
clust$merge

# Make test dataset
#
# a 2-dimensional example
X <- rbind(matrix(rnorm(30, sd = 0.3), ncol = 2),
           matrix(rnorm(30, mean = 1, sd = 0.3), ncol = 2))
colnames(X) <- c("M1", "M2")
rownames(X) <- paste("gene", 1:30, sep="")
View(X)
kclust <- kmeans(X, centers=2, nstart=5)
clust <- makeKMeansDendrogram(kclust, dbg=T)
plot(clust, hang=-1, ann=F, axes=F)
clust$merge

}