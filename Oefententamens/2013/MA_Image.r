##############################################################################
#
# This file contains several functions to help imaging microarray data.
#
# Emile Apol, feb-april 2013 
#
# ----------------------------------------------------------------------------
#
# 1.  myPlotMA              Plots a matrix/dataframe of microarray data as a nice image with 
#                           gene and sample names in logical order.
# 2.  MA.colors             Makes a standard MA Green-Black-Red color palette
# 3.  MA.colors.2           Makes a standard MA Green-Black-Red color palette
#
#################################################################


#################################################################
#
# 1. function myPlotMA
#
# Emile Apol, 15 march 2013
# 
# Plots a matrix/dataframe of microarray data as a nice image with 
# gene and sample names in logical order.
#
# x = dataframe
# coll = color palette (opt)
# scaleLabels = F (opt: T): scale the size of the labels according to the 
#         number of labels?
#
################################################################

myPlotMA <- function(x, coll=heat.colors(12), scaleLabels=T){
  # store old graphics parameters and set new margins:
  oldpar <- par(no.readonly=T)
  par(mar=c(5,5,2,7))
  # change dataframe into matrix:
  M <- as.matrix(x)
  # calculate scaling of labels (cex.lab)
  nX <- ncol(M); nY <- nrow(M) # nr of X and Y points
  nMax <- max(nX, nY); axisScale <- 1.0*nMax^(-0.1)
  # first make an image without axes:
  if(scaleLabels) par(cex.axis=axisScale)
  image(t(M)[, nY:1], col=coll,
        xlab="Sample", ylab="Gene", axes=F)
  # add x and y axes:
  axis(1, at = seq(0, 1, length.out=nX), 
       labels=colnames(M), tick=F, las=1)
  axis(4, at = seq(0, 1, length.out=nY), 
       labels=rownames(M)[nY:1], tick=F, las=1)
  par(oldpar) # reset graphics parameters
}

##################################################################
#
# 2. function MA.colors
#
# Emile Apol
# 15 march 2013
#
# Makes a standard MA Green-Black-Red color palette in RGB style: 0 (min) .. 255 (max)
# negative M-values:  green
# zero M-values:      black
# positive M-values:  red
#
##################################################################

MA.colors <- function(n=12){
  
  halve.1 <- trunc(n/2)
  halve.2 <- n - halve.1
  RGBMat <- data.frame(red=c( rep(0, halve.1), trunc(seq(0, 255,length.out=halve.2)) ),
                       green=c( trunc(seq(255, 0, length.out=halve.1)), rep(0, halve.2) ),
                       blue=rep(0, n))
  return( rgb(RGBMat, maxColorValue=255) )
}

##################################################################
#
# MA.colors.2
#
# Marcel Kempenaar
# march 2013
#
# Make a standard MA color palette in RGB style: 0 (min) .. 255 (max)
# negative M-values:  green
# zero M-values:      black
# positive M-values:  red
#
##################################################################

MA.colors.2 <- function(n=12){
  colorRampPalette(c("green", "black", "red"), space="rgb")(n)
}

cat("Sources: MA_Image_v3.r\n")
