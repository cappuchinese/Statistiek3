##########################################################################################
#
# Functions to perform post-hoc tests after 1-way ANOVA analysis
#
##########################################################################################


#########################################################################################
#
# function postHocHomSubsets(y, g, p.adjust.method, alpha)
#
# Perform a post-hoc test on data y and a SINGLE grouping factor g, 
# using the function pairwise.t.test with p.adjust.method as method for 
# multiple testing correction (for Fisher's LSD, Bonferroni, Holm or FDR methods), 
# or evaluating the Tukey-Kramer method explicitly (for Tukey/Tukey-Kramer method).
# Next, create homogeneous subsets based on the resulting (adjusted) p-values 
# similar to SPSS output, selecting "significant" differences based on alpha.
#
# Emile Apol
# April 11, 2014
#
# Current options of multiple-testing corrections:
# "none" = "LSD"    - Fisher's Least Significant Difference method (default)
# "bonferroni"      - Bonferroni method
# "holm"            - Holm method (more general than Bonferroni)
# "fdr" = "BH"      - False Discovery Rate method of Benjamini & Hochberg (1995)
# "tukey" = "HSD"   - Tukey-Kramer Honest Significant Difference method
#
#
#########################################################################################

postHocHomSubsets <- function(y, g, p.adjust.method="none", alpha=0.05){
  
  if(p.adjust.method=="none"){
    METH <- "PAIRWISE"
  }
  else if(p.adjust.method=="LSD"){
    METH <- "PAIRWISE"
    p.adjust.method <- "none"
  }
  else if(p.adjust.method=="bonferroni"){
    METH <- "PAIRWISE"
  }
  else if(p.adjust.method=="holm"){
    METH <- "PAIRWISE"
  }
  else if(p.adjust.method=="fdr"){
    METH <- "PAIRWISE"
  }
  else if(p.adjust.method=="BH"){
    METH <- "PAIRWISE"
    p.adjust.method=="fdr"
  }
  else if(p.adjust.method=="tukey"){
    METH <- "TUKEY"
  }
  else if(p.adjust.method=="HSD"){
    METH <- "TUKEY"
  }
  
  n.g <- nlevels(g)                   # number of levels
  means <- tapply(y, g, mean)         # mean per group
  ORDER <- order(means)
  levels.ordered <- levels(g)[ORDER]
  means.ordered <- means[ORDER]
  
  switch(METH,
         
         # PAIRWISE correction methods (LSD, Bonferroni, Holm, FDR)
         
         PAIRWISE={
           # perform post-hoc test with appropriate multiple-testing correction
           res.PH <- pairwise.t.test(y, g, p.adjust.method=p.adjust.method, pool.sd=T)
           # and print result to screen
           print(res.PH)
           
           # create full p-values matrix
           pMat <- matrix(rep(1, n.g^2), ncol=n.g)
           rownames(pMat) <- levels(g)
           colnames(pMat) <- levels(g)
           # fill the p-matrix
           for(r in rownames(res.PH$p.value)){
             for(k in colnames(res.PH$p.value)){
               z <- as.numeric(res.PH$p.value[r,k])
               if(!is.na(z)){
                 pMat[r,k] <- z
                 pMat[k,r] <- z
               }
             }
           }
           # reorder p-matrix
           pMat.ordered <- pMat[levels.ordered, levels.ordered]
           rm(pMat)
         },
         
         # TUKEY correction methods (Tukey)
         
         TUKEY={
           vars <- tapply(y, g, var)           # var per group
           ns <- tapply(y, g, length)          # n per group
           var.p <- sum((ns-1)*vars)/sum(ns-1) # pooled variance
           s.p <- sqrt(var.p)                  # pooled standard deviation
           df.2 <- sum(ns-1)                   # df's of s.p
           df.1 <- n.g
           
           if(0){
             cat("ns   = ",ns,"\n")
             cat("vars = ",vars,"\n")
             cat("s.p  = ",s.p,"\n")
             cat("df.1 = ",df.1," df.2 = ",df.2,"\n")
           }
           
           absDiffsMat <- abs(outer(means, means, FUN="-"))
           SEMat <- s.p*sqrt( (1/2) * outer(1/ns, 1/ns, FUN="+") )
           qMat <- absDiffsMat / SEMat
           
           if(0){
             cat("absDiffsMat = \n")
             print(absDiffsMat)
             cat("SEMat =\n")
             print(SEMat)
             cat("qMat = \n")
             print(qMat)
           }
           
           # 1-sided p-values
           pMat <- apply(qMat, c(1,2), ptukey, nmeans=df.1, df=df.2, lower.tail=F)
           if(0){
             # 2-sided p-values
             pMat <- 2*pMat
             pMat <- apply(pMat, c(1,2), min, 1)
           }
           rownames(pMat) <- levels(g)
           colnames(pMat) <- levels(g)
           
           # print to screen
           cat("\n\n Tukey-Kramer post-hoc test:\n\n")
           print(pMat)
           
           # reorder p-matrix
           pMat.ordered <- pMat[levels.ordered, levels.ordered]
           rm(pMat)
           
           # stop("Tukey-Kramer HSD method: Not yet implemented...")
         },
         
         # ERROR
         
         stop("Error on switch in method METH ...")
         )
  
  # create homogeneous subsets
  n.homsubs <- 1
  bFirst <- T
  data.homsub <- data.frame(subset.1=rep(NA, nlevels(g)))
  rownames(data.homsub) <- levels.ordered
  for(r in rownames(pMat.ordered)){
    if(bFirst==T){
      for(k in colnames(pMat.ordered)){
        if(pMat.ordered[r, k] >= alpha){
          data.homsub[k, n.homsubs] <- means.ordered[k]
        } # if p > alpha
      } # colomns k
      bFirst <- F
      first.k <- 1
    } # bFirst
    else{ 
      # only make new homogeneous subset if start at new column!
      if(pMat.ordered[r, first.k] < alpha){
        # make new subset
        n.homsubs <- n.homsubs + 1
        newName <- paste("subset.", n.homsubs, sep="")
        data.homsub <- data.frame(data.homsub,
                                  rep(NA, nlevels(g)))
        colnames(data.homsub)[n.homsubs] <- newName
        # first set new first.k value
        k.nr <- 0
        for(k in colnames(pMat.ordered)){
          k.nr <- k.nr + 1
          if(pMat.ordered[r, k] >= alpha){
            first.k <- k.nr
            break 
          } # if p > alpha
        } # colomns k
        # and then fill dataframe
        for(k in colnames(pMat.ordered)){
          if(pMat.ordered[r, k] >= alpha){
            data.homsub[k, n.homsubs] <- means.ordered[k]  
          } # if p > alpha
        } # colomns k      
      } # if make new subset    
    } # bFirst
  } # for rows r
  
  # and print result to screen
  cat("\n\nHomogeneous subsets:\n\n")
  print(as.matrix(data.homsub), na.print="")
  cat("\n(Means per level per subset)\n")
}

cat("Sourced: Homogeneous_subsets_v1.r\n")

test <- F

if(test==T){
  ( yy <- c(1, 2, 1, 3, 1, 2, 3, 4, 5, 5, 3, 4, 3, 6, 2, 2, 3) )
  ( gg <- factor(c(rep("Algae", 6), rep("Bears", 4), rep("Crocs", 4), rep("Dinos", 3))) )
  
  postHocHomSubsets(yy, gg, p.adjust.method="none")
  postHocHomSubsets(yy, gg)
  postHocHomSubsets(yy, gg, p.adjust.method="bonferroni")
  postHocHomSubsets(yy, gg, p.adjust.method="fdr")
  postHocHomSubsets(yy, gg, p.adjust.method="holm")
  TukeyHSD(aov(yy ~ gg))
  postHocHomSubsets(yy, gg, p.adjust.method="HSD")
}


