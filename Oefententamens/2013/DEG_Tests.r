##############################################################################
#
# This file contains several functions to apply statistical tests
# on rows (columns) of matrices containing microarray data (i.e., 
# logtransformed expression values). The methods are parametric 
# (t-tests, 1-way ANOVA) als well as non-parametric (Wilcoxon tests, 
# Kruskal-Wallis test).
#
# All functions return p-values, some functions return effect sizes
# as well (eta2, r2). For non-parametric tests, the effect size is
# bases on Fields (2005) - Discovering Statistics using SPSS.
#
# Emile Apol, feb-april 2013 
#
# ----------------------------------------------------------------------------
#
# Statistical test returning p-values:
#
# 1.  matrix1SampleTTest0       1-sample t-test for rows x of matrix X, 
#                               with mu = 0
# 2.  matrix1SampleTTest        1-sample t-test for rows x of matrix X, 
#                               with mu to be specified as m=...
# 3.  matrixPairedTTest         paired t-test for rows x of matrix X, 
#                               with order: n/2 level "1", n/2 level "2"
# 4.  matrix1SampleWilcoxTest0  1-sample Wilcoxon-test for rows x of matrix X, 
#                               with mu = 0
# 5.  matrix1SampleWilcoxTest   1-sample Wilcoxon-test for rows x of matrix X, 
#                               with mu to be specified as m=...
# 6.  matrixPairedWilcoxTest    paired Wilcoxon-test for rows x of matrix X, 
#                               with order: n/2 level "1", n/2 level "2"
# 7.  matrix1WayANOVA           1-way ANOVA-test for rows x of matrix X, 
#                               with order: g = grouping vector along row
# ----------------------------------------------------------------------------
#
# Statistical tests returning p-values and effect sizes:
#
# 8.  matrixTTest               (2-sample/Welch/paired) t-test for rows x of matrix X, 
#                               with order: g = grouping vector along row
# 9.  matrix1WayANOVATest       1-way ANOVA-test for rows x of matrix X, 
#                               with order: g = grouping vector along row
#
##############################################################################



#########################################################################
# 1. function matrix1SampleTTest0
#
# Emile Apol
# 25 feb 2013
#
# 1-sample t-test for rows x of matrix X, with mu = 0
#
# Usage:
# apply(X, 1, matrix1SampleTTest0)
#
#########################################################################

matrix1SampleTTest0 <- function(x){
  t.test(x, mu=0)$p.value
}

#########################################################################
# 2. function matrix1SampleTTest
#
# Emile Apol
# 27 feb 2013
#
# 1-sample t-test for rows x of matrix X, 
# with mu to be specified as m=...
#
# Usage:
# apply(X, 1, matrix1SampleTTest, m=5)
# apply(X=X, MARGIN=1, FUN=matrix1SampleTTest, m=5)
#
#########################################################################

matrix1SampleTTest <- function(x, m){
  
    t.test(x, mu=m)$p.value
  
}

#########################################################################
#
# 3. function matrixPairedTTest 
#
# Emile Apol
# 5 march 2013
#
# paired t-test for rows x of matrix X, 
# with order: n/2 level "1", n/2 level "2"
#
# Usage:
# apply(X, 1, matrixPairedTTest)
# 
#########################################################################

matrixPairedTTest <- function(x){
  
  N.elem <- length(x)
  # number of paired = half of the length
  N.pairs <- N.elem / 2 
  # make a factor with groups
  mySample <- factor(c(rep(1, N.pairs), rep(2, N.pairs)))
  t.test(x ~ mySample, paired=T)$p.value  

}

#########################################################################
# 4. function matrix1SampleWilcoxTest0
#
# Emile Apol
# 6 march 2013
#
# 1-sample Wilcoxon-test for rows x of matrix X, with mu = 0
#
# Usage:
# apply(X, 1, matrix1SampleWilcoxTest0)
#
#########################################################################

matrix1SampleWilcoxTest0 <- function(x){
  wilcox.test(x, mu=0)$p.value
}

#########################################################################
# 5. function matrix1SampleWilcoxTest
#
# Emile Apol
# 6 march 2013
#
# 1-sample Wilcoxon-test for rows x of matrix X, 
# with mu to be specified as m=...
#
# Usage:
# apply(X, 1, matrix1SampleWilcoxTest, m=5)
# apply(X=X, MARGIN=1, FUN=matrix1SampleWilcoxTest, m=5)
#
#########################################################################

matrix1SampleWilcoxTest <- function(x, m){
  
  wilcox.test(x, mu=m)$p.value
  
}

#########################################################################
#
# 6. function matrixPairedWilcoxTest 
#
# Emile Apol
# 5 march 2013
#
# paired Wilcoxon-test for rows x of matrix X, 
# with order: n/2 level "1", n/2 level "2"
#
# Usage:
# apply(X, 1, matrixPairedWilcoxTest)
# 
#########################################################################

matrixPairedWilcoxTest <- function(x){
  
  N.elem <- length(x)
  # number of paired = half of the length
  N.pairs <- N.elem / 2 
  # make a factor with groups
  mySample <- factor(c(rep(1, N.pairs), rep(2, N.pairs)))
  wilcox.test(x ~ mySample, paired=T)$p.value  
  
}

#########################################################################
#
# 7. function matrix1WayANOVA 
#
# Emile Apol
# 28 march 2013
#
# 1-way ANOVA-test for rows x of matrix X, 
# with order: g = grouping vector along row
#
# Usage:
# apply(X, 1, matrix1WayANOVA, g=)
# 
#########################################################################

matrix1WayANOVA <- function(x, g){
  
  g <- as.factor(g)
  summary(aov(x ~ g))[[1]]$Pr[1]
  
}


#########################################################################
#
# 8. function matrixTTest 
##
# Emile Apol
# 3 april 2013
#
# (2-sample/Welch/paired) t-test for rows x of matrix X, 
# with order: g = grouping vector along row
#
# Usage:
# apply(X, 1, matrixTTest, g= ,[var.equal=T/F, paired=T,...])
#
# Output: vector a with
# a[1] = p.value (statistical significance)
# a[2] = eta^2 (effect size)
# 
#########################################################################

matrixTTest <- function(x, g, ...){
  
  g <- as.factor(g)
  Q <- t.test(x ~ g, ...)
  p.value <- Q$p.value
  eta2 <- Q$statistic^2/(Q$statistic^2 + Q$parameter)
  a <- c()
  a[1] <- p.value
  a[2] <- eta2
  names(a) <- c("p-value","eta2")
  return(a)
}

#########################################################################
#
# 9. function matrix1WayANOVATest 
#
# Emile Apol
# 3 april 2013
#
# 1-way ANOVA-test for rows x of matrix X, 
# with order: g = grouping vector along row
#
# Usage:
# apply(X, 1, matrix1WayANOVATest, g=)
#
# Output: vector a with
# a[1] = p.value (statistical significance)
# a[2] = eta^2 (effect size)
# 
#########################################################################

matrix1WayANOVATest <- function(x, g){
  
  g <- as.factor(g)
  Q <- summary(aov(x ~ g))[[1]]
  p.value <- Q$Pr[1]
  SS.g <- Q$Sum[1]; SS.tot <- sum(Q$Sum); eta2 <- SS.g/SS.tot
  a <- c()
  a[1] <- p.value
  a[2] <- eta2
  names(a) <- c("p-value","eta2")
  return(a)
  
}


cat("Sourced: DEG_Tests_v10.r")

test=F
if(test){
  
# Some test data sets
  ( x <- c(1.0,1.1,1.2,0.9,0.8,0.7) )
  ( X <- matrix(c(1.0,1.1,1.2,0.9,0.8,0.7,
                0.5,0.4,0.6,0.5,0.5,0.7,
                0.2,0.3,0.1,0.1,0.2,0.3),
              byrow=T, nrow=3) )
  rownames(X) <- paste("Gene", LETTERS[seq_along(X[,1])], sep="-")
  colnames(X) <- c("M.1.1","M.1.2","M.1.3","M.2.1","M.2.2","M.2.3")
  X
  ( samples <- factor(rep(c(1,2), each=3)) )
  
  t.test(x ~ samples, var.equal=F)
  apply(X, 1, matrixTTest, g=samples, var.equal=F) # Welch t-test
  apply(X, 1, matrixTTest, g=samples, var.equal=T) # 2-sample t-test
  apply(X, 1, matrixTTest, g=samples, paired=T)    # paired t-test
  wilcox.test(x ~ samples)
  apply(X, 1, matrixWilcoxTest, g=samples)            # 2-sample Wilcox-test
  apply(X, 1, matrixWilcoxTest, g=samples, paired=T)  # paired Wilcox-test
  
  # Pitman example (see Kruskall & Wallis, 1952)
  x <- c(0,11,12,20,16,19,22,24,29)
  samples <- factor(c(1,1,1,1,2,2,2,2,2))
  X <- matrix(x,byrow=T, nrow=1)
  apply(X, 1, matrixWilcoxTest, g=samples)
  wilcox_test(x ~ samples)
  
  # Brownlee example (see Kruskall & Wallis, 1952)
  x <- c(95.6, 94.9, 96.2, 95.1, 95.8, 96.3,
         93.3, 92.1, 94.7, 90.1, 95.6, 90.0, 94.7)
  samples <- factor(c(1,1,1,1,1,1,2,2,2,2,2,2,2))
  X <- matrix(x,byrow=T, nrow=1)
  apply(X, 1, matrixWilcoxTest, g=samples)
  wilcox_test(x ~ samples)

}
