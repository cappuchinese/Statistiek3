##############################################################################
#
# This file contains several functions for Gene Ontology & classification.
#
# Emile Apol, feb-april 2013 
#
# ----------------------------------------------------------------------------
#
# 1.  GO.binom.test       Calculates the binomial approximation of the Fisher exact test.
# 2.  GO.tests            Performs several tests (Fisher's exact test, Chisquared test 
#                         and binomial test) on occurance of trait F in N genes, with 
#                         K genes classified as DEG (or subcluster), x the number of DEG genes 
#                         with trait F and M the number of genes (out of the N genes) with trait F.
# 3.  print.GOtest        Define a new S3 print method for a GOtest object, as created by GO.tests().
#
##############################################################################


################################################################################
#
# 1. function GO.binom.test
#
# Emile Apol
# 27 march 2013
#
# Calculates the binomial approximation of the Fisher exact test.
#
# sided = (optional) sidedness of test ("two.sided" (def), "greater", "less")
#
################################################################################

GO.binom.test <- function(X, sided=c("two.sided", "greater", "less"), dbg=F){
  x <- X[1,1]
  n <- sum(X[1,])
  M <- sum(X[,1])
  N <- sum(X)
  p <- M/N
  sided <- match.arg(sided) # checks if correct option has been given
  if(dbg) cat("x = ", x, "n = ", n, "M = ", M, "N = ", N, "p = ", p, "\n")
  if(sided=="two.sided"){
    bin <- binom.test(x, n, p, alternative="two.sided")
  } else if(sided=="greater"){
    bin <- binom.test(x, n, p, alternative="greater")
  } else if(sided=="less"){
    bin <- binom.test(x, n, p, alternative="less")
  }
  return(bin)
}

################################################################################
#
# 2. function GO.tests
#
# Emile Apol
# 27 march 2013
#
# Performs several tests (Fisher's exact test, Chisquared test 
# and binomial test) on occurance of trait F in N genes, with 
# K genes classified as DEG (or subcluster), x the number of DEG genes 
# with trait F and M the number of genes (out of the N genes) with trait F.
#
# All tests are performed 2-sided (Enrichment/Depletion) and - if
# appropriate - 1-sided (Enrichment or Depletion), yielding the appropriate
# 2- or 1- sided p-values. Note that the R-functions binom.test() and
# fisher.test() provide 2-sided "minimum likelihood" p-values, 
# not mid-p-values! 
# See also: Rivals et al. (2007) Bioinformatics 23(4), 401-407. 
#
# The datamatrix X is assumed to have the following layout:
#
# counts    |   F         non F     |   tot
# -----------------------------------------
# DEG       |   x         (K-x)     | K
# non DEG   |   (M-x)     (N-M-K-x) | (N-K)
# -----------------------------------------
# tot       |   M         (N-M)     | N
#
################################################################################

GO.tests <- function(X, dbg=F){
  x <- X[1,1]
  n <- sum(X[1,])
  M <- sum(X[,1])
  N <- sum(X)
  p <- M/N
  OR <- (X[1,1]/X[1,2])/(X[2,1]/X[2,2])
  log2OR <- log2(OR)
  a <- list()
  fish <- fisher.test(XMat) 
  fish.G <- fisher.test(XMat, or=1, alternative="greater")
  fish.L <- fisher.test(XMat, or=1, alternative="less")
  chi.T <- chisq.test(XMat, correct=T)   
  chi.F <- chisq.test(XMat,correct=F)  
  bin <- binom.test(x, n, p, alternative="two.sided")
  bin.G <- binom.test(x, n, p, alternative="greater")
  bin.L <- binom.test(x, n, p, alternative="less")
  p.values <- c(fish$p.value, fish.G$p.value, fish.L$p.value, 
                chi.T$p.value, chi.F$p.value, 
                bin$p.value, bin.G$p.value, bin.L$p.value)
  names(p.values) <- c("Fisher", "Fisher (enrich)", "Fisher (deplete)", 
                       "Chisq test (Yates corr)", "Chisq test", 
                       "Binomial", "Binomial (enrich)", "Binomial (deplete)")
  pMat <- matrix(rep(NA,9), nrow=3)
  pMat[1,1] <- fish$p.value; pMat[1,2] <- fish.G$p.value; pMat[1,3] <- fish.L$p.value
  pMat[2,1] <- chi.T$p.value
  pMat[3,1] <- bin$p.value;  pMat[3,2] <- bin.G$p.value; pMat[3,3] <- bin.L$p.value
  rownames(pMat) <- c("Fisher", "Chisq", "Binomial")
  colnames(pMat) <- c("enrich/deplete (2-sided)", "enrich (1-sided)", "deplete (1-sided)")
  a$p.values <- pMat
  a$OR <- OR
  a$log2OR <- log2OR
  a$observed <- chi.F$observed
  a$expected <- chi.F$expected # E > 5 to use chi.sq approximation
  a$residuals <- chi.F$residuals
  chisq.check <- min(chi.F$expected)
  names(chisq.check) <- "Must be > 5 (10) to use chisq.test..."
  a$chisq.check <- chisq.check # E > 5 to use chi.sq approximation
  binom.check <- N/n
  names(binom.check) <- "Must be > 20 to use binom.test..." 
  a$binom.check <- binom.check # N/K must be > 20 to use binomial approximation
  a$method <- "GO Enrichment and/or Depletion tests (2- and 1-sided):\n   Fisher, Chisq (with Yates correction), Binomial"
  a$data.name <- "Fisher (E/D), Chisq (E/D), Binom (E/D), Fisher (E), -, Binom (E), Fisher (D), -, Binom(D)"
  class(a) <- "GOtest"
  return(a)
}

################################################################################
#
# 3. function print.GOtest
#
# Emile Apol, 27 march 2013
#
# Define a new S3 print method for a GOtest object, as created by GO.tests().
#
################################################################################

print.GOtest <- function(x){
  cat("\n", x$method,"\n\n")
  cat("P-values using various tests:\n")
  print(x$p.values, na.print="-")
  cat("\nlog2(OR) = ",x$log2OR,"\n")
  cat("\nAsymptotic distributions (Chisq, Binom): diagnostics\n")
  cat("Chisq: min(E) = ",x$chisq.check,"(must be > 5...)\n")
  cat("Binom: N/K = ",x$binom.check,"(must be > 20...)\n")
}


cat("Sourced: GO_Classification_v3.r\n")

test=F
if(test){
####################################################################
#
# Test data sets
#
####################################################################

# Ia. Test data set of number of DEG's for a certain trait
# counts    |   F       non F   |   tot
# ---------------------------------------
# DEG       |   5       10      | 15
# non DEG   |  95      1390     | 1485
# ---------------------------------------
# tot       |  100     1400     | 1500
#
XMat <- matrix(c(5, 10, 
                 95, 1390), nrow=2, byrow=T)
rownames(XMat) <- c("DEG", "non DEG")
colnames(XMat) <- c("F", "non F")
View(XMat)

# Ib. Test data set of number of DEG's for a certain trait
# counts    |   F       non F   |   tot
# ---------------------------------------
# DEG       |  10       5       | 15
# non DEG   |  1390     95      | 1485
# ---------------------------------------
# tot       |  1400     100     | 1500
#
XMat <- matrix(c(10, 5, 
                 1390, 95), nrow=2, byrow=T)
rownames(XMat) <- c("DEG", "non DEG")
colnames(XMat) <- c("F", "non F")
View(XMat)

# II. Test data set of number of DEG's for a certain trait
# counts    |   F       non F   |   tot
# ---------------------------------------
# DEG       |   15       0      | 15
# non DEG   |   95      1390    | 1485
# ---------------------------------------
# tot       |  110     1390     | 1500
#
XMat <- matrix(c(15, 0, 
                 95, 1390), nrow=2, byrow=T)
rownames(XMat) <- c("DEG", "non DEG")
colnames(XMat) <- c("F", "non F")
View(XMat)

# III. Test data set of number of DEG's for a certain trait
# counts    |   F       non F   |   tot
# ---------------------------------------
# DEG       |   300     100     | 400
# non DEG   |   900     200     | 1100
# ---------------------------------------
# tot       |  1200     300     | 1500
#
XMat <- matrix(c(300, 100, 
                 900, 200), nrow=2, byrow=T)
rownames(XMat) <- c("DEG", "non DEG")
colnames(XMat) <- c("F", "non F")
View(XMat)

# IVa. Test data set of number of DEG's for a certain trait
# Rivals et al. (2007) Bioinformatics 23(4), pp 401-407, $6.1. 
# counts    |   F       non F   |   tot
# ---------------------------------------
# DEG       |   4       3       | 7
# non DEG   |   2       11      | 13
# ---------------------------------------
# tot       |   6       14      | 20
#
XMat <- matrix(c(4, 3, 
                 2, 11), nrow=2, byrow=T)
rownames(XMat) <- c("DEG", "non DEG")
colnames(XMat) <- c("F", "non F")
View(XMat)

# IVb. Test data set of number of DEG's for a certain trait
# Rivals et al. (2007) Bioinformatics 23(4), pp 401-407, $6.2. 
# counts    |   F       non F   |   tot
# ---------------------------------------
# DEG       |   10      30      | 40
# non DEG   |   90      670     | 760
# ---------------------------------------
# tot       |   100     700     | 800
#
XMat <- matrix(c(10, 30, 
                 90, 670), nrow=2, byrow=T)
rownames(XMat) <- c("DEG", "non DEG")
colnames(XMat) <- c("F", "non F")
View(XMat)

# IVc. Test data set of number of DEG's for a certain trait
# Rivals et al. (2007) Bioinformatics 23(4), pp 401-407, $6.3. 
# counts    |   F       non F   |   tot
# ---------------------------------------
# DEG       |   100     900     | 1000
# non DEG   |   1900    22100   | 24000
# ---------------------------------------
# tot       |   2000    23000   | 25000
#
XMat <- matrix(c(100, 900, 
                 1900, 22100), nrow=2, byrow=T)
rownames(XMat) <- c("DEG", "non DEG")
colnames(XMat) <- c("F", "non F")
View(XMat)

# VI. Test data set of number of DEG's for a certain trait
# ppt lesson 12
# counts    |   F       non F   |   tot
# ---------------------------------------
# DEG       |   90      10      | 100
# non DEG   |   1810    90      | 1900
# ---------------------------------------
# tot       |   1900    100     | 2000
#
XMat <- matrix(c(90, 10, 
                 1810, 90), nrow=2, byrow=T)
rownames(XMat) <- c("DEG", "non DEG")
colnames(XMat) <- c("F", "non F")
View(XMat)

######################################
# data analysis
######################################

chisq.test(XMat, correct=T)

( res <- GO.tests(XMat) )
res$p.values
res$chisq.check
res$binom.check
res$log2OR
res$OR

fisher.test(XMat)

fisher.test(XMat, or=1, alternative="greater")

fisher.test(XMat, or=1, alternative="less")


}
