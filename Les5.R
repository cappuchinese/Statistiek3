#
# Statistiek 3 BIN 2015-2016
#
# Les 05
#

# Make 10 000 random normally distributed numbers
N <- 10000
mu <- 5.0
sigma <- 1.5
y <- rnorm(N, mean = mu, sd = sigma)

# Make a histogram
hist(y, breaks=40, freq = F)
curve(dnorm(x, mean = mu, sd = sigma), from=mu-6*sigma, to=mu+6*sigma, col="red", add=T)


# Make 1000 datasets of 10 points, and evaluate the 1000 averages
n <- 10

M <- matrix(y, ncol = n, byrow = T)
means <- apply(M, 1, mean)

# OR

n.sets <- N/n
g <- factor(rep(1:n.sets, each=n))
means <- tapply(y,g,mean)


hist(means, freq = F, breaks = 40)
curve(dnorm(x, mean = mu, sd = sigma/sqrt(n)), from = 3.5, to = 6.5, col="red", add = T)

stdevs <- tapply(y,g,sd)
s.mean <- stdevs/sqrt(n)

t <- (means - mu)/s.mean

hist(t, freq = F, breaks = 40)
curve(dt(x, df = n-1), from =-5, to = 5, col="red", add=T)
curve(dnorm(x, mean=0, sd=1), from =-5, to = 5, col="blue", add=T)
