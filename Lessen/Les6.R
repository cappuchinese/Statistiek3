# Statistiek 3 (BIN)
#
# Les 06
#
# 2014-2015

# Maak zelf dataserie van normaal verdeelde waarden

mu    <- 10       # echte gemiddelde
sigma <- 2        # echte standaarddeviatie
n     <- 4        # metingen per dataset
n.col <- 1000     # aantal datasets
n.tot <- n*n.col  # totaal aantal metingen

# Maak n.tot random normaalverdeelde variabelen
y <- rnorm(n = n.tot, mean = mu, sd = sigma)
# Bereken experimentele gemiddelde en standaarddeviatie
( y.mean <- mean(y) )
( s <- sd(y) )


# Verdeel de data in matrix met n rijen (= metingen); elke kolom is dus 1 dataset:
M <- matrix(y, nrow = n)

# Gemiddelden per kolom (= n waarden)
( y.mean.s <- apply(M, 2, mean) )
# Bereken experimentele fout in gemiddelde op basis van n metingen:
( s.e <- s/sqrt(n) )

# bereken t-waarde voor elke dataset van n metingen
( sd.s <- apply(M, 2, sd) )                # vector met standaarddeviaties per kolom
( n.s <- apply(M, 2, length) )             # vector met aantallen per kolom
( se.s <- sd.s/sqrt(n.s) )                 # vector met standaardfouten per kolom 
( t.s <- (y.mean.s - mu)/se.s )            # vector met t-waarden per kolom
( span <- ceiling(max(abs(min(t.s)), abs(max(t.s)))) )  # afgerond grootste uitwijking
hist(t.s, breaks = seq(-span, span, 0.5), freq = F, xlab="t-waarde")
# plot t-verdeling
x <- seq(-span, span, 0.1)                 # dummy variabele voor plotten
p.t <- dt(x, df = n-1)                     # kansen volgens t-verdeling
points(p.t ~ x, type="l", col="red")       # voeg punten toe aan grafiek
# en de corresponderende normaalverdeling (klopt NIET!):
p.norm <- dnorm(x, mean = 0, sd = s.e)     # kansen volgens normaalverdeling
points(p.norm ~ x, type="l", col="blue")   # voeg punten toe aan grafiek

alpha <- 0.05  # kans op foute conclusie
( CI.s <- qt(p = 1-alpha/2, df = n.s-1) * se.s )  # CI = t * fout in gemiddelde
( lower <- y.mean.s - CI.s )
( upper <- y.mean.s + CI.s )

# Maak een scatterplot met foutenstreepjes (dit is NIET standaard in R!!!)
ind <- 1:n.col                             # nummer van dataset
plot(y.mean.s ~ ind, ylim = range(lower, upper), pch = 19, 
     col = "black", xlab = "Dataset", ylab = "Gemiddelde",
     xlim = c(0, 50))
arrows(ind, lower, ind, upper, length = 0.05, angle = 90, code = 3)
abline(mu,0, col = "red")

# in hoeveel gevallen ligt mu NIET binnen het betrouwbaarheidsinterval?
( n.verkeerd <- sum(((lower < mu) & (upper < mu)) | ((lower > mu) & (upper > mu))) )
( kans.verkeerd <- n.verkeerd/n.col )


### Herhaling vorige les: 1-sample t-toets, 2-sample t-toets,
### Welch-t-toets, F-toets


##########################################################################
#
# t-toetsen : 1-sample t-toets
#
##########################################################################

# Gemeten lengte van 6 meisjes:
y.1 <- c(165, 167, 161, 168, 171, 177)
( y.av <- mean(y.1) )

# Q: Is het gemiddelde ANDERS dan 170 cm?

( s <- sd(y.1) )
( n <- length(y.1) )
( mu <- 170 )
( t <- (y.av - mu)/(s/sqrt(n)))
( t.krit <- qt(p = 0.975, df = n-1))

t.test(y.1, mu = 170)

# Q: Is het gemiddelde LAGER dan 170 cm?

t.test(y.1, mu = 170, alternative = "less")

# Q: Is het gemiddelde HOGER dan 170 cm?

t.test(y.1, mu = 170, alternative = "greater")

###############
# Opdracht
###############

# Lengte van 7 jongens:
y.2 <- c(175, 177, 184, 169, 186, 178, 182)

# Q: Is het gemiddelde ANDERS dan 170 cm?


# Q: Is het gemiddelde LAGER dan 170 cm?


# Q: Is het gemiddelde HOGER dan 170 cm?


##########################################################################
#
# t-toetsen : 2-sample t-toets
#
##########################################################################

# Q: Is de gemiddelde lengte van meisjes ANDERS dan de gemiddelde lengte van jongens?

mean(y.1)
mean(y.2)

boxplot(y.1, y.2)

var(y.1)
var(y.2)

# De varianties zijn behoorlijk gelijk -> 2-sample t-toets (var.eqal = T)

t.test(y.1, y.2, var.equal = T)

# Q: Zijn meisjes gemiddeld KLEINER dan jongens?

t.test(y.1, y.2, var.equal = T, alternative = "less")

# Q: Zijn meisjes gemiddeld GROTER dan jongens?

t.test(y.1, y.2, var.equal = T, alternative = "greater")



# Alternatieve syntax: met een groepsaanduiding = factor

y <- c(y.1, y.2)
g <- factor(c(rep("meisje", length(y.1)), 
              rep("jongen", length(y.2))))

plot(y ~ g)

# Q: Is de gemiddelde lengte van meisjes ANDERS dan de gemiddelde lengte van jongens?

t.test(y ~ g, var.equal = T)

##################
# Opdracht
##################

# Q: Zijn meisjes gemiddeld KLEINER dan jongens? Voer de t-test uit met y ~ g



# Q: Zijn meisjes gemiddeld GROTER dan jongens? Voer de t-test uit met y ~ g



############################################################################
#
# F-toets
#
############################################################################

# Zijn de varianties (standaarddeviaties) van 2 datasets gelijk?

var(y.1)
var(y.2)

var.test(y.1, y.2)
var.test(y ~ g)

# Je ziet dat bijna alle opdrachten dezelfde syntax hebben:
plot(y ~ g)
var.test(y ~ g)
t.test(y ~ g, ...)



### Verder met les 6
###

# Klap de matrix M om, zodat we nu 1000 regels met 4 datapunten hebben
# (transponeneren)

M
M <- t(M)
M

# We willen nu per regel (= dataset) een t-toets uitvoeren om te onderzoeken
# of het gemiddelde significant afwijkt van mu = 10.

# Per regel in matrix: apply(M, 1, functie)

# Definieer daarom eerst eigen functie die dit doet... De functie,
# bijv. my.t.test levert per regel de p-waarde

# Stel:
y <- c(1, 3, 2, 2, 3)
test <- t.test(y, mu = 1)
test

# Q: Hoe krijg ik uit test de p-waarde?


# Maken van de functie my.t.test

my.t.test <- function(x, mu0){
  test <- t.test(x, mu=mu0)
  return(test$p.value)
}

pVals <- apply(M, 1, my.t.test, mu0=10)

# Hoeveel datasets wijken af van mu = 10?

sum(pVals < alpha)


# Q: Hoeveel datasets wijken af van mu = 10.2?



# Q: Hoeveel datasets wijken af van mu = 12?




sum(apply(M, 1, my.t.test, mu0=10.2) < alpha)
sum(apply(M, 1, my.t.test, mu0=12) < alpha)



# Q: Maak een functie die voor 2 datasets zelf bepaald of je een 
# 2-sample t-toets of een Welch t-toets moet gebruiken:

# Datasets:
y.1 <- c(1.1, 1.3, 1.4, 1.1, 0.9)
y.2 <- c(1.5, 1.8, 0.7, 1.9, 1.1, 0.4)
# of
y <- c(y.1, y.2)
g <- factor(c(rep(1, length(y.1)), rep(2, length(y.2))))
var.test(y.1, y.2)
var.test(y ~ g)
tapply(y, g, mean)   # tapply voert functie mean uit per level van g
























