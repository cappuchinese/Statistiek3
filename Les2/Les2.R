#################################################
#
# Statistiek 3 BIN
# Les 2
#
# Sort, Dataframes, Lists en Functions in R
#
#################################################


#################################################
# Sort
#################################################

a <- c(0, 5, 3, 4, 7, 10)

# Q: Welke elementen van a zijn groter dan 4?

# Q: Vervang alle elementen in a die groter zijn dan 4 door 0

# Q: Vervang elk element in a dat 0 is door 20

# Q: Sorteer de nieuwe vector op grootte

# Q: Vervang elk element in a dat 20 is door NA (d.w.z. "not available")

# Q: Sorteer de nieuwe vector opnieuw



a <- c(0, 5, 3, 4, 7, 10)

# Q: Geef de eerste 3 elementen van de gesorteerde vector




#################################################
# Dataframes
#################################################

geneNames <- c("P53", "BRCA1", "VAMP1")
sig <- c(T, T, F)
meanExp <- c(4.5, 7.3, 5.4)

genes <- data.frame(geneNames, sig, meanExp)
rm(geneNames, sig, meanExp)

genes
geneNames
sig
meanExp

genes[1,]
genes[,1]
genes[,"geneNames"]
genes$geneNames

attach(genes)
geneNames
sig

geneNames[1:2]

detach(genes)
geneNames

geneNames <- c("P53", "BRCA1", "VAMP1")
sig <- c(T, T, F)
meanExp <- c(4.5, 7.3, 5.4)

# Je kunt ook de "variabelenamen" binnen een dataframe veranderen:

genes <- data.frame("names"= geneNames, 
                    "significance" = sig, 
                    "meanExp" = meanExp)
genes$names

genes[1:2,]   # eerste 2 rijen
genes[,1:2]   # eerste 2 kolommen
genes[1:2]    # OOK de eerste 2 kolommen!

str(genes)

# Sorteren op meerdere levels binnen dataframe

geneNames <- c("P53", "BRCA1", "VAMP1", "BRCA1", "VAMP1", "VAMP1")
sig <- c(T, T, F, T, F, F)
meanExp <- c(4.5, 7.3, 5.4, 5.5, 3.5, 7.9)
genes <- data.frame("names"= geneNames, 
                    "significance" = sig, 
                    "meanExp" = meanExp)
genes
genes[order(genes$names),]   # sorteer dataframe op gen-namen
genes[order(names),]         # waarom werkt dit niet?
genes[order(genes$names, genes$meanExp),]   # sorteer dataframe op naam, dan op expressie

search()        # wat is het zoekpad voor variabelen?


#################################################
# Lists
#################################################

geneName <- "P53"
annotations <- c("Onco gene", "Cell Cycle", "Apoptosis")
gene <- list("name" = geneName, 
             "annotation" = annotations)
gene
str(gene)

gene[[1]]    # 1e element uit list = gen-naam
gene[[2]]    # 2e element uit list = vector met annotaties
gene[[3]]    # 3e element uit list, bestaat niet!

gene[[2]][3] # 3e element uit vector van annotaties, die zelf het 2e element uit de list is 
gene$name
gene$annotation
gene$annotation[3]

gene[1]      # dit is zelf weer een list, met alleen 1 element, nl. name
gene[2]      # ook dit is zelf een list, met 1 element, nl. annotation

#################################################
# Functions
#################################################

addFour <- function(x){
  x <- x + 4
  return(x)
}

a <- c(0, 5, 3, 4, 7, 10)
addFour(a)
b <- addFour(a)
plot(a,b)

# Gebruik van sapply op vector a:

b <- sapply(a, addFour)
b

# Gebruik van apply op matrix M:

M <- matrix(1:12, nrow=4, ncol=3, byrow=T)    # vouw de vector 1:12 per rij tot een matrix!
M

apply(M, 1, mean)  # vector met gemiddelden per rij (= dimensie 1)
apply(M, 2, mean)  # vector met gemiddelden per kolom (= dimensie 2)

apply(M, 1, addFour)      # Wat gebeurt hier?
apply(M, 2, addFour)
apply(M, 1:2, addFour)
M + 4


