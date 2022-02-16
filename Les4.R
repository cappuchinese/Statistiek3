# Statistiek 3 BIN (2015-2016)
# Les 4
#

# Loopjes in R

# MATRIX
( M <- matrix(1:12, nrow = 4, byrow = T) )

# Q: Bereken het gemiddelde van M per rij:



# Q: Bereken het gemiddelde van M per kolom:




# VECTOR
( v <- 1:12 )

# Q: Bereken het gemiddelde van v:



# Q: Vervang in v alle elementen die groter zijn dan 6 door het getal 21




# Q: Bereken van elk element van v het kwadraat:



# Q: Bereken van elk element van v de wortel:



# Q: Schrijf een functie 'pow' die van een getal een willekeurige macht n neemt:



# Q: Pas deze functie 'pow'toe op v met macht 4:



# DATAFRAME

v <- c(4, 5, 4, 8, 6, 5, 5, 9, 6, 7, 6, 9)
f <- factor(rep(c("male", "female"), 6), levels = c("female", "male"))
( myData <- data.frame(score = v, gender = f) )

# Q: Maak een boxplot van score als functie van gender:



# Q: Bereken per geslacht de gemiddelde score:



# Q: Geef een samenvatting van de inhoud van myData:



# Q: Selecteer uit dataframe myData alle regels voor "female"


# Q: Sorteer het dataframe myData op grootte van score (laag naar hoog)


# Q: Voor welke regels uit myData geldt dat gender = male?



# Q: Verwijder alle mannen uit myData:



# LIST

# Een list is een datastructuur waar je verschillende datatypes (zoals
# getallen en strings) in kunt stoppen. Dit is handig als uitvoer van
# verschillende statistische toetsen in R (zoals t-toetsen), waarbij
# je zowel getallen als uitvoer hebt (zoals gemiddelde, p-waarde, etc.) maar ook
# strings (zoals naam van de toets, naam van de dataset, etc.).

L <- list()
L[[1]] <- 3.14
L[[2]] <- "Pythagoras"
L[[3]] <- c(1.1, 2.2, 3.3, 4.4)
L

L[[2]]
L[[3]]
L[[3]][4]

# OF

L <- list()
L$pi <- 3.14
L$wiskundige <- "Pythagoras"
L$waarde <- c(1.1, 2.2, 3.3, 4.4)
L

L$wiskundige
L$waarde
L$waarde[4]

L[[3]]   # Hoe werkt dit?

# OF

L <- list(pi = 3.14,
          wiskundige = "Pythagoras",
          waarde = c(1.1, 2.2, 3.3, 4.4))
L

L$wiskundige
L$waarde
L$waarde[4]

# Om te toetsen of de gemiddelde score van vrouwen verschilt van die van mannen,
# kunnen we een t-toets uitvoeren:

plot(myData$score ~ myData$gender)
t.test(myData$score ~ myData$gender)

# Hoe ziet de structuur van de uitvoer eruit?

results <- t.test(myData$score ~ myData$gender)
str(results)

# Dit is een LIST met verschillende elementen, zoals
results$statistic   # de t-waarde
results$p.value     # de p-waarde
results$estimate    # de gemiddelden van vrouwen en mannen
results$estimate[1] # het gemiddelde van vrouwen
