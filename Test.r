library(adegraphics)
library(ade4)

##### On construit le premier tableau

# C'est une simulation des abondances d'espèces dans plusieurs stations
# Ces stations sont au nombre de 16
# Elles sont réparties régulièrement le long d'un gradient de température
# On suppose que les abondances des espèces sont données par le modèle dit des températures cardinales
# Il y a en tout 10 espèces observées

# Fonction d'abondance des poissons
CTMI <- function(T, param)
{
  Tmin <- param[1]
  Topt <- param[2]
  Tmax <- param[3]
  Muopt <- param[4]
  if( T <= Tmin || T >= Tmax )
  # Si on est en dessous de la température min ou au dessus de la température max, le poisson est absent
  {
    return(0)
  }
  else
  # Si on est entre Tmin et Tmax
  {
      Num <- (T-Tmax)*(T-Tmin)^2
      Den <- (Topt-Tmin)*((Topt-Tmin)*(T-Topt)-(Topt-Tmax)*(Topt+Tmin-2*T))
      return(Muopt*Num/Den)
  }
}

# 16 stations, avec leurs températures moyennes
Tstations <- seq(from = 5, to = 80, by = 5)
Tstations
length(Tstations)

nsp <- 10 # Nombre d'espèces
topt <- seq(from = 5, to = 100, length = nsp) # Températures optimales, une par espèce
tmin <- 0.953*topt - 28.913 # Température minimale, une par espèce
tmax <- 1.101*topt + 3.203 # Température maximale, une par espèce

# On construit le tableau :
# 10 lignes correspondant aux 10 espèces, 16 colonnes correspndant aux 16 stations
data <- data.frame(matrix(nrow = nsp, ncol = length(Tstations)))
colnames(data) <- paste("st",1:length(Tstations), sep = "")
rownames(data) <- paste("sp", round(topt, 0), sep = "")

for( i in 1:nsp )
{
  for( j in 1:length(Tstations))
  {
    # Pour chaque case de data, on simule la valeur à partir de la température de la station, de la
    # température optimale de l'espèce, de la température maximale, de la température minimale, de Muopt,
    # et de la fonction CTMI
    data[i, j] <- round(CTMI(T = Tstations[j], param = c(tmin[i], topt[i],tmax[i], 1000)), 0)
  }
}

# On va visualiser tout ça
Taxis <- 1:100
cols <- rev(rainbow(nsp, end = 5/6))
plot(x = Taxis, y = sapply(Taxis, CTMI, param = c(tmin[1], topt[1], tmax[1], 1000)), type = "l",
     col = cols[1], main = paste("Répartition des", nsp, "espèces"), xlab = "Température [C]",
     ylab = "Nombre d'individus", lwd = 2)
for( i in 2:nsp )
{
  lines(x = Taxis, y = sapply(Taxis, CTMI, param = c(tmin[i], topt[i], tmax[i], 1000)),
        col = cols[i], lwd = 2)
}

# La structure des données est très simple, on a une succession des espèces le long du gradient
table.value(data, clegend = 0)
# En abscisse les stations, en ordonnée les espèces

# Dans la vraie vie, on ne connait pas dès le début l'ordre des espèces et des stations qui mettent
# en évidence cette structure des données
# Elles ressembleraient plutôt à ça :
data.vv <- data[sample(nrow(data)), sample(ncol(data))] # On permute aléatoirement lignes et colonnes
table.value(data.vv, clegend = 0)

# On fait l'AFC (car données d'abondance) du tableau
afc <- dudi.coa(data.vv, scann = FALSE)
scatter(afc)
# On visualise bien l'effet Guttman (les données ont une forme de parabole)
# Cet effet témoigne de liaisons fortes

# Cela permet de ré-ordonner le tableau pour mettre en évidence le gradient
data.vv.ord <- data.vv[order(afc$li[, 1]), order(afc$co[, 1])]
table.value(data.vv.ord, clegend = 0)

##### On construit le deuxième tableau

# Mesure de la température dans chaque station, pour chaque heure de la journée
mesureT <- sapply(Tstations, function(x){10*sin(seq(0, pi, le = 24)) + x})
# Donc tableau avec 24 lignes et 16 colonnes
colnames(mesureT) <- paste("st",1:16,sep="")
rownames(mesureT)<-paste("h",1:24,sep="")
plot(0, type = "n", xlim = c(1,24), ylim = range(mesureT), xlab = "Heures", ylab = "Temperature",
     las = 1, main = "Température des stations") # Cadre vide
apply(mesureT, 2, lines) # 2 donc application en colonne, donc on affiche une courbe par station

# On ajoute un bruit Gaussien, pour plus de réalisme
set.seed(1)
mesureT <- sapply(Tstations, function(x){10*sin(seq(0, pi,le = 24)) + x + rnorm(n = 24, sd = 5)})
colnames(mesureT) <- paste("st", 1:16, sep = "")
rownames(mesureT) <- paste("h", 1:24, sep = "")
plot(0, type = "n", xlim = c(1,24), ylim = range(mesureT), xlab = "Heures",ylab = "Temperature",
     las = 1, main = "Température des stations") # Cadre vide
apply(mesureT, 2, lines)
# On a donc un gradient des stations froides aux stations chaudes
# Et une composante circadienne
# C'est difficile à voir vu l'intensité du bruit
table.value(mesureT, clegend = 0)

# On fait l'ACP (car données quanti) du tableau
acp <- dudi.pca(t(mesureT), scan=FALSE) # t pour mettre les stations en ligne donc en tant qu'individus
scatter(acp)
# Le premier facteur représente la quasi-totalité de la variabilité, c'est simplement le gradient
# thermique entre les stations
plot(Tstations, acp$li[,1], ylab = "F1", xlab = "Température",
     main = "Le premier facteur de l'ACP est\nle gradient thermique", pch = 19)


##### Couplage des deux tableaux

# On fait les vérifications nécessaires
# 1. Les individus communs aux deux tableaux sont-ils bien en ligne dans les deux tableaux ?
dim(data)
dim(mesureT)
tab1 <- as.data.frame(t(data))
tab2 <- as.data.frame(t(mesureT))
dim(tab1)
dim(tab2)
# 2. Les individus communs aux deux tableaux sont-ils bien dans le même ordre dans les deux tableaux ?
all.equal(rownames(tab1),rownames(tab2))
# 3. La même pondération a-t-elle bien été utilisée pour les deux analyses ?
coa <- dudi.coa(df = tab1, scannf = FALSE, nf = 2)
pca <- dudi.pca(df = tab2, row.w = coa$lw, scannf = FALSE, nf = 2)
all.equal(pca$lw,coa$lw)

# Tout est bon, on peut utiliser l'analyse de co-inertie pour coupler les deux tableaux
cia <- coinertia(dudiX = pca, dudiY = coa, scannf = FALSE, nf = 2)
cia$eig[1]/sum(cia$eig)
# Le premier facteur extrait 99.3% de la variabilité, c'est le gradient thermique (commun aux 2 tableaux)

# Significativité de la co-structure entre les deux tables ?
cia$RV
# On fait un test basé sur la comparaison entre la valeur de la co-inertie totale (coefficient RV entre
# les 2 tables) et la distribution empirique du coefficient RV qd on détruit la co-structure entre les 2
# tables en permutant les lignes de l'une
ciatest <- randtest(cia, nrepet =  999, fixed = 2)
plot(ciatest)
# Donc ici co-structure très significative entre les 2 tables, car RV est bien supérieur à ce qu'il
# devrait être si c'était une distribution aléatoire (= si pas de corrélation) !

# Représentation graphique des trois tableaux analysés :
plot.new()
par(mfrow = c(2, 2))
# Données de température (tab2, ACP)
par(mfg = c(1, 2))
table.value(pca$tab, clabel.row = 0.7, clabel.col = 0.7)
# Données d'abondance (tab1, AFC)
par(mfg = c(2, 1))
table.value(t(coa$tab), clabel.row = 0.7, clabel.col = 0.7)
# Analyse de co-inertie
par(mfg = c(2, 2))
table.value(cia$tab, clabel.row = 0.7, clabel.col = 0.7)
