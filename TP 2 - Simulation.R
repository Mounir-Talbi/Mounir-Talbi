#R code pour le TP2 - Simulation, Inférence MCMC - M1 Maths UFR MIM

# This version: 06/03/2022

#Exercice 1

#La seule différence entre les questions 1. et 2. ci-dessous est le graphique donné en output des fonctions.
#Dans le premier cas, on renvoie et on trace la suite de points.
#Dans le deuxième cas, on ajoute l'option type="l" à la fonction plot pour obtenir une interpolation linéaire

#Question 1. 

MarcheAleatoire1 <- function(nb_pas,P){
    pas <- c(1, -1)
    X<-c()
    X[1]<-0;
    for (k in seq(1,nb_pas)){
      proba <- c(P,1-P)
      X[k+1]<-X[k]+sample(pas,1,prob = proba)
    }
    plot(seq(0,nb_pas),X) #on commence bien l'axe des abcisses à 0 pour représenter le point (0,0)
}

MarcheAleatoire1(100,0.5)


#Question 2.

MarcheAleatoire1IL <- function(nb_pas,P){
  pas <- c(1, -1)
  X<-c()
  X[1]<-0;
  for (k in seq(1,nb_pas)){
    proba <- c(P,1-P)
    X[k+1]<-X[k]+sample(pas,1,prob = proba)
  }
  plot(seq(0,nb_pas),X, type="l") #on commence bien l'axe des abcisses à 0 pour représenter le point (0,0)
}

MarcheAleatoire1IL(10000,0.5) #Essayer avec P=0.48 et comparer!



#Question 3

#L'écart-type de la loi eta définie dans l'énoncé est égal à 1.
sigma<-1;

T<-1 #temps final
n<-1000
temps <- seq(0,T,0.001)

#On commence par simuler toute la trajectoire de X
pas <- c(1, -1)
P<-0.5
X<-c()
X[1]<-0;
for (k in seq(1,floor(n*T))){
  proba <- c(P,1-P)
  X[k+1]<-X[k]+sample(pas,1,prob = proba)
}

#On en déduit les valeurs de Y (qui est juste X, à un changement d'échelle près)
Y<-c()
for (t in temps){
    u=1+floor(n*t); #R n'acceptant pas l'indice 0, on décale tout de 1 (à la fois pour X et Y)
    Y <- c(Y,X[u]/(sigma*sqrt(n)))
}
plot(temps,Y, type="l")


#Question 4

#En reproduisant M fois les opérations de la question précédente, on obtient un échantillon de taille M de Y_T^n. On compare ensuite l'histogramme empirique des M valeurs simulées, avec celui de la loi normale centrée réduite.
M<-1000
Z<-c()
for (m in 1:M){
  X<-c()
  X[1]<-0;
  for (k in seq(1,floor(n*T))){
    proba <- c(P,1-P)
    X[k+1]<-X[k]+sample(pas,1,prob = proba)
  }
  Y<-c()
  for (t in temps){
    u=1+floor(n*t); #R n'acceptant pas l'indice 0, on décale tout de 1 (à la fois pour X et Y)
    Y <- c(Y,X[u]/(sigma*sqrt(n)))
  }
  Z[m]=Y[1+floor(n*T)];
}

#On représente maintenant graphiquement un histogramme des valeurs de Z
hist(Z, 30, freq = FALSE)

#On superpose ensuite la densité gaussienne standard
x<-seq(-4,4,length.out = 100)
lines(x,dnorm(x), type = "l", xlab ="", ylab = "", col="red", lwd=1.5)


#Exercice 3


MarkovChain <- function(Q,i,Npas){
  n <- nrow(Q) #cardinal de l'espace d'états
  X<-c()
  X[1]<-i
  for (j in 1:Npas){
    X[j+1] <- sample(1:n,1,prob = Q[X[j],])
  }
  return(X)
}

Q <- matrix(c(0.4,0.7,0,0,0.6,0.2,0.2,0,0,0.1,0.2,0.7,0,0,0.6,0.3),4,4)
MarkovChain(Q,4,30)


#Exercice 4

#"Au hasard" n'était pas clairement défini dans l'énoncé. Il fallait donc se donner une mesure de son choix pour tirer "au hasard" une matrice de transition.
# Dans ce qui suis, nous allons simuler (n-1) variables aléatoires indépendantes et de loi uniforme sur [0,1].
#Nous prenons ensuite la statistique d'ordre croissante associée à ces n-1 points.
#On obtient alors une subdivision aléatoire de [0,1], dont les longeurs des intervalles nous donneront une mesure de probabilité sur {1,...,n}.
#En répétant cette opération n fois, nous avons une matrice de transition.

RandomTransition <- function(n){
  Q <- matrix(rep(0,n*n),n,n)
  for (j in 1:n){
    U <- sort(runif(n-1))
    Q[j,] <- c(U[1], diff(U), 1-U[n-1])
  }
  return(Q)
}

#Simulation des 30 premiers pas d'une chaîne de Markov, de matrice de transition Q tirée au hasard (4 états), partant de l'état 1 :
MarkovChain(RandomTransition(4),1,30)


