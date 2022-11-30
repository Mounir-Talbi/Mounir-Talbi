#R code pour le TP1 - Simulation, Inférence MCMC - M1 Maths UFR MIM

# This version: 05/11/2022

#Exercice 1

#1
unif_utile <- function(n,a,b){
  u<-runif(n,0,1)
  u <- (b-a)*u+a
  return(u)
}

#2 
a <-1
b <-3
c <- 5
d <-7
n<- 10000

x <- matrix(rep(0,2*n),n,2)

for (i in 1:n){
  x[i,1] <- (b-a)*runif(1)+a
  x[i,2] <- (d-c)*runif(1)+c
}
plot(x[,1],x[,2])

#3 - Il y avait ici deux types de réponses possibles, résumées dans les deux versions ci-dessous : (On peut également utiliser directement la fonction sample de R)
#Version 1
quantile_unif_discrete <- function(t,N){
  i <- 1
  while(t>i/N){
    i=i+1
  }
  return(i)
}

#Version 2
quantile_unif_discreteBis <- function(t,N){
  return(floor(t*N)+1)
}



#Exercice 2
#Après inversion des fonctions de répartition des mesures de probabilités données dans l'énoncé, on obtient :

#1
simulExpon <- function(n,a){
  -log(runif(n))/a
}

#2
simulCauchy <- function(n,a){
  a*tan(pi*(runif(n)-0.5))
}

#3
simulPareto <- function(n,a){
  1/runif(n)^(1/a)
}


#Exercice 3

#2
n=10000

I1 <- function(n){
  U <- runif(n)
  I_1 <- sum(4*sqrt(1-U^2))/n
  return(I_1)
}

print(I1(n))

I2 <- function(n){
  U <- runif(n,-1,1)
  V <- runif(n,-1,1)
  I_2 <- 4*sum(U^2+V^2<=1)/n
  return(I_2)
}

print(I2(n))

I3 <- function(n){
  U <- runif(n,-1,1)
  V <- runif(n,-1,1)
  W <- runif(n,-1,1)
  I_3 <- 2*3*sum(U^2+V^2+W^2<=1)/n
  return(I_3)
}

print(I3(n))

#3
N <- floor(4*10^4*qnorm(0.975)^2*5/3)+1
print(N)

#4
#Nous pouvons illustrer la convergence, en faisant par exemple varier le nombre n de simulations :

n_vec <- seq(100,100000,500)

#En utilisant la fonction R vapply, nous allons appliquer les fonctions I1, I2 et I3 le long du vecteur n_vec :

I1_vec <- vapply(n_vec, I1, FUN.VALUE = sqrt(2))
I2_vec <- vapply(n_vec, I2, FUN.VALUE = sqrt(2))
I3_vec <- vapply(n_vec, I3, FUN.VALUE = sqrt(2))
#Le dernier paramètre sqrt(2) sert simplement à indiquer que les fonctions I1, I2 et I3 renvoient des réels.

#Nous pouvons maintenant faire les graphiques 
plot(n_vec, I1_vec, type = "l", xlab ="", ylab = "", col="red", lwd=1.5)
lines(n_vec, I2_vec, type = "l", xlab ="", ylab = "", col="blue", lwd=0.7)
lines(n_vec, I3_vec, type = "l", xlab ="", ylab = "", col="green", lwd=0.3)
#On constate graphiquement que la méthode I1 (de plus faible dimension) semble converger plus rapidement


#Exercice 4

#Nombre de simulations
n<-10000

#Données du problème
m<-0.6947
sigma<-0.6563
P<-166275000

# Question 1. 
X <- rnorm(n,m,sigma)
Z <- P*exp(X)
S <- seq(500000000,2000000000,10000)
MC_n <- c()
exact <- c()
i=1
for (K in S){
  MC_n[i] = 1/n*sum(Z>K)
  Ktilde = log(K/P)
  exact[i]=1-pnorm(Ktilde,m,sigma)
  i=i+1
}

# Question 2.
#Dans la boucle ci-dessus, on a calculé à la fois les valeurs estimées par Monte-Carlo de la 
#probabilité qui nous intéresse, ainsi que la valeur "exacte" (approchée numériquement de manière
# précise par R) de cette même probabilité.
#Nous pouvons maintenant représenter ces deux courbes (en fonction de K) sur un même graphique :
plot(S, MC_n, type = "l", xlab ="", ylab = "", col="red", lwd=1.5)
lines(S, exact, type = "l", xlab ="", ylab = "", col="blue", lwd=1.5)

#On constate que pour les grandes valeurs de K, l'approximation de Monte-Carlo donne une valeur
#nulle, alors que la courbe bleue est strictement positive.

# Question 3.
#La densité g donnée dans l'énoncé correspond à une gaussienne centrée en log(K/P). Nous aurons donc
#environ la moitié des valeurs simulées qui seront au dessus de log(K/P), même si K est grand.

# Question 4.
S <- seq(2000000000,6000000000,50000000) #de 2*10^9 à 6*10^9 par pas de 5*10^7

Imp_sampling <- function(m,sigma,P,n){
  i=1
  alpha_g <- c()
  exact <- c()
  for (K in S){
    Ktilde <- log(K/P)
    Y <- rnorm(n,Ktilde,sigma) 
    alpha_g[i] <- 1/n*sum(exp(-1/(2*sigma^2)*(2*Y-m-Ktilde)*(Ktilde-m))*(Y>Ktilde))
    exact[i] <- 1-pnorm(Ktilde,m,sigma)
    i=i+1
  } 
  A_renvoyer <- list("Approx" = alpha_g, "Exact" = exact)
  return(A_renvoyer)
}

#graphe de la valeur "exacte" VS valeur approximée par échantillonage préférentiel
Imp <- Imp_sampling(0.6947,0.65,166275000,1000)
plot(S,Imp$Approx, type = "l", xlab ="", ylab = "", col="red", lwd=1.5)
lines(S,Imp$Exact, type = "l", xlab ="", ylab = "", col="blue", lwd=1.5)





