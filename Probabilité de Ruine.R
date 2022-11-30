#Paramètre de la fonction

u=25 #niveau de réserve de l'entreprise
M=2000 #Nombre d'assurés
m=0.5
sigma2=1
s=1
alpha=1
p1=0.3;p2=0.2;p3=0.5
lambda=matrix(rep(0,3),1)
lambda[1,1]=0.001; lambda[1,2]=0.005; lambda[1,3]=0.01 #Probabilité d'accident en fonction de la météo
H0=1 #Temps au jour précèdent l'année 2022
Q=matrix( c(0.35,0.6,0.05,0.1,0.7,0.2,0.1,0.4,0.5), 3, byrow = TRUE) #Matrice de transition météo
c=0.001 #Taux de la prime par semaine
nbSimul=1000 #Nombre de simulation
#______________________________________

#Constantes

T=52 #Nombre de semaines en 2022

#______________________________________

#Loi montant sinistre

MontantSinistre=function(M,m,sigma2,s,alpha,p1,p2,p3){
  
  Gauss=runif(M,m,sigma2)  #On simule la loi Gaussienne
  
  
  Z=exp(Gauss) #Loi Log-Normale
  
  
  #  On inverse la fonction qui simule la loi de Pareto, on obtient : s/x^(1/alpha)*1(x>=s)+1(x<s)
  
  Y=as.numeric(Gauss>=s)*(s/Gauss^(1/alpha))+as.numeric(Gauss<s) #Loi de Pareto
  
  p1*0+p2*Y+p3*Z #Loi du montant des sinistres
  
}

#_______________________________

#Chaine de MARKOV (Météo)

MarkovChain = function(Q, H0,T){
  res = 1:T
  res[1] = sample( x = 1:length(Q[1,]),size = 1, replace = TRUE, prob = Q[H0,] )
  if(T > 1){
    for(i in 2:T){
      res[i] = sample( x = 1:length(Q[1,]),size = 1, replace = TRUE, prob = Q[res[i-1],] )
    }
  }
  res
}

#_____________________________________

#Loi du nombres de sinistres

NombreSinistre =function(Q,H0,T,M){
  
  Meteo=matrix(MarkovChain(Q,1,T),1) 
  N=matrix(rep(0,T),1)
  for(i in 1:T){
    N[1,i]=lambda[1,Meteo[1,i]]*M
    
  }
  return(N)
}

#________________________________

# Simulation de la variation de l'argent en réserve durant une année


Simule=function(u,m,M,s,alpha,p1,p2,p3,lambda,H0,Q,c){

N=NombreSinistre(Q,H0,T,M)
plan=rep(c(0),T)#Pour garder les infos chaque semaine
  for (i in 1:T){
        cout=MontantSinistre(N[1,i],m,sigma2,s,alpha,p1,p2,p3)
      plan[i]=u-sum(cout)+c*M
  }
plot(plan,xlab="Semaine de l'année 2022",ylab="Argent en réserve",pch=4, col = "red")
lines(plan, col = "black") # Ligne de ruine
return(plan) #Argent qu'il reste à l'entreprise chaque semaine

}


#______________________________

# Faillite de l'entreprise durant une simulation ?

Ruine=function(plan,T){
  faillite=0 #Variable qui enregistre la faillite de la boite
  test=rep(FALSE,T)
  for (i in 1:T){
    test[i]=(plan[i]<=0)
    if(test[i] ==TRUE){
      faillite=1
    }
  }
  return(faillite)
}

#_______________________________

# Approximation de la probabilité de ruine

ProbaRuine=function(u,m,M,s,alpha,p1,p2,p3,lambda,H0,Q,c,nbSimul){
  
  nbRuine=0
  for (i in 1:nbSimul){
    plan=Simule(u,m,M,s,alpha,p1,p2,p3,lambda,H0,Q,c)
    faillite=Ruine(plan,T)
    nbRuine=nbRuine+faillite
  }
  proportion=nbRuine/nbSimul
    return(proportion)
}

ProbaRuine(u,m,M,s,alpha,p1,p2,p3,lambda,H0,Q,c,nbSimul)