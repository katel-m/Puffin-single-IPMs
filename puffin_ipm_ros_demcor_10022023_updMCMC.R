# Load packages 
#library(arm)
#library(Rfast)
library("coda")
#library("bayesplot")
#library("igraph")
library("jagsUI")
#library("devtools")
library("pracma")
library("matrixStats")

# Set seed
mySeed <- 0
set.seed(mySeed)

#***************#
# DATA ASSEMBLY #
#***************#
# Røst
m=read.csv("data/02_processed/WP4/data puffin ipms/puff_ros_marray.csv",sep=";",na.strings = "NA")
m=m[,-1]
row.names(m)=NULL
colnames(m)=NULL
rsum <- rowSums(as.matrix(m),na.rm=TRUE)

# Load breeding data 
breed_dat=read.csv("data/02_processed/WP4/data puffin ipms/puff_ros_breeding.csv",sep=";")
breed_dat=breed_dat[which(breed_dat$Year >= 1990),]
breed_dat=breed_dat[which(breed_dat$Year <= 2019),]

# Subset for correct correlation
#breed_dat=breed_dat[-1,] # surv-postBS
breed_dat=breed_dat[-30,] # preBS-surv

F.dat = breed_dat$nfled
E.dat = breed_dat$N

# Load population counts 
pop_dat=read.csv("data/02_processed/WP4/data puffin ipms/puff_ros_counts.csv",sep=";")
n=pop_dat[13:41,"N"]
n = round(n/100)
sd(n)
Tmax = dim(m)[2]+1

Data <- list(Tmax = Tmax, 
             d = 5,
             mean.Ntmp = n[1:5],
             m = as.matrix(m), # m-array
             rel = rsum, # m-array row sums
             Fdat = F.dat, # dataF = fledgling count
             Edat = E.dat, # E = breeder count
             n = round(n,0)) # n = population count

str(Data)

#************#
# MODEL CODE #
#************#
sink("dem_corr_ipm.txt")
cat("
    model{
  #------------------------------------------------------------------------
  # PRODUCTIVITY
  #------------------------------------------------------------------------
  
  # Priors ---------------------------------
  for (i in 1:(Tmax-2)){
    logit(rho[i]) <- b0.rho + eps[i,1]
  }
  b0.rho ~ dunif(-5,5)

  # Likelihood ------------------------------
  for (i in 1:(Tmax-2)){
    Fdat[i] ~ dbin(rho[i], Edat[i])
  }
  
  #------------------------------------------------------------------------
  # BREEDING ADULTS MARK-RESIGHT DATA
  #------------------------------------------------------------------------
  
  # Priors ----------------------------------
  for (i in 1:(Tmax-2)){
    logit(sa[i])<- b0.sa + eps[i,2]

    # Year-dependent component of p
    pt[i] ~ dnorm(0, 1.0E-4)
    
    # p when not captured the previous occasion
    logit(p[i]) <- pt[i]
    
    # p when captured the previous occasion
    logit(pdep[i]) <- pt[i] + a
  }
  
  a ~ dunif(-5, 5) # Prior for trap-dependence in p
  b0.sa ~ dunif(-5,5)
  
  # Multinomial likelihood --------------------
  for (i in 1:(Tmax-2)){
    m[i,1:(Tmax-1)] ~ dmulti(q[i,1:(Tmax-1)],rel[i])
  }
  
  # m-array cell probabilities ----------------
  chi[Tmax-1]<-1 # chi (recursion for never seen)
  
  for (j in 1:(Tmax-2)){
    chi[Tmax-1-j]<- 1 - sa[Tmax-1-j]*(1-(1-p[Tmax-1-j])*chi[Tmax-1-j+1])
  }
  
  # Cell probabilities
  for (i in 1:(Tmax-2)){
    # Cells in diagonal
    q[i,i]<-pdep[i]*sa[i]
    
    # Cells above diagonal
    for (j in (i+1):(Tmax-2)){
      logprod[i,j,i]<-log(sa[i]*(1-pdep[i]))
      for (k in (i+1):(j-1)){
        logprod[i,j,k]<-log(sa[k]*(1-p[k]))
      }
      q[i,j]<-p[j]*sa[j]*exp(sum(logprod[i,j,i:(j-1)]))
    }# Empty cells below diagonal
    for (j in 1:(i-1)){
      q[i,j]<-0}
    # Probability of an animal never seen
    q[i,Tmax-1] <- 1-sa[i]*(1-(1-pdep[i])*chi[i+1])
  }
  
   # Temporal random effects
  for (j in 1:2) {
    E[j, j] ~ dnorm(0.0, 1.0)T(0.0,)
    DeltaE[j, j] <- 1/tauE[j]
    tauE[j] ~ dgamma(1.5, 1.5)
    LE[j, j] <- 1.0
  }
  
  LE[1, 2] <- 0.0; E[1, 2] <- 0.0; DeltaE[1, 2] <- 0.0
  LE[2, 1] ~ dnorm(0.0, 4.0); E[2, 1] <- 0.0; DeltaE[2, 1] <- 0.0
  
  # covariance matrix
  Lambda <- E %*% LE %*% DeltaE %*% t(LE) %*% E
  
  for(i in 1:(Tmax-2)){
    eps[i, 1] <- E[1, 1] * (LE[1, 1] * xi_e[i, 1])
    eps[i, 2] <- E[2, 2] * (LE[2, 1] * xi_e[i, 1] + LE[2, 2] * xi_e[i, 2])
    
    for(j in 1:2){
      xi_e[i, j] ~ dnorm(0.0, tauE[j])
    }
  }
  
  # Derived quantities
  sigma.sa <- sqrt(Lambda[1, 1])
  sigma.rho <- sqrt(Lambda[2, 2])
  cor.sa.rho <- Lambda[1, 2] / sqrt(Lambda[1, 1] * Lambda[2, 2])
  
  #------------------------------------------------------------------------
  # POPULATION MODEL
  #------------------------------------------------------------------------

  # Priors ----------------------------------
  phicomb ~ dunif(0.3, 0.9)
  sigma.N ~ dunif(0, 1500) 
  tau.N <- pow(sigma.N,-2)
  
  # Initialisation of population model (time: 1 to d)----------
  for(i in 1:d){
    Ntmp[i] ~ dnorm(mean.Ntmp[i], tau.N)T(0.0, 100000)
  } 

  N[1:d] <- round(Ntmp[1:d])
  
  # Population model (time: d+1 to Tmax-1) ----------------------------------
  for (i in (d+1):(Tmax-2)){
    
    ## Survival
    S[i] ~ dbin(sa[i-1], N[i-1])
    
    ## Expected per-capita recruitment
    combinedPr[i] <- rho[i-d]*0.5*phicomb*sa[i-1]
    
    ## Realized recruitment (population-level)
    R[i] ~ dbin(combinedPr[i], N[i-d])
    
    ## Population size
    N[i] <- S[i] + R[i] 
  }
  
  # Likelihood for population census data:
  for (i in 1:(Tmax-2)){
    n[i] ~ dnorm(N[i], tau.N)
  }
  
    }
    ",fill = TRUE)
sink()

Inits <- function(){list(
  # survival
  b0.sa = rnorm(1, 0, 1),
  a = runif(1, -2, 2),
  pt = runif(Tmax-2, -1, 1),
  
  # productivity
  b0.rho = rnorm(1, 0, 1),
  
  # population counts 
  Ntmp = Data$mean.Ntmp
)}

parameters <- c('phicomb','N','sigma.N','tau.N','tauE',
                'b0.sa','sa','R', 'S', 'pt', 'a', 'b0.rho','rho',
                'eps', 'Lambda','cor.sa.rho','sigma.sa','sigma.rho')

testInits <- Inits() 

# One line running of model (test run) 
mod_out <- jags(model.file ='dem_corr_ipm.txt',
                data=Data,
                parameters.to.save = parameters,
                inits = Inits,
                n.chains=1, 
                n.iter = 1000, 
                n.burnin = 100, n.thin = 1,
                parallel=T)

# One line running of model (long run)
mod_out <- jags(model.file ='dem_corr_ipm.txt',
                data=Data,
                parameters.to.save = parameters,
                inits = Inits,
                n.chains = 3, 
                n.iter = 500000, 
                n.burnin = 50000, 
                n.thin = 10,
                n.adapt = 2000, 
                parallel=T)

View(mod_out$summary)

save(mod_out,file="outputs/WP4/Puffin IPMs/lunde_ros_ipm_demcorr_preBS_13022023.Rdata")


