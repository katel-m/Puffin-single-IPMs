library(coda)
library(matrixStats)
library(reshape)
library(bayesplot)
library(reshape)
library(ggplot2)

# Load rates
load(file="outputs/WP4/Puffin IPMs/lunde_iom_ipm_demcorr_postBS_1302.Rdata")
str(mod_out)

mod_out$mean$cor.sa.rho
## Set sample and year number
nosamples <- 10000
noyears <- Tmax-2

###############
#### SETUP ####
###############

## Prepare matrices to rearrange samples - Vital rates & population sizes

# Non-varying vital rates
si <- rep(NA,nosamples)

# Time-varying vital rates, population sizes and growth rates
sa <- rho <- Na <- Na2 <- lambda <- matrix(NA, nrow = nosamples, ncol = noyears)

## Fill posterior samples into vectors and matrices
for(i in 1:nosamples){
 for(t in 1:noyears){
  
  # Time-invariant vital rates 
  si[i] <- mod_out$sims.list$phicomb[i]
  
  # Time-varying vital rates
  rho[i,t] <- mod_out$sims.list$rho[i,t]
  sa[i,t] <- mod_out$sims.list$sa[i,t]
  
  # Time-varying population sizes
  Na[i,t] <- mod_out$sims.list$N[i,t]
  Na2[i,t] <- 1
 }
 
 for(t in 1:(noyears-1)){
  # Population growth rate
  lambda[i,t] <- Na[i,t+1]/Na[i,t]
 }
}

## Make time-average population sizes
Na_mean <- rowMeans(Na[,1:(noyears-1)])

## Make time-average vital rates
si_mean <- si # already time-invariant
sa_mean <- rowMeans(sa, na.rm = T)
rho_mean <- rowMeans(rho, na.rm = T) # need to double check indexing is correct [,2:noyears]

# NOTE: I am averaging across indeces 1:22 for juvenile & adult survival/mortality (and population sizes), 
# and 2:23 for all reproductive parameters (and immigration) because in the matrix, the former appear with index t and the latter with index t+1

## Make time-average population growth rate
lambda_mean <- rowMeans(lambda, na.rm = T)

################################################
#### CALCULATION OF TRANSIENT SENSITIVITIES ####
################################################

## Calculate transient sensitivities for vital rates and population size/structure (evaluated at the temporal mean)
sens_si <- rep(NA, nosamples)
sens_sa <- rep(NA, nosamples)
sens_rho <- rep(NA, nosamples)
sens_Na <- rep(NA, nosamples)

## Get partial derivatives
gr.n <- expression((si*sa*rho*Na + sa*Na)/(Na))
sens.rho <- D(gr.n,'rho')
sens.si <- D(gr.n,'si')
sens.sa <- D(gr.n,'sa')
sens.Na <- D(gr.n,'Na')

for(i in 1:nosamples){
 
 sens_si[i] <- sa_mean[i] * rho_mean[i] * Na_mean[i]/Na_mean[i]
 
 sens_sa[i] <- (si_mean[i] * rho_mean[i] * Na_mean[i] + Na_mean[i])/(Na_mean[i])
 
 sens_rho[i] <- si_mean[i] * sa_mean[i] * Na_mean[i]/(Na_mean[i])
 
 sens_Na[i] <- (si_mean[i] * sa_mean[i] * rho_mean[i] + sa_mean[i])/(Na_mean[i]) - (si_mean[i] * sa_mean[i] * rho_mean[i] * Na_mean[i] + sa_mean[i] * Na_mean[i])/(Na_mean[i])^2
}

###############################################
#### CALCULATION OF TRANSIENT ELASTICITIES ####
###############################################

## Calculate transient elasiticities for vital rates and population size/structure (evaluated at the temporal mean)

elas_si <- sens_si*(si_mean/lambda_mean)
elas_sa <- sens_sa*(sa_mean/lambda_mean)
elas_rho <- sens_rho*(rho_mean/lambda_mean)
elas_Na <- sens_Na*(Na_mean/lambda_mean)

###############################################################
#### CALCULATE LTRE CONTRIBUTIONS ####
###############################################################

## Prepare vectors to store results
cont_sa <- rep(NA, nosamples)
cont_rho <- rep(NA, nosamples)
cont_sa_rho <- rep(NA, nosamples)
cont_Na <- rep(NA, nosamples)
cont_tot <- rep(NA, nosamples)
est_var <- matrix(NA, nrow = 3, ncol = nosamples)
est_covar <- matrix(NA, nrow = 3, ncol = nosamples)

## Calculate LTRE contributions (random design LTRE)

for(i in 1:nosamples){
 
 ## Make matrix of vital rates and population sizes
 dp_stoch <- cbind(sa[i,1:(noyears)],
                   rho[i,1:noyears],
                   Na[i,1:(noyears)]
 )
 
 ## Derive process variances and covariances
 dp_varcov <- var(dp_stoch)
 
 ## Save total estimated (co)variance per parameter
 est_var[,i] <- diag(dp_varcov)
 est_covar[,i] <- rowSums(dp_varcov)
 
 ## Make a vector of sensitivities
 sensvec <- c(sens_sa[i],
              sens_rho[i],
              sens_Na[i])
 
 ## Calculate demographic contributions
 # NOTE: Here we multiply sensitivities and (co)variances
 
 cont.mat <- matrix(NA, nrow = length(sensvec), ncol = length(sensvec))
 for(k in 1:length(sensvec)){
  for(m in 1:length(sensvec)){
   cont.mat[k,m] <- dp_varcov[k,m]*sensvec[k]*sensvec[m]
  }
 }
 rownames(cont.mat) = colnames(cont.mat) = c("S","F","N")
 ## Summarise contributions (sum of variances and covariances)
 cont_sa[i] <- cont.mat[1,1]
 cont_rho[i] <- cont.mat[2,2]
 cont_sa_rho[i] <- cont.mat[2,1]
 cont_Na[i] <- cont.mat[3,3]
 cont_tot[i] <- sum(cont.mat)
}

tot_cont=matrix(c(cont_rho,cont_sa_rho,cont_sa,cont_tot),nrow = 10000,byrow = F)
colnames(tot_cont)=c("rho","sa_rho","sa","sum")
tot_cont_prop=tot_cont
tot_cont_prop[,1]=tot_cont_prop[,1]/tot_cont_prop[,4]*100
tot_cont_prop[,2]=tot_cont_prop[,2]/tot_cont_prop[,4]*100
tot_cont_prop[,3]=tot_cont_prop[,3]/tot_cont_prop[,4]*100

res_df=data.frame(par=c("F","rho_sa","sa"),
                  mean=c(mean(tot_cont_prop[,"rho"]), mean(tot_cont_prop[,"sa_rho"]), mean(tot_cont_prop[,"sa"])),
                  lcr=c(quantile(tot_cont_prop[,"rho"],0.05),quantile(tot_cont_prop[,"sa_rho"],0.05),quantile(tot_cont_prop[,"sa"],0.05)),
                  ucr=c(quantile(tot_cont_prop[,"rho"],0.95),quantile(tot_cont_prop[,"sa_rho"],0.95),quantile(tot_cont_prop[,"sa"],0.95)))


# Extra columns for merging
res_df$col="iom"
res_df$lag="post"

p1=ggplot(res_df, aes(x = par, y = mean)) + 
 geom_bar(stat = "identity", width = 0.7,show.legend = T,colour="black",position="dodge")  +
 geom_errorbar(aes(ymin=lcr, ymax=ucr), width=.2,
               position=position_dodge(.9)) +
 xlab(" ")+
 ylab(expression(paste("% contribution to variance of log(", italic(lambda)[t], ")"))) +
 theme_bw() + 
 theme(panel.border = element_blank(), 
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       strip.background = element_blank(),
       axis.line = element_line(colour = "black"),
       axis.text=element_text(family="arial",colour="black",size=14),
       axis.title =element_text(family="arial",colour="black",size=14),
       strip.text = element_text(family="arial",size=14))


p1

write.csv(res_df,file = "outputs/WP4/Puffin IPMs/ltre_dat_iom_postBS_1402.csv")

