#------------------------------------------------------------
#
# Script to fit dynamic ZIP N-mixture model to helicopter counts
# of Tahr and calculate finite population totals for each 
# Management unit (MU) using 'hybrid' designed-based and model-
# based inference
#
# Dave Ramsey (7/3/2022)
#----------------------------------------------------------------
options(tidyverse.quiet = TRUE)
options(dplyr.summarise.inform = FALSE)
library(tidyverse)
library(sf)
library(nimble)
library(nngeo)
library(MCMCvis)
library(bayesplot)
library(gridExtra)
library(progress)
library(ggspatial)
library(ggvoronoi)

source("r/misc_functions.r")

# Data input ----

tahr_data<- readRDS("Data/tahr_data.rds")

plots<- readRDS("Data/plots.rds")     # spatial (sf object) of plot boundaries
MU<- readRDS("Data/MU.rds")           # spatial (sf object) of management unit boundaries 
PCL<- readRDS("Data/PCL.rds")         # Spatial (sf object) of Public Conservation land
nzsouth<- readRDS("Data/nzsouth.rds") # spatial (sf object) of south island of New Zealand

plotID<- tahr_data$PlotID

## Collate data for Nimble ----

y<- as.matrix(tahr_data %>% select(y1,y2,y3))

nsites<- dim(y)[1]
T<- dim(y)[2]
# grab time between successive counts (occasions). convert days to weeks
tt<- tahr_data %>% select(t1,t2) %>% mutate(t1 = as.numeric(t1/7), t2 = as.numeric(t2/7))  
tt<- as.matrix(tt)

mu<- as.numeric(tahr_data$MU) # Grab MU numeric levels 
nmu<- max(mu)

form<- as.formula(~ MU)

lam_mat<- model.matrix( ~MU, data=tahr_data)  # Get design matrix
nbeta<- ncol(lam_mat)

AA<- tahr_data$Area_3d

# Assemble data for Nimble 
constants<- list(nsites=nsites, nyears=T, Area=AA, nb=nbeta) 
Data <- list(y=y, tt=tt, X=lam_mat)

## DM ZIP survival/recruitment model----

code<- nimbleCode({
  for(j in 1:nsites) {
    Nbar[j]<- mean(N[j,1:nyears])
    Nbar2[j]<- N[j,1]
    Dbar[j]<- Nbar[j]/Area[j]
    log(eta[j]) <- inprod(X[j, 1:nb], beta[1:nb])
  }
  
  omega ~ dbeta(1, 1)
  gamma ~ dgamma(0.1, 0.1)
  theta ~ dbeta(1, 1)
  p ~ dbeta(1, 1)
  
  for(j in 1:nb) {
    beta[j] ~ dnorm(0, sd=5)
  }
  
  for(i in 1:nsites) {
    z[i] ~ dbern(theta)
    mu[i,1]<- eta[i] 
    N[i,1] ~ dpois(mu[i,1] * z[i])
    D[i,1]<- N[i,1]/Area[i]
    y[i,1] ~ dbin(p, N[i,1])
    yrep[i,1] ~ dbin(p, N[i,1])
    for(t in 2:nyears) {
      S[i,t-1] ~ dbin(omega^tt[i,t-1], N[i, t-1])
      G[i,t-1] ~ dpois(N[i,t-1] * (gamma*tt[i,t-1])) 
      N[i,t] <- S[i,t-1] + G[i,t-1]
      y[i,t] ~ dbin(p, N[i,t])
      D[i,t]<- N[i,t]/Area[i]
      yrep[i,t] ~ dbin(p, N[i,t])
    }
  }
})


ymax<- apply(y, 1 ,max)
yinit<- ymax + 200
linit<- cbind(yinit,yinit,yinit)
zi<- ifelse(ymax > 0, 1, 0)

inits = function(){
  list(gamma=runif(1), omega=runif(1), eta = yinit, mu=linit, N=linit, beta=rnorm(nbeta),
              z=zi, theta=runif(1), p=runif(1))
}

parameters<- c("p","omega","gamma","beta","N","Nbar","D","Dbar","yrep","theta")

## create the model object
Rmodel <- nimbleModel(code=code, constants=constants, data=Data, inits=inits())
ModSpec <- configureMCMC(Rmodel)
ModSpec$resetMonitors()
ModSpec$addMonitors(parameters)
ModSpec$removeSamplers('beta', print=FALSE)
ModSpec$addSampler(target='beta', type='AF_slice')
Rmcmc <- buildMCMC(ModSpec)
Cmodel<- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)


# MCMC settings
ni <- 10000   # increase for serious inference
nt <- 1
nb <- 5000
nc <- 3

samp<- runMCMC(Cmcmc, niter=ni, nburnin=nb,thin=nt,nchains=nc,inits=inits,
               samplesAsCodaMCMC = TRUE)

## Summarise results ----

MCMCsummary(samp, params=c("beta","omega","gamma","theta","p"), n.eff=TRUE, round=3)



win.graph(10,7)
PR<- rbeta(1e4, 1, 1)
MCMCtrace(samp, params="p", priors=PR,pdf=F,post_zm=F)


# Print betas
tmp<- MCMCsummary(samp, params=c("beta"), n.eff=T, round=3)
dimnames(tmp)[[1]]<- dimnames(lam_mat)[[2]]
tmp


## Design-based predictions ----

# Find total number of potential plots within PCL in each MU
A<- 4  # set area of plots
ncells<- calc_nplots(plots, MU, PCL, res = sqrt(A) * 1000)  # res in m

# retrieve plot level abundance estimates
summ<- MCMCsummary(samp, params = c("Nbar"))
results<- tibble(Nhat=summ[,"mean"], varN=summ[,"sd"]^2)
results<- bind_cols(select(tahr_data, PlotID, MU), results)

est_list<- split(results, results$MU)
ncells_list<- split(ncells, ncells$MU)

tmp<- map2_dfr(est_list, ncells_list, calc_strata_nhat, A=A, type="abundance", .id="MU")

pred_n<- bind_rows(tmp, calc_total_nhat(tmp))

# round to nearest 50
pred_n<- pred_n %>% mutate(Nhat = round_any(Nhat,50), SD=round_any(SD,1),
                           lcl = round_any(lcl,50), ucl = round_any(ucl,50))
# Abundance estimates
pred_n

# Density estimates for each MU
pred_d<- map2_dfr(est_list, ncells_list, calc_strata_nhat, A=A, type="density", .id="MU")
pred_d


## Model based  Predictions ----

# Create df for prediction 
df<- split(ncells, ncells$MU) %>% map_dfr(function(x){tibble(MU=rep(x$MU,x$nplots))})
df<- df %>% mutate(MU = factor(MU))
pred_mat<- model.matrix(~MU, df)

pred_mu<- predict_model(samp, pred_mat, df, nsamp=5000)

pred_mu<- pred_mu %>% mutate(Nhat= round_any(Nhat,50),se=round_any(se,1),lcl=round_any(lcl, 50),
                             ucl=round_any(ucl, 50))

pred_mu


## Posterior predictive checks ----

pp<- MCMCpstr(samp, params = c("yrep"), type="chains")
pp<- t(apply(pp$yrep, c(1,3), sum))

yy<- apply(y, 1, sum)

prop_zero<- function(x) mean(x == 0) 

p1<- ppc_stat(yy, pp, stat=prop_zero, binwidth=0.01) + 
  ggtitle("Proportion zeros") + xlim(0.1, 0.4) + theme_bw() + 
  theme(legend.position = "none", plot.title = element_text(size=10))

p2<- ppc_stat(yy, pp, stat=mean) + ggtitle("Mean") + xlim(30, 40)  + 
  theme_bw() + theme(legend.position = "none", plot.title = element_text(size=10))

p3<- ppc_stat(yy, pp, stat=sd) + ggtitle("Standard deviation") + xlim(40, 60) + 
  theme_bw() + theme(legend.position = "none",plot.title = element_text(size=10))

p4<- ppc_stat(yy, pp, stat=max) + xlim(200, 350) + ggtitle("Maximum") + 
  theme_bw() + theme(legend.position = "none", plot.title = element_text(size=10))

win.graph(8,8)
grid.arrange(p1,p2,p3,p4, ncol=2)

# Figure 2
win.graph(8,8)
pp_check(sqrt(yy), sqrt(pp), fun="scatter_avg") + theme_bw()



## Plot densities ----

# Average density

summ<- MCMCsummary(samp, params = "Dbar")
results<- tibble(Density=summ[,"mean"], sd=summ[,"sd"], lcl=summ[,"2.5%"], ucl=summ[,"97.5%"])
results<- bind_cols(tahr_data, results)
results<- results %>% mutate(ytot=rowMeans(y)/Area_3d) 


win.graph(10,5)
results %>% mutate(Plot = fct_reorder(PlotID, Density, .desc = TRUE)) %>%
  ggplot(aes(x=Plot, y=Density))+
  geom_errorbar(aes(ymin=lcl,ymax=ucl), color="blue", width=0.15, size=1) +
  geom_point(colour="blue", size=3, shape=21, fill="white")+
  geom_point(aes(Plot,ytot), colour="red", size=1,shape=16) +
  scale_y_continuous(breaks=seq(0,40,10),lim=c(0, 40))+
  ylab(expression(paste("Density (tahr k",m^-2,")"))) +
  xlab("Plot") +
  theme_bw() +
  theme(text = element_text(size=15),
        axis.text.x = element_blank())

# now do density for each occasion

summ<- MCMCsummary(samp, params = "D")

allN<-tibble(
  Plot=rep(plotID,3),
  Count_num=rep(c(1,2,3),each=nsites),  
  D=as.vector(summ[,"mean"]),
  lcl=as.vector(summ[,"2.5%"]),
  ucl=as.vector(summ[,"97.5%"]),
  ytot=as.vector(y/AA))

win.graph(15,15)
allN %>% 
  ggplot(aes(x=Count_num, y=D)) +
  geom_line(colour="black", size=1.2)+
  geom_ribbon(aes(ymin=lcl,ymax=ucl, colour=NULL),alpha=0.25, fill="darkorange3")+
  geom_point(colour="blue", size=3, shape=21, fill="white")+
  ylab(expression(paste("Tahr k",m^-2,""))) +
  xlab("Sampling occasion")+
  #scale_y_continuous(breaks=seq(0,80,20),lim=c(0, 80))+
  scale_x_continuous(breaks=c(1,2,3),lim=c(1, 3))+
  geom_point(aes(x=Count_num, y=ytot), size=1.5) +
  facet_wrap(~Plot,  ncol=10, nrow=12, scales="free_y") +
  theme_bw()+
  theme(strip.background = element_blank(), 
        strip.text.x=element_text(hjust=0.05),
        panel.border = element_rect(colour = "black"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x=element_text(size=10))
        

## Density surface tessellation ----

summ<- MCMCsummary(samp, params = "Dbar")
results<- tibble(plotID=plotID, Density=summ[,"mean"])

xy<- as.data.frame(st_coordinates(st_centroid(MU)))
xy$MU<- c("MU1","MU2","MU3","MU4","MU5","MU6","MU7","EZ1","EZ2")

MU1<- MU %>% mutate(UNIT_NO=lvls_revalue(UNIT_NO,c("MU1","MU2","MU3","MU4","MU5","MU6","MU7","EZ1","EZ2")))

tmp<- plots
tmp<- as.data.frame(st_coordinates(st_centroid(tmp)))
tmp<- tmp %>% mutate(Density=results$Density)
MU1<- as_Spatial(MU1)

win.graph(10,10)
MU %>% ggplot() +
  geom_sf(fill=NA) +
  annotation_scale(location = "br", width_hint = 0.25) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  geom_voronoi(aes(X,Y,fill=Density),outline=MU1,data=tmp) +
  geom_sf(fill=NA, size=1, color="black") +
  geom_point(aes(X,Y, fill = Density), data=tmp) +
  scale_fill_distiller(palette = "OrRd", direction=1, limits=c(0,35)) +
    geom_label(aes(X, Y, label=MU), data=xy, size=3,label.padding=unit(0.1, "lines"),hjust="inward",
               vjust="inward") +
  labs(x="Longitude",y="Latitude",fill=expression(paste("Density (",km^2,")"))) +
  theme_bw() +
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))

