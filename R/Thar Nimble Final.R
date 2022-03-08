
# Data input ----
options(tidyverse.quiet = TRUE)
options(dplyr.summarise.inform = FALSE)
library(tidyverse)
library(readxl)
library(lubridate)
library(sf)
library(terra)
library(forcats)
library(nngeo)
library(nimble)
library(MCMCvis)
library(bayesplot)
library(mcmcOutput)
library(ggspatial)
library(ggvoronoi)
library(viridis)
library(progress)
source("r/misc_functions.r")

ung16<- read_excel("data/Tier1_Tahr_2016.xlsx",sheet="Aerial_and_GPS")
ung17<- read_excel("data/Tier1_Tahr_2017.xlsx",sheet="Aerial_and_GPS")
ung18<- read_excel("data/Tier1_Tahr_2018.xlsx",sheet="Aerial_and_GPS")
ung19<- read_excel("data/Tier1_Tahr_2019.xlsx",sheet="Aerial_and_GPS")

ung16<- dplyr::select(ung16, PlotID,Count_Number,Plot_Date,AnimalCode,NumberOfAnimal,POINT_X,POINT_Y)
ung17<- dplyr::select(ung17, PlotID,Count_Number,Plot_Date,AnimalCode,NumberOfAnimal,POINT_X,POINT_Y)
ung18<- dplyr::select(ung18, PlotID,Count_Number,Plot_Date,AnimalCode,NumberOfAnimal,POINT_X,POINT_Y)
ung19<- dplyr::select(ung19, PlotID,Count_Number,Plot_Date,AnimalCode,NumberOfAnimal,POINT_X,POINT_Y)

ung16<- dplyr::mutate(ung16, Year=1)
ung17<- dplyr::mutate(ung17, Year=2)
ung18<- dplyr::mutate(ung18, Year=3)
ung19<- dplyr::mutate(ung19, Year=4) %>% filter(!(PlotID %in% c("AT133", "AZ124"))) # two plots not in MU

ung<- bind_rows(ung16, ung17, ung18, ung19)
survey.times<- ung %>% group_by(PlotID, Count_Number) %>% summarise(Date=min(Plot_Date),Year=min(Year)) 
survey.times<- group_by(survey.times, PlotID) %>% summarise(t1=nth(Date,2)-nth(Date,1),
                                                t2=nth(Date,3)-nth(Date,2), Year=min(Year))
survey.times<- survey.times %>% arrange(PlotID)
plots.in.units<- st_nn(st_centroid(plots), MU, sparse=T)
survey.times<- survey.times %>% mutate(MU = MU$UNIT_NO[unlist(plots.in.units)])

thar.codes<- c("TAM","TYM","TAF","TAU","TU","TY","TK")

thar<- filter(ung, AnimalCode %in% thar.codes) %>%
  dplyr::select(PlotID, Count_Number, AnimalCode, NumberOfAnimal, POINT_X, POINT_Y) %>%
  group_by(PlotID, Count_Number, POINT_X, POINT_Y) %>% summarize(total=sum(NumberOfAnimal))

# Number of plots in each MU
# Designate plot area A

A<- 4

### Assign each sampled plot to a MU and calculate number of potential plots within each MU
### including only plots within PCL

ncells<- calc_ncells(plots, MU, docsouth, res = sqrt(A) * 1000)  # res in m

# 3D plot area
# 3D plot area of sampled cells 

Area<- plot_area(plots, covar,A=A)

# Compile counts in plot area44

ytot<- compile_counts(thar, plots, maxcount=3, A=A)


## Assemble covariate data from each plot ----

habitat_raster<- readRDS("data/habitat_raster.rds")
tmp<- raster::extract(habitat_raster, plots)
tmp<- do.call('rbind',tmp)
tmp<- as.data.frame(tmp)

Xvar<- tmp %>% mutate(Forest = Indigenous.Forest + Broadleaved.Indigenous.Hardwoods,
                      Tussock = Tall.Tussock.Grassland,
                      Grass = Low.Producing.Grassland + Depleted.Grassland + High.Producing.Exotic.Grassland,
                      Herbfield = Alpine.Grass.Herbfield,
                      Shrubs = Sub.Alpine.Shrubland + Manuka.and.or.Kanuka+Matagouri.or.Grey.Scrub+Fernland,
                      Snow = Permanent.Snow.and.Ice,
                      Water = River + Lake.or.Pond,
                      Nonhab = Gravel.or.Rock + Landslide)

Xvar<- Xvar %>% dplyr::select(Forest,Tussock,Grass,Herbfield,Shrubs,Snow,Water,Nonhab,Elevation,slope)


## Construct data for predictions over entire MU ----

XX<- as.data.frame(habitat_raster)
XX<- XX %>% mutate(Forest = Indigenous.Forest + Broadleaved.Indigenous.Hardwoods,
                   Tussock = Tall.Tussock.Grassland,
                   Grass = Low.Producing.Grassland + Depleted.Grassland + High.Producing.Exotic.Grassland,
                   Herbfield = Alpine.Grass.Herbfield,
                   Shrubs = Sub.Alpine.Shrubland + Manuka.and.or.Kanuka + Matagouri.or.Grey.Scrub + Fernland,
                   Snow = Permanent.Snow.and.Ice,
                   Water = River + Lake.or.Pond,
                   Nonhab = Gravel.or.Rock + Landslide)

XX<- XX %>% dplyr::select(Forest,Tussock,Grass,Herbfield,Shrubs,Snow,Water,Nonhab,Elevation,slope)

# transform Elevation
XX<- XX %>% mutate(Elevation = (Elevation - mean(Xvar$Elevation))/sd(Xvar$Elevation), 
                   slope = (slope - mean(Xvar$slope))/sd(Xvar$slope),
                   cell=1:ncell(habitat_raster))
ext<- st_bbox(MU)
target_crs<- st_crs(MU)$proj4string
tmprast<- raster(xmn=ext[1],ymn=ext[2],xmx=ext[3],ymx=ext[4],resolution=sqrt(A)*1000, crs=target_crs) 
tmprast<- setValues(tmprast, 0)

pred_rast<- tmprast
pred_rast<- mask(pred_rast, MU)
pred_rast<- mask(pred_rast, docsouth)
cells<- extract(pred_rast, MU, cellnumbers=TRUE)

new_data<- list()
mu_id<- list()
n<- nrow(MU)

for(i in 1:n) {
  cellnums<- cells[[i]]
  inds<- complete.cases(cellnums)
  cellnums<- cellnums[inds,]
  ncells<- nrow(cellnums)
  tmp<- subset(XX, cell %in% cellnums[,"cell"]) 
  tmp$MU<- as.character(MU$UNIT_NO[i])
  new_data[[i]]<- tmp
  mu_id[[i]]<- cbind(cellnums[,"cell"],i)
}  


new_data<- do.call('rbind',new_data)
mu_id<- do.call('rbind', mu_id)

new_data<- as.data.frame(new_data)
new_data<- new_data %>% mutate(MU = factor(MU))
mu_id<- as.data.frame(mu_id)
names(mu_id)<- c("cell","MU")

mu_inds<- list()
for(i in 1:n)
  mu_inds[[i]]<- range(which(mu_id$MU == i))

mu_inds<- do.call(rbind,mu_inds)


## Covariate models ------------------------------------------------

nsites<- dim(ytot)[1]
T<- dim(ytot)[2]
t1<- as.numeric(survey.times$t1)
t2<- as.numeric(survey.times$t2)
tt<- as.matrix(cbind(t1,t2))/7  #weeks
mu<- as.numeric(survey.times$MU)
nmu<- max(mu)

myscale<- function(x) return((x - mean(x))/sd(x))
dat<- Xvar
dat<- dat %>% mutate_at(vars(Elevation,slope), .funs=myscale)
dat<- dat %>% mutate(MU = survey.times$MU, year = factor(survey.times$Year))

form<- as.formula(~ MU)

form<- as.formula(~ Tussock + Grass + Herbfield + Shrubs + Elevation + MU)

lam_mat<- model.matrix(form, data=dat)
nbeta<- ncol(lam_mat)

pred_mat<- model.matrix(form, data=new_data)

## -- DM survival/recruitment model----

code<- nimbleCode({
  for(j in 1:nsites) {
    Nbar[j]<- mean(N[j,1:nyears])
    Dbar[j]<- Nbar[j]/Area[j]
    log(eta[j]) <- inprod(X[j, 1:nb], beta[1:nb])
  }
  
  omega ~ dbeta(1, 1)
  gamma ~ dgamma(0.1, 0.1)
  p ~ dbeta(1, 1)
  
  for(j in 1:nb) {
    beta[j] ~ dnorm(0, sd=5)
  }
  
  for(i in 1:nsites) {
    mu[i,1]<- eta[i] 
    N[i,1] ~ dpois(mu[i,1])
    D[i,1]<- N[i,1]/Area[i]
    y[i,1] ~ dbin(p, N[i,1])
    yrep[i,1] ~ dbin(p, N[i,1])
    for(t in 2:nyears) {
      S[i,t-1] ~ dbin(omega^tt[i,t-1], N[i, t-1])
      G[i,t-1] ~ dpois(gamma * N[i,t-1] * tt[i,t-1]) 
      N[i,t] <- S[i,t-1] + G[i,t-1]
      y[i,t] ~ dbin(p, N[i,t])
      D[i,t]<- N[i,t]/Area[i]
      yrep[i,t] ~ dbin(p, N[i,t])
    }
  }
})

constants<- list(nsites=nsites, nyears=T, Area=Area$Area_3d, nb=nbeta) 
Data <- list(y=ytot, tt=tt, X=lam_mat)

ymax<- apply(ytot, 1 ,max)
yinit<- as.vector(apply(ytot,1,max) + 200)
linit<- cbind(yinit,yinit,yinit)
zi<- ifelse(ymax > 0, 1, 0)

inits1 = list(p=runif(1), gamma=1, omega=0.5, eta = yinit, mu=linit, N=linit, 
              beta=rnorm(nbeta))


parameters<- c("p","omega","gamma","beta","N","Nbar","D","Dbar","yrep")

## create the model object
Rmodel <- nimbleModel(code=code, constants=constants, data=Data, inits=inits1)
ModSpec <- configureMCMC(Rmodel, enableWAIC = FALSE)
ModSpec$resetMonitors()
ModSpec$addMonitors(parameters)
ModSpec$removeSamplers('beta', print=FALSE)
ModSpec$addSampler(target='beta', type='AF_slice')
Rmcmc <- buildMCMC(ModSpec)
Cmodel<- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)


# MCMC settings
ni <- 60000
nt <- 4
nb <- 20000
nc <- 3

inits = function(){inits1}

samp<- runMCMC(Cmcmc, niter=ni, nburnin=nb,thin=nt,nchains=nc,inits=inits,
               samplesAsCodaMCMC = TRUE, WAIC = TRUE)

Cmcmc$getWAIC()
samp<- samp$samples

MCMCsummary(samp, params=c("beta","omega","gamma","p"), n.eff=TRUE, round=3)

MCMCsummary(samp, params=c("p"), n.eff=TRUE, round=3)

saveRDS(samp,"papers/results/mu_nozip_surv.rds")

## -- DM ZIP survival/recruitment model----

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

constants<- list(nsites=nsites, nyears=T, Area=Area$Area_3d, nb=nbeta) 
Data <- list(y=ytot, tt=tt, X=lam_mat)

ymax<- apply(ytot, 1 ,max)
yinit<- as.vector(apply(ytot,1,max) + 200)
linit<- cbind(yinit,yinit,yinit)
zi<- ifelse(ymax > 0, 1, 0)

inits1 = list(gamma=1, omega=0.5, eta = yinit, mu=linit, N=linit, beta=rnorm(nbeta),
              z=zi, theta=runif(1), p=runif(1))

parameters<- c("p","omega","gamma","beta","N","Nbar","D","Dbar","yrep","theta","Nbar2")

## create the model object
Rmodel <- nimbleModel(code=code, constants=constants, data=Data, inits=inits1)
ModSpec <- configureMCMC(Rmodel, enableWAIC = TRUE)
ModSpec$resetMonitors()
ModSpec$addMonitors(parameters)
ModSpec$removeSamplers('beta', print=FALSE)
ModSpec$addSampler(target='beta', type='AF_slice')
Rmcmc <- buildMCMC(ModSpec)
Cmodel<- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)


# MCMC settings
ni <- 60000
nt <- 4
nb <- 20000
nc <- 3

inits = function(){inits1}

samp<- runMCMC(Cmcmc, niter=ni, nburnin=nb,thin=nt,nchains=nc,inits=inits,
               samplesAsCodaMCMC = TRUE, WAIC = TRUE)

Cmcmc$getWAIC()
samp<- samp$samples

MCMCsummary(samp, params=c("beta","omega","gamma","theta","p"), n.eff=TRUE, round=3)



saveRDS(samp,"papers/results/mu_zip_surv.rds")

# Summarise outputs ----

samp<- readRDS("papers/results/mu_zip_surv.rds")


win.graph(10,7)
PR<- rbeta(1e4, 1, 1)
MCMCtrace(samp, params="p", priors=PR,pdf=F,post_zm=F)


# Print betas
tmp<- MCMCsummary(samp, params=c("beta"), n.eff=T, round=3)
dimnames(tmp)[[1]]<- dimnames(lam_mat)[[2]]
tmp

MCMCsummary(samp, params=c("beta","omega","gamma","theta","p"), n.eff=TRUE, round=3)

write.csv(tmp, file="papers/betas.csv",row.names=TRUE)


# PPD
pp<- MCMCpstr(samp, params = c("yrep"), type="chains")
pp<- rbind(pp$yrep[,1,], pp$yrep[,2,], pp$yrep[,3,])
pp<- t(pp)
yy<- c(ytot[,1], ytot[,2], ytot[,3])

# Some more
tiff("papers/figure2.tif", width=6, height=6, units="in",res=300, compression="lzw")
win.graph(6,6)
pp_check(sqrt(yy), sqrt(pp), fun="scatter_avg") + theme_bw()
dev.off()

bayes_R2_res(yy, pp)
## Design-based predictions ----


A<- 4  # set area of plots
ncells<- calc_ncells(plots, MU, docsouth, res = sqrt(A) * 1000)  # res in m

summ<- MCMCsummary(samp, params = c("Nbar"))
results<- tibble(Nhat=summ[,"mean"], varN=summ[,"sd"]^2)
results<- bind_cols(survey.times, results)

tmp<- calc_nhat2(results, ncells, A=A, CI.level=0.95)

pred_n<- tmp$thar.mu %>% dplyr::select(Nhat,SD,CV,lcl,ucl,n=U)
Ntot<- tmp$Ntot %>% mutate(n=sum(ncells$ncells))

pred_n<- bind_rows(pred_n, Ntot) %>% mutate(Nhat= round_any(Nhat,50),SD=round_any(SD,1),
                                            lcl = round_any(lcl,50),
                                            ucl = round_any(ucl,50))
pred_n

pred_d<- tmp$thar.D

write_csv(pred_n, file="papers/pred_n.csv")
write_csv(pred_d, file="papers/pred_d.csv")


## Model based  Predictions ----

pred_mu<- predict_model(samp, pred_mat, new_data, nsamp=5000)

pred_mu<- pred_mu %>% mutate(Nhat= round_any(Nhat,50),se=round_any(se,1),lcl=round_any(lcl, 50),
                             ucl=round_any(ucl, 50))

pred_mu


write_csv(pred_mu, file="papers/pred_mu.csv")


## Plot densities ----
summ<- MCMCsummary(samp, params = "Dbar")
results<- tibble(ID=ID, Density=summ[,"mean"], sd=summ[,"sd"], lcl=summ[,"2.5%"], ucl=summ[,"97.5%"])

allN<-tibble(
  Plot=ID,
  D=as.vector(summ[,"mean"]),
  lcl=as.vector(summ[,"2.5%"]),
  ucl=as.vector(summ[,"97.5%"]),
  ytot=rowMeans(ytot)/Area$Area_3d)

tiff("papers/figure3.tif",width=10,height=5,units="in",res=300, compression="lzw")
#win.graph(10,5)
allN %>% mutate(Plot = fct_reorder(Plot, D, .desc = TRUE)) %>%
  ggplot(aes(x=Plot, y=D))+
  geom_errorbar(aes(ymin=lcl,ymax=ucl), color="blue", width=0.15, size=1) +
  geom_point(colour="blue", size=3, shape=21, fill="white")+
  geom_point(aes(Plot,ytot), colour="red", size=1,shape=16) +
  scale_y_continuous(breaks=seq(0,40,10),lim=c(0, 40))+
  ylab(expression(paste("Density (tahr k",m^-2,")"))) +
  xlab("Plot") +
  theme_bw() +
  theme(text = element_text(size=15),
        axis.text.x = element_blank())
#  theme(text = element_text(size=15),
#        axis.text.x = element_text(angle=90, hjust=1, size=8))
dev.off()

summ<- MCMCsummary(samp, params = "D")

allN<-tibble(
  Plot=rep(ID,3),
  Count_num=rep(c(1,2,3),each=nsites),  
  D=as.vector(summ[,"mean"]),
  lcl=as.vector(summ[,"2.5%"]),
  ucl=as.vector(summ[,"97.5%"]),
  ytot=as.vector(ytot/Area$Area_3d))

tiff("papers/figureS2.tif",width=15,height=15,units="in",res=300, compression="lzw")
#win.graph(15,15)
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
        
dev.off()

# mean naive densities in each MU

summ<- MCMCsummary(out, params = "Dbar")
results<- tibble(ID=ID, Density=summ[,"mean"], sd=summ[,"sd"], lcl=summ[,"2.5%"], ucl=summ[,"97.5%"])
results<- left_join(results, survey.times, by = c("ID" = "PlotID"))
results<- results %>% mutate(ymean = rowMeans(ytot)/Area$Area_3d, ytot=rowSums(ytot))

results %>% group_by(MU) %>% summarise(Dhat=mean(Density),Dnaive=mean(ymean),total=sum(ytot))

tmp<- cbind(survey.times,ytot) %>% arrange(MU)
tmp<- tmp %>% dplyr::select(PlotID,MU,t1,t2,Year,'1','2','3')
write_csv(tmp,"Reports/2019 report/counts.csv")

## Density surface tessellation ----
#
summ<- MCMCsummary(samp, params = "Dbar")
results<- tibble(ID=ID, Density=summ[,"mean"], sd=summ[,"sd"], lcl=summ[,"2.5%"], ucl=summ[,"97.5%"])

xy<- as.data.frame(st_coordinates(st_centroid(MU)))
xy$MU<- c("MU1","MU2","MU3","MU4","MU5","MU6","MU7","EZ1","EZ2")

MU1<- MU %>% mutate(UNIT_NO=lvls_revalue(UNIT_NO,c("MU1","MU2","MU3","MU4","MU5","MU6","MU7","EZ1","EZ2")))

tmp<- plots
tmp<- as.data.frame(st_coordinates(st_centroid(tmp)))
tmp<- tmp %>% mutate(Density=results$Density)
MU1<- as(MU1, "Spatial")

tiff("papers/figure4.tif", width=10, height=10, units="in",res=300, compression="lzw")
#win.graph(10,10)
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
  #  geom_label(aes(X, Y, label=MU), data=xy, size=3,label.padding=unit(0.1, "lines"),hjust="inward",
  #             vjust="inward") +
  labs(x="Longitude",y="Latitude",fill=expression(paste("Density (",km^2,")"))) +
  theme_bw() +
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))
dev.off()

## prediction maps ----
preds<- predict_model(samp, pred_mat, new_data, nsamp=5000, all_cells = TRUE)
preds<- bind_cols(preds$cells, mu_id)

preds %>% group_by(MU) %>% summarise(N = sum(cell_pred))
pred_rast<- tmprast
pred_rast<- mask(pred_rast, MU)
pred_rast<- mask(pred_rast, docsouth)
pred_rast[preds$cell]<- preds$cell_pred/4
pred_rast<- trim(pred_rast)

#tiff("papers/figure5.tif",width=10,height=10,units="in",res=300, compression="lzw")
win.graph(10,10)
MU %>% ggplot() +
  geom_sf(fill="grey80") +
  annotation_scale(location = "br", width_hint = 0.25) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  layer_spatial(pred_rast) +
  scale_fill_viridis(direction=-1, lim=c(0, 30),na.value = "transparent") +
  geom_sf(fill=NA) +
  labs(x="Longitude",y="Latitude",fill=expression(paste("Density (",km^2,")"))) +
  theme_bw() +
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))
dev.off()

## Posterior predictive checks ----

library(gridExtra)

samp<- readRDS("papers/results/mu_nozip_surv.rds")

pp<- MCMCpstr(samp, params = c("yrep"), type="chains")
pp<- t(apply(pp$yrep, c(1,3), sum))

yy<- apply(ytot, 1, sum)

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

m1<- arrangeGrob(p1, p2, p3, p4, ncol=1, nrow=4, top="Non-ZIP model")

samp<- readRDS("papers/results/mu_zip_surv.rds")

pp<- MCMCpstr(samp, params = c("yrep"), type="chains")
pp<- t(apply(pp$yrep, c(1,3), sum))

yy<- apply(ytot, 1, sum)

p1<- ppc_stat(yy, pp, stat=prop_zero, binwidth=0.01) + 
  ggtitle("Proportion zeros") + xlim(0.1, 0.4) + theme_bw() + 
  theme(legend.position = "none", plot.title = element_text(size=10))

p2<- ppc_stat(yy, pp, stat=mean) + ggtitle("Mean") + xlim(30, 40)  + 
  theme_bw() + theme(legend.position = "none", plot.title = element_text(size=10))

p3<- ppc_stat(yy, pp, stat=sd) + ggtitle("Standard deviation") + xlim(40, 60) + 
  theme_bw() + theme(legend.position = "none",plot.title = element_text(size=10))

p4<- ppc_stat(yy, pp, stat=max) + xlim(200, 350) + ggtitle("Maximum") + 
  theme_bw() + theme(legend.position = "none", plot.title = element_text(size=10))

m2<- arrangeGrob(p1, p2, p3, p4, ncol=1, nrow=4, top="ZIP model")

tiff("papers/figureS1.tif", width=6, height=10, units="in",res=300, compression="lzw")
#win.graph(6,10)
grid.arrange(m1,m2, ncol=2)
dev.off()

# Some more
#tiff("papers/figure2.tif", width=6, height=6, units="in",res=300, compression="lzw")
win.graph(6,6)
pp_check(sqrt(yy), sqrt(pp), fun="scatter_avg") + theme_bw()
#dev.off()


win.graph(10,10)
ppc_stat_grouped(yy, pp, stat=prop_zero, group=mu)  + theme_bw()

win.graph(10,10)
ppc_stat_grouped(yy, pp, stat=mean, group=mu) + theme_bw()

win.graph(10,10)
ppc_stat_grouped(yy, pp, stat=sd, group=mu) + theme_bw()

win.graph(10,10)
ppc_stat_grouped(yy, pp, stat=max, group=mu) + theme_bw()

win.graph(8,5)
ppc_dens_overlay(yy, pp[1:1000,]) + theme_bw()

# Some more
#tiff("papers/figure2.tif", width=6, height=6, units="in",res=300, compression="lzw")
win.graph(8,8)
ppc_ecdf_overlay_grouped(yy, pp[1:1000,], group=mu, discrete=TRUE) + theme_bw()

bayes_R2_res(yy, pp)

## Loo-----------------
library(loo)

tmp<-mcmcOutput(samp)
p<- tmp$p
N<- tmp$N
niter<- dim(N)[1]


ll<- matrix(N, nrow=niter, ncol=nsites)
for(i in 1:niter) {
  for(j in 1:nsites)
    ll[i,j]<- sum(dbinom(ytot[j,], N[i,j,], p[i], log=TRUE))
}


loo(ll)
waic(ll)

## Scratch -----------------------
N=100
S = N * 0.49^4
R = N * 0.25 * 4
R + S

summ<- mcmcOutput(samp)

jj<- apply(summ$N, c(2,3), mean)


xx<- MCMCsummary(samp, params="beta")
dm<- model.matrix(~UNIT_NO, data=MU)

lam<- exp(dm %*% xx$mean)

pred_mu<- predict_model(samp, pred_mat, new_data, nsamp=5000, all_cells = TRUE)

jj<- pred_mu$cells
jj<- jj %>% mutate(MU = mu_id$MU)
jj %>% group_by(MU) %>% summarise(N=sum(cell_pred),se=sqrt(sum(sd^2)), cv= se/N)
