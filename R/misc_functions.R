
#--------------------------------------------------
logCI<- function(N, varN, alpha=0.95) {
  # log transformed (1-alpha) confidence intervals
  cv<- sqrt(varN)/N
  cv2<- cv^2
  asymp.ci <- exp(qnorm((1-alpha)/2) * sqrt(log(1+cv2)))
  list(lcl=N*asymp.ci, ucl=N/asymp.ci)
}

calc_nplots<- function(plots, MU, PCL, res=2000) {
  # calc number of plots that could be sampled in each MU
  # for given resolution
  # including only cells in PCL
  
  tmprast<- raster::raster(MU, ext=raster::extent(MU), resolution=res)
  murast<- raster::rasterize(as_Spatial(MU), tmprast, field="UNIT_NO")
  murast<- raster::mask(murast, as_Spatial(PCL))
  plots.in.units<- st_nn(st_centroid(plots), MU, sparse=FALSE)
  sampled.plots<- apply(plots.in.units,2,sum)
  nplots<- as.integer(table(raster::values(murast)))
  nunits<- tibble(MU=levels(MU$UNIT_NO),sampled=sampled.plots, nplots=nplots)
  nunits
}


#--------------------

calc_strata_nhat<- function(est, ncells, A, CI.level=0.95, type = c("abundance","density")){
  # design-based estimates of abundance for stratum (MU)   
  type<- match.arg(type)
  m<- nrow(est)
  if(m > 1) {
    M<- ncells$nplots
    ip<- m/M # equal probability sampling
    TA<- M*A # Total area of MU
    nhat<- est$Nhat
    varn<- est$varN
    ti <- m * nhat/ip
    N<- sum(nhat/ip)  # Alternatively- mean(ti)
    S2N<- sum((ti - N)^2)/(m-1)
    varN<- (1-m/M) * S2N/m + sum(varn/ip)
    SD = sqrt(varN)
    CV = SD/N
    D<- N/TA
    varD<- varN/TA^2
    climitsN<- logCI(N, varN, alpha=CI.level)
    climitsD<- logCI(D, varD, alpha=CI.level)
    thar.mu<- tibble(Nhat=round(N), SD=round(SD), CV=round(CV, 2), lcl=round(climitsN$lcl), 
                     ucl=round(climitsN$ucl),m=m, M=M)
    thar.D<- tibble(Dhat=round(D,2), SD=round(sqrt(varD),2), CV=round(sqrt(varD)/D, 2), 
                   lcl=round(climitsD$lcl,2), ucl=round(climitsD$ucl,2), m=m, M=M)
  }
  
  if(identical(type, "abundance")) ans<- thar.mu
  else if(identical(type, "density")) ans<- thar.D
  else stop("unknown type")
  return(ans)
}

calc_total_nhat<- function(x, CI.level=0.95) {
  Ntot<- round(sum(x$Nhat, na.rm=T))
  varN<- round(sum(x$SD^2, na.rm=T))
  SD<- sqrt(varN)
  CL<- logCI(Ntot, varN, alpha=CI.level)
  Ntot<- tibble(MU="Total", Nhat=Ntot, SD=SD, CV=round(SD/Ntot,2),
                lcl=round(CL$lcl), ucl=round(CL$ucl),m=sum(x$m), M=sum(x$M))
  Ntot
}

round_any<- function(x, r){
  round(x/r)*r
}

#---------------------

predict_model<- function(samples, X, pred_df, nsamp=1000, all_cells = FALSE) {
  post <- MCMCpstr(samples, params = "all", type="chains")
  nn <- dim(post$beta)[2]
  ii <- sample(1:nn, size=nsamp, replace=FALSE)
  beta <- post$beta[,ii]
  if(!is.null(post$theta)) {
    theta <- post$theta[,ii]
    ZIP<- TRUE
  }
  else ZIP<- FALSE
  npred <- nrow(X)
  out <- matrix(NA, npred, nsamp)
  pb <- progress_bar$new(total = nsamp)
  pb$message("Generating predictions...")
  for(i in 1:nsamp) {
    if(ZIP)
      P <- rbinom(npred, 1, theta[i])
    else
      P<- 1
    lambda <- exp(X %*% beta[,i])
    out[,i] <- rpois(npred, lambda * P)
    pb$tick()
  }
  # do management units 
  strata<- levels(pred_df$MU)
  strata.out <- matrix(0, nrow = length(strata), ncol=6)
  pb2 <- progress_bar$new(total = length(strata))
  pb2$message("Calculating strata totals...")
  for(i in 1:length(strata)) {
    inds <- which(pred_df$MU == strata[i])
    submat <-  out[inds,]
    Nsum <- apply(submat,2,sum)
    Ntot <- mean(Nsum)
    se <- sd(Nsum)
    cv <- ifelse(!is.finite(se/Ntot), 0, se/Ntot)
    lcl <- quantile(Nsum, 0.025)
    ucl <- quantile(Nsum, 0.975)
    strata.out[i,] <- c(round(Ntot), round(se), round(cv,2),round(lcl),round(ucl),length(inds))
    pb2$tick()
  }
  strata.out <- data.frame(strata.out)
  names(strata.out) <- c("Nhat","se","cv","lcl","ucl","n")
  strata.out <- cbind(strata, strata.out)
  # now do overall
  Nsum <- apply(out,2,sum)
  Ntot <- median(Nsum)
  se <- mad(Nsum)
  cv <- ifelse(!is.finite(se/Ntot), 0, se/Ntot)
  lcl <- quantile(Nsum, 0.025)
  ucl <- quantile(Nsum, 0.975)
  df <- data.frame(strata="Total", Nhat = round(Ntot), se=round(se), cv=round(cv,2),
                   lcl=round(lcl),ucl=round(ucl),n=nrow(X))
  strata.out<- bind_rows(strata.out, df)
  row.names(strata.out) <- NULL
  # All cells
  if(all_cells) {
    cell_preds<- apply(out, 1, mean)
    cell_sd<- apply(out, 1, sd)
    cell_preds<- tibble(cell_pred=cell_preds, sd=cell_sd)
    strata.out<- list(strata=strata.out, cells=cell_preds)
  }
  strata.out
}

