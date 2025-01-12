#Mean excess plots
mean.excess.plot.grid<-function(df, gpd.fits=NULL, gpd.fits2=NULL, lower.excesses=FALSE){
  par(mfrow=c(ceiling(ncol(df)/3), 3))
  tickers<-colnames(df)
  sapply(1:ncol(df), function(i){
    #Take observations to the left of the distribution
    if(lower.excesses){
      data<-unlist(df[, i][df[,i]<0]*(-1))
      threshold<-gpd.fits[[i]]$Threshold.lower*(-1)
      if(!is.null(gpd.fits2)){
        threshold<-c(threshold, gpd.fits2[[i]]$Threshold.lower*(-1))
      }
    }
    else{
      data<-unlist(df[, i][df[,i]>0])
      threshold<-gpd.fits[[i]]$Threshold.upper
      if(!is.null(gpd.fits2)){
        threshold<-c(threshold, gpd.fits2[[i]]$Threshold.upper)
      }
    }
    mean_excess_plot(x = data, main = tickers[i], ylab=expression(e[n](u)), xlab="u")
    abline(v=threshold, lty=1:2)
  })
  par(mfrow=c(1,1))
}
mean.excess.plot.grid(df = scenario.set, gpd.fits = fit_double_tail, lower.excesses = F)
mean.excess.plot.grid(df = scenario.set, gpd.fits = fit_double_tail_2, lower.excesses = F)
mean.excess.plot.grid(df = scenario.set, gpd.fits = fit_double_tail, gpd.fits2 = fit_double_tail_2, lower.excesses = F)

#Shape parameter plot
GPD.shape.plot.grid<-function(df, gpd.fits, gpd.fits2=NULL, lower.excesses=FALSE){
  par(mfrow=c(ceiling(ncol(df)/3), 3))
  tickers<-colnames(df)
  for(i in 1:ncol(df)){
      if(lower.excesses){
        data<-unlist(df[, i]*(-1))
        threshold<-gpd.fits[[i]]$Threshold.lower*(-1)
        if(!is.null(gpd.fits2)){
          threshold<-c(threshold, gpd.fits2[[i]]$Threshold.lower*(-1))
        }
      }
      else{
        data<-unlist(df[, i][df[,i]>0])
        threshold<-gpd.fits[[i]]$Threshold.upper
        if(!is.null(gpd.fits2)){
          threshold<-c(threshold, gpd.fits2[[i]]$Threshold.upper)
        }
      }
    tryCatch({
      GPD_shape_plot(x = data, estimate.cov = T, ylab = expression(xi), 
                     xlab = expression(paste("u and ", n[u])), xlab2 = tickers[i])
    }, error=function(e){
      GPD_shape_plot(x = data, estimate.cov = F, ylab = expression(xi), 
                     xlab = expression(paste("u and ", n[u])), xlab2 = tickers[i])
    })
    abline(v=threshold, lty=1:2)
  }
}

GPD.shape.plot.grid(df = scenario.set, gpd.fits = fit_double_tail, gpd.fits2 = fit_double_tail_2)
GPD.shape.plot.grid(df = scenario.set, gpd.fits = fit_double_tail, gpd.fits2 = fit_double_tail_2, lower.excesses = T)

#Table of fitted parameters
fitted.param.tab<-function(gpd.fits, tickers){
  double.tail<-ifelse(length(gpd.fits)>4, TRUE, FALSE)
  if(double.tail){
    shape.lower<-as.numeric(lapply(gpd.fits, function(i) round(i$Shape.lower, 3)))
    se.shape.lower<-as.numeric(lapply(gpd.fits, function(i) round(i$se.shape.lower, 3)))
    scale.lower<-as.numeric(lapply(gpd.fits, function(i) round(i$Scale.lower, 3)))
    se.scale.lower<-as.numeric(lapply(gpd.fits, function(i) round(i$se.scale.lower, 3)))
    threshold.lower<-as.numeric(lapply(gpd.fits, function(i) round(i$Threshold.lower, 3)))
    shape.upper<-as.numeric(lapply(gpd.fits, function(i) round(i$Shape.upper, 3)))
    se.shape.upper<-as.numeric(lapply(gpd.fits, function(i) round(i$se.shape.upper, 3)))
    scale.upper<-as.numeric(lapply(gpd.fits, function(i) round(i$Scale.upper, 3)))
    se.scale.upper<-as.numeric(lapply(gpd.fits, function(i) round(i$se.scale.upper, 3)))
    threshold.upper<-as.numeric(lapply(gpd.fits, function(i) round(i$Threshold.upper, 3)))
    tab<-cbind(shape.lower, se.shape.lower, scale.lower, se.scale.lower,
               threshold.lower, shape.upper, se.shape.upper, scale.upper, 
               se.scale.upper, threshold.upper)
  }
  else{
    shape<-as.numeric(lapply(gpd.fits, function(i) round(i$Shape, 3)))
    se.shape<-as.numeric(lapply(gpd.fits, function(i) round(i$se.shape, 3)))
    scale<-as.numeric(lapply(gpd.fits, function(i) round(i$Scale, 3)))
    se.scale<-as.numeric(lapply(gpd.fits, function(i) round(i$se.scale, 3)))
    threshold<-as.numeric(lapply(gpd.fits, function(i) round(i$Threshold, 3)))
    tab<-cbind(shape, se.shape, scale, se.scale, threshold)
  }
  rownames(tab)<-tickers
  as.data.frame(tab)
}
param_tab_150<-fitted.param.tab(gpd.fits = fit_double_tail, tickers = colnames(scenario.set))
print(xtable(fitte))