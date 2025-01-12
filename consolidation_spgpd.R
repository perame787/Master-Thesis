#Testing the new semiparametric distribution
tickers<-c("IVV", "EZU", "JPXN", "EPP", "CEC.PA", "EEM", "VCIT", "EMB", "TRET.L", "GSG",
           "DE10Y", "USFI10Y", "NEIXCTAT")
prices_ss<-fread("C:/Users/Pedro/Documents/Studium/WU Dokummente/SoSe 23/Thesis/Scripts/Data/prices_ss.csv")
scenario.set<-loss.dist(price_df = prices_ss, names = names(prices_ss)[-1], 
                        ret_type = "log", neg_rets = T, no_zeroes = F)
scenario_set<-loss.dist(price_df = prices_ss, names = names(prices_ss)[-1], 
                        ret_type = "log", neg_rets = F, no_zeroes = F)

spgpd.fit.list<-function(scenario.set, threshold.q, method="MLE", lower.tail=FALSE){
  tickers<-colnames(scenario.set)
  fits<-vector(mode="list", length(ncol(scenario.set)))
  fits<-lapply(1:length(tickers), function(i){
    loss_dist<-scenario.set[, tickers[i]]
    gpd.fit(data = loss_dist, threshold.q, method, lower.tail)
  })
  names(fits)<-paste0(colnames(scenario.set), ".fit")
  fits
}
undebug(spgpd.fit.list)
fits.spgpd<-spgpd.fit.list(scenario.set = scenario.set, threshold.q = 0.95, lower.tail = F)  

copula.emp<-function(scenario.set, fits, method=c("SPGPD", "spd")){
  copula<-matrix(nrow = nrow(scenario.set), ncol = ncol(scenario.set))
  if(method=="SPGPD"){
    for(i in 1:ncol(scenario.set)){
      fit_i<-fits[[i]]
      ss_i<-scenario.set[,i]
      copula[,i]<-sapply(ss_i, function(x) pSPGPD(x, data = ss_i, fit = fit_i))
    }
  }
  else if(method=="spd"){
    copula<-as.data.frame(sapply(1:ncol(scenario.set), function(i){
      copula[,i]<-pspd(scenario.set[, i], fit = fits[[i]])
    })) 
  }
  colnames(copula)<-paste0(colnames(scenario.set), ".q")
  copula
}
debug(copula.emp)
emp.copula<-copula.emp(scenario.set = scenario.set, fits = fits.spgpd,method = "SPGPD")

inverse.copula<-function(copula_df, scenario.set=NULL, fits, method=c("SPGPD", "spd"), bootstrap=FALSE){
  inverse_df<-matrix(nrow = nrow(copula_df), ncol = ncol(copula_df))
  #Create bootstraped copula_df
  if(bootstrap==TRUE){
    copula_df<-copula_df[sample(x = 1:nrow(copula_df), size = nrow(copula_df), replace = T),] 
  }
  #Apply inverse transform for each column of the bootstraped copula
  if(method=="SPGPD"){
    for(i in 1:ncol(copula_df)){
      fit_i<-fits[[i]]
      q_i<-copula_df[,i]
      ss_i<-scenario.set[,i]
      inverse_df[,i]<-sapply(q_i, function(x) qSPGPD(prob = x, data = ss_i, fit = fit_i))
    }
  }
  else if(method=="spd"){
    inverse_df<-as.data.table(sapply(1:ncol(copula_df), function(i){
      inverse_df[,i]<-qspd(copula_df[, i], fit = fits[[i]])
    }))
  }
  names(inverse_df)<-sub("q", "boot", names(copula_df))
  inverse_df
}
undebug(inverse.copula)
set.seed(123)
ss_2<-inverse.copula(copula_df = emp.copula, scenario.set = scenario.set,
                         fits = fits.spgpd, method = "SPGPD", bootstrap = T)

#Attempt to generalize
gpd.fit.ss<-function(scenario.set, upper.threshold=NULL, lower.threshold=NULL, 
                     lower.tail=FALSE, method=c("SPGPD", "spd", "GNP"), kernelfit=NULL){
  tickers<-colnames(scenario.set)
  fits<-vector(mode="list", length(ncol(scenario.set)))
  if(method=="SPGPD"){
    if(lower.tail==TRUE){
      upper.threshold<-lower.threshold
    }
    fits<-lapply(1:length(tickers), function(i){
      loss.dist<-scenario.set[, tickers[i]]
      gpd.fit(data = loss.dist, threshold.q = upper.threshold, method = "MLE" ,lower.tail = lower.tail)
    })
  }
  fits<-lapply(1:length(tickers), function(i){
    loss.dist<-scenario.set[, tickers[i]]
    spdfit(data = loss.dist, upper = upper.threshold, lower = lower.threshold, tailfit = "GPD", type="mle", kernelfit)
  })
  names(fits)<-paste0(colnames(scenario_set), ".fit")
  fits
}
debug(gpd.fit.ss)
gpd.fit.ss(scenario.set = scenario.set, upper.threshold = NULL, lower.threshold = 0.1, lower.tail = TRUE, method = "SPGPD")

#Comparison one-tailed SPGPD vs. two-tailed SPGPD  vs. spd
#Fits for one-tail SPGPD
spgpd.fit.one.tail<-spgpd.fit.list(scenario.set = scenario.set, upper = 0.9, method = "MLE")
spgpd.fit.double.tail<-spgpd.fit.list(scenario.set = scenario.set, lower = 0.1, 
                                      upper = 0.9, method = "MLE", double.tail = TRUE)
spd.fit<-semi.gpd.ss(scenario.set = scenario.set, upper = 0.9, lower = 0.1, 
                     tailfit = "GPD", type = "mle", kernelfit = "epanech")

emp.copula.one.tail<-copula.emp(scenario.set = scenario.set, fits = spgpd.fit.one.tail, 
                                method = "SPGPD")
emp.copula.double.tail<-copula.emp(scenario.set = scenario.set, fits = spgpd.fit.double.tail, 
                                method = "SPGPD")
emp.copula.spd<-copula.emp(scenario.set = scenario.set, fits = spd.fit, method = "spd")

#Revision of similar values in copular
#Create bootstraped copulas
set.seed(123)
boot_1t<-emp.copula.one.tail[sample(x = 1:nrow(emp.copula.one.tail), size = nrow(emp.copula.one.tail), replace = T),] 
ss_one_tail<-inverse.copula(copula_df = boot_1t, scenario.set = scenario.set, fits = spgpd.fit.one.tail, method = "SPGPD", bootstrap = F)
set.seed(123)
boot_2t<-emp.copula.double.tail[sample(x = 1:nrow(emp.copula.double.tail), size = nrow(emp.copula.double.tail), replace = T),] 
ss_double_tail<-inverse.copula(copula_df = boot_2t, scenario.set = scenario.set, fits = spgpd.fit.double.tail, method = "SPGPD", bootstrap = F)

#Create scenario set using each method with a fixed seed1
#Problem: bootstrap draws a lot from the center, leading to very similar sets
set.seed(123)
ss_one_tail<-inverse.copula(copula_df = emp.copula.one.tail, scenario.set = scenario.set, 
                                fits = spgpd.fit.one.tail, method = "SPGPD", bootstrap = TRUE)
set.seed(123)
ss_double_tail<-inverse.copula(copula_df = emp.copula.double.tail, scenario.set = scenario.set, 
                                fits = spgpd.fit.double.tail, method = "SPGPD", bootstrap = TRUE)
set.seed(123)
ss_spd<-inverse.copula(copula_df = emp.copula.spd, scenario.set = scenario_set, 
                                fits = spd.fit, method = "spd", bootstrap = TRUE)

stoch.optim(scenario.set = ss_one_tail*(-1), method = "ES", conf_level = 0.95)
stoch.optim(scenario.set = ss_double_tail*(-1), method = "ES", conf_level = 0.95)
stoch.optim(scenario.set = ss_spd, method = "ES", conf_level = 0.95)

spgpd1_eff<-eff.frontier(scenario.set = ss_one_tail*(-1), mu_sigma_space = "ES", conf_level = 0.01)
spgpd2_eff<-eff.frontier(scenario.set = ss_double_tail*(-1), mu_sigma_space = "ES", conf_level = 0.01)
spgpd2_eff_sd<-eff.frontier(scenario.set = ss_double_tail*(-1), mu_sigma_space = "sd", conf_level = 0.01)
spd_eff<-eff.frontier(scenario.set = ss_spd, mu_sigma_space = "ES", conf_level = 0.01)
emp_eff<-eff.frontier(scenario.set = scenario_set, mu_sigma_space = "ES", conf_level = 0.1)
emp_eff_sd<-eff.frontier(scenario.set = scenario_set, mu_sigma_space = "sd", conf_level = 0.01)

plot(spgpd1_eff$frontier_es99, spgpd1_eff$frontier_mean, type="l", col="blue2",
     xlim=c(0, max(spgpd1_eff$frontier_es99)+0.01), 
     ylim=c(0, max(spgpd1_eff$frontier_mean)+0.0001), xlab="ES", ylab=expression(mu),
     main="Semiparametric scenario sets")
lines(spgpd2_eff$frontier_es99, spgpd2_eff$frontier_mean, 
      type="l", col="red2",lty=2)
lines(spd_eff$frontier_es99, spd_eff$frontier_mean, 
      type="l", col="green3",lty=2)
lines(emp_eff$frontier_es99, emp_eff$frontier_mean, 
      type="l", col="orange2",lty=2)

plot(spgpd2_eff$frontier_es90, spgpd2_eff$frontier_mean, type="l", col="blue2",
       xlim=c(0, max(spgpd2_eff$frontier_es90)+0.01), 
       ylim=c(0, max(spgpd2_eff$frontier_mean)+0.0001), xlab="ES", ylab=expression(mu),
       main="ES vs. Markowitz")
lines(spgpd2_eff_sd$frontier_es90, spgpd2_eff_sd$frontier_mean, 
        type="l", col="red2",lty=2)

#Eff. frontier with scenario from copula
set.seed(123)
ss_vc<-inverse.copula(vinecop = cop_4, scenario.set = scenario.set, 
                      fits = spgpd.fit.double.tail, method = "SPGPD")

spgpd2_eff_vc<-eff.frontier(scenario.set = ss_vc*(-1), mu_sigma_space = "ES", 
                            conf_level = 0.01)
spgpd2_eff_sd_vc<-eff.frontier(scenario.set = ss_vc*(-1), mu_sigma_space = "sd", 
                               conf_level = 0.01)

plot(spgpd2_eff_vc$frontier_es95, spgpd2_eff_vc$frontier_mean, type="l", col="blue2",
     xlim=c(0, max(spgpd2_eff_vc$frontier_es95)+0.01), 
     ylim=c(0, max(spgpd2_eff$frontier_mean)+0.0001), xlab="ES 95%", ylab=expression(mu),
     main="Copula vs. Bootstrap vs. Empirical scenario sets")
lines(spgpd2_eff$frontier_es95, spgpd2_eff$frontier_mean, 
      type="l", col="red2",lty=2)
lines(emp_eff$frontier_es95, emp_eff$frontier_mean, 
      type="l", col="green3",lty=2)
legend("bottomright", legend = c("Copula", "Bootstrap", "Empirical"), 
       lty = c(1,2,2), col = c("blue2", "red2", "green3"))
