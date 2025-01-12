#Table of summary statistics for each asset
library(PerformanceAnalytics)
library(tsoutliers)
scenario.set.p<-loss.dist(price_df = prices_ss, names = colnames(prices_ss)[-1], 
                          ret_type = "log", neg_rets = F, no_zeroes = F)
scenario.set.n<-loss.dist(price_df = prices_ss, names = colnames(prices_ss)[-1], 
                          ret_type = "simple", neg_rets = F, no_zeroes = F)
exp(colSums(scenario.set.p))
tail(prices_ss[,-1],1)/head(prices_ss[,-1],1)-1
cagr<-function(scenario.set, neg_rets=FALSE){
  if(neg_rets){
    scenario.set<-scenario.set*(-1)
  }
  time<-nrow(scenario.set)/260
  cagr<-(exp(colSums(scenario.set))^(1/time)-1)*100
  cagr
}
cagr(scenario.set.p)
exp(colSums(scenario.set.p))^(1/11.27)-1

ann_sd<-function(scenario.set, log_rets=TRUE){
  if(log_rets){
    scenario.set<-exp(scenario.set)-1
  }
  std_dev<-apply(scenario.set, MARGIN = 2, sd)
  ann_sd<-std_dev*sqrt(260)
  ann_sd*100
}
ann_sd(scenario.set)
ann_sd(scenario.set.n)

skew<-function(scenario.set, neg_rets=FALSE){
  if(neg_rets){
    scenario.set<-scenario.set*(-1)
  }
  skew<-apply(scenario.set, MARGIN = 2, FUN = skewness)
  skew
}
skew(scenario.set, neg_rets = T)

kurt<-function(scenario.set){
  kurt<-apply(scenario.set, MARGIN = 2, FUN = kurtosis)
  kurt
}
kurt(scenario.set)

VaR_ss<-function(scenario.set, alpha, neg_rets=FALSE, log_rets=TRUE){
  if(neg_rets==FALSE){
    scenario.set<-scenario.set*(-1)
  }
  var<-numeric(ncol(scenario.set))
  var<-sapply(1:ncol(scenario.set), function(x) 
    quantile(scenario.set[,x], alpha))
  if(log_rets){
    unname((exp(var)-1)*100)
  }
  else{
    unname(var*100 )
  }
}
VaR_ss(scenario.set = scenario.set, alpha = 0.99, neg_rets = T, log_rets = T)

ES_ss<-function(scenario.set, alpha, neg_rets=FALSE, log_rets=TRUE){
  if(neg_rets==FALSE){
    scenario.set<-scenario.set*(-1)
  }
  ES<-numeric(ncol(scenario.set))
  for(i in 1:ncol(scenario.set)){
    ES[i]<-ES.discrete(loss = scenario.set[,i], alpha = alpha, tail = "right")
  }
  if(log_rets){
    unname((exp(ES)-1)*100)
  }
  else{
    unname(ES*100)
  }
}
ES_ss(scenario.set = scenario.set, alpha = 0.99, neg_rets = T)

#Summary table
summary_stats<-function(scenario.set, neg_rets, log_rets, alpha){
  sum_tab<-as.data.frame(matrix(data = NA, nrow = ncol(scenario.set), ncol = 7))
  colnames(sum_tab)<-c("Asset", "CAGR", "Ann. SD", "Skewness",
                       "Kurtosis", paste(alpha, "VaR"), paste(alpha, "ES"))
  sum_tab[,1]<-colnames(scenario.set)
  sum_tab[,2]<-cagr(scenario.set, neg_rets)
  sum_tab[,3]<-ann_sd(scenario.set, log_rets)
  sum_tab[,4]<-skew(scenario.set, neg_rets)
  sum_tab[,5]<-kurt(scenario.set)
  sum_tab[,6]<-VaR_ss(scenario.set, alpha, neg_rets, log_rets)
  sum_tab[,7]<-ES_ss(scenario.set, alpha, neg_rets, log_rets)
  sum_tab
}
summary_stats(scenario.set = scenario.set, neg_rets = T, log_rets = T, alpha = 0.95)

#Table with optimal portfolios for different set-ups
opt_port_tab<-function(scenario.set.list, method=c("ES", "sd"), conf_level=NULL){
  optimals_tab<-as.data.frame(matrix(NA, nrow = length(scenario.set.list), 
                                     ncol = ncol(scenario.set.list[[1]])+2))
  for(i in 1:nrow(optimals_tab)){
    optimals_tab[i, 3:ncol(optimals_tab)]<-stoch.optim(scenario.set = scenario.set.list[[i]], 
                                                       method = method, conf_level = conf_level, 
                                                       mean = 0)
  }
  colnames(optimals_tab)<-c("Scenario set", "Method", colnames(scenario.set.list[[1]]))
  optimals_tab[,1]<-names(scenario.set.list)
  optimals_tab[,2]<-rep(method, nrow(optimals_tab))
  optimals_tab
}

ss_list<-list("scenario.set"=scenario.set*(-1), "ssboot_double_tail"=ssboot_double_tail, 
              "ssvine_double_tail"=ssvine_double_tail)
debug(opt_port_tab)
es_tab<-opt_port_tab(scenario.set.list = ss_list, method = "ES", conf_level = 0.05)
sd_tab<-opt_port_tab(scenario.set.list = ss_list, method = "sd")
opt_ports<-rbind(es_tab, sd_tab)
opt_ports<-opt_ports[order(opt_ports[,1]),]
opt_ports

#Risk report for optimal portfolios
risk_report<-function(opt_weights, scenario.set){
  loss<-scenario.set %*% as.numeric(opt_weights)
  time<-length(loss)/260
  cagr<-(exp(sum(loss))^(1/time)-1)*100
  ann_sd<-(sd(exp(loss)-1)*sqrt(260))*100
  var95<-unname(quantile(loss*(-1), 0.95)*100)
  es95<-ES.discrete(loss = loss*(-1), alpha = 0.95, tail = "right")*100
  es99<-ES.discrete(loss = loss*(-1), alpha = 0.99, tail = "right")*100
  c("CAGR"=cagr, "Ann. SD"=ann_sd, "VaR 95%"=var95, "ES 95%"=es95, "ES 99%"=es99)
}

risk_report(opt_ports[1, 3:15], scenario.set = scenario.set*(-1))
risk_report(opt_ports[2, 3:15], scenario.set = scenario.set*(-1))
risk_report(opt_ports[3, 3:15], scenario.set = scenario.set*(-1))
risk_report(opt_ports[4, 3:15], scenario.set = scenario.set*(-1))
risk_report(opt_ports[5, 3:15], scenario.set = scenario.set*(-1))
risk_report(opt_ports[6, 3:15], scenario.set = scenario.set*(-1))

risk_report_tab<-function(opt_weights_tab, scenario.set.list){
  risk_rep_tab<-as.data.frame(matrix(NA, nrow = nrow(opt_weights_tab), ncol = 7))
  risk_rep_tab[1, 3:7]<-risk_report(opt_weights_tab[1, 3:15], scenario.set.list[[1]]*(-1))
  risk_rep_tab[2, 3:7]<-risk_report(opt_weights_tab[2, 3:15], scenario.set.list[[1]]*(-1))
  risk_rep_tab[3, 3:7]<-risk_report(opt_weights_tab[3, 3:15], scenario.set.list[[2]]*(-1))
  risk_rep_tab[4, 3:7]<-risk_report(opt_weights_tab[4, 3:15], scenario.set.list[[2]]*(-1))
  risk_rep_tab[5, 3:7]<-risk_report(opt_weights_tab[5, 3:15], scenario.set.list[[3]]*(-1))
  risk_rep_tab[6, 3:7]<-risk_report(opt_weights_tab[6, 3:15], scenario.set.list[[3]]*(-1))
  risk_rep_tab[,1]<-opt_weights_tab[,1]
  risk_rep_tab[,2]<-opt_weights_tab[,2]
  colnames(risk_rep_tab)<-c("Scenario Set", "Method", "CAGR", "Ann. SD", 
                            "VaR 95%", "ES 95%", "ES 99%")
  risk_rep_tab
}
risk_report_tab(opt_weights_tab = opt_ports, scenario.set.list = ss_list)
