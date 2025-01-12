#Generate list of scenario sets
scen.generator.list<-function(copula_df=NULL, vinecop=NULL, scenario.set=NULL, fits, 
                              method="SPGPD", qrng=TRUE, n=1000, t=nrow(scenario.set)){
  names<-colnames(scenario.set)
  ncol_ss<-length(names)
  #Generate list with empty matrices
  scenario_set_list<-lapply(1:n, matrix, data=NA, nrow=t, ncol=ncol_ss)
  if(is.null(vinecop)){
    scenario_set_list<-lapply(scenario_set_list, function(i){ 
      scen.generator(copula_df = copula_df, scenario.set = scenario.set, fits = fits, 
                     method = method, qrng = qrng, t = t)
    })
  }
  else{
    scenario_set_list<-lapply(scenario_set_list, function(i){ 
      scen.generator(vinecop = vinecop, scenario.set = scenario.set, fits = fits, 
                     method = method, qrng = qrng, t = t)
    })
  }
  scenario_set_list
}

#start<-Sys.time()
#scenario.set.list<-lapply(1:3, matrix, data=NA, nrow=2929, ncol=13)
#scenario.set.list<-lapply(scenario.set.list, function(i){
#  scen.generator(vinecop = vine_dtail, scenario.set = scenario.set, fits = fit_double_tail, t = 2929)
#})
#end<-Sys.time()
#end-start
scenario.set.list<-scen.generator.list(vinecop = vine_dtail, 
                                       scenario.set = scenario.set, fits = fit_double_tail, n = 3)
start<-Sys.time()
set.seed(123)
scenario.set.list2<-scen.generator.list(copula_df = emp_cop_dtail, 
                                        scenario.set = scenario.set, 
                                        fits = fit_double_tail, n = 1000)
end<-Sys.time()
end-start

#Monte Carlo simulation for portfolio optimzation - Eff. Frontier with CI
eff.portfolio.mc<-function(scenario_set_list, method=c("ES", "var"), 
                           conf_level=NULL, min_weight=NULL, max_weight=NULL, 
                           pos_limit=NULL, mean=NULL){
  names<-colnames(scenario.set)
  opt.weight.matrix<-matrix(data = NA, nrow = length(scenario_set_list), ncol = length(names))
  colnames(opt.weight.matrix)<-names
  opt.weight.matrix<-t(sapply(scenario_set_list, function(i){
    stoch.optim(scenario.set = i, method = method, conf_level = conf_level, 
                min_weight = min_weight, max_weight = max_weight, 
                pos_limit = pos_limit, mean = mean)
  }))
  weights_mc<-colMeans(opt.weight.matrix)
  std_errors<-apply(opt.weight.matrix, MARGIN = 2, function(i) sd(i)/sqrt(length(i)))
  rbind("Weight"=weights_mc, "Std. error"=std_errors)
}

start<-Sys.time()
s.optim3<-eff.portfolio.mc(scenario_set_list = scenario.set.list2, method = "var", conf_level = 0.95)
end<-Sys.time()
end-start

#Efficient frontier with confidence interval
frontier1<-eff.frontier(scenario.set = scenario.set*(-1), mu_sigma_space = "var", conf_level = 0.95)
frontier2<-eff.frontier(scenario.set = scenario.set*(-1), mu_sigma_space = "ES", conf_level = 0.95)
frontier3<-eff.frontier.b(scenario.set = scenario.set*(-1), mu_sigma_space = "ES", conf_level = 0.95)

#Function for MC estimate plus columns with lower and upper values for the CI
mc_ci<-function(param_matrix, conf_level_mc, param_name){
  mc_matrix<-matrix(data = NA, nrow = nrow(param_matrix), ncol = 3)
  mc_matrix[,2]<-rowMeans(param_matrix)
  sd<-apply(param_matrix, MARGIN = 1, FUN = sd)
  bounds<-(qnorm((1-(conf_level_mc/2)))*sd)/sqrt(nrow(param_matrix))
  mc_matrix[,1]<-mc_matrix[,2]-bounds
  mc_matrix[,3]<-mc_matrix[,2]+bounds
  colnames(mc_matrix)<-c(paste0(param_name, "_lb"), paste0(param_name, "_mc"), 
                         paste0(param_name, "_ub"))
  mc_matrix
}

eff.frontier.mc<-function(scenario_set_list, mu_sigma_space, conf_level, conf_level_mc=0.05){
  n<-length(scenario_set_list)
  frontier_size<-ncol(scenario_set_list[[1]])+5
  mu_matrix<-matrix(data = NA, nrow = frontier_size, ncol = n)
  sd_matrix<-matrix(data = NA, nrow = frontier_size, ncol = n)
  es95_matrix<-matrix(data = NA, nrow = frontier_size, ncol = n)
  es99_matrix<-matrix(data = NA, nrow = frontier_size, ncol = n)
  eff_frontier_list<-lapply(1:n, matrix, data=NA, nrow=frontier_size, ncol=4)
  eff_frontier_list<-lapply(scenario_set_list, function(i){
    eff.frontier.b(scenario.set = i*(-1), mu_sigma_space = mu_sigma_space, conf_level = conf_level)
  })
  mu_matrix<-sapply(eff_frontier_list, function(i) i[,1])
  sd_matrix<-sapply(eff_frontier_list, function(i) i[,2])
  es95_matrix<-sapply(eff_frontier_list, function(i) i[,3])
  es99_matrix<-sapply(eff_frontier_list, function(i) i[,4])
  mu_mc<-mc_ci(param_matrix = mu_matrix, conf_level_mc = conf_level_mc, param_name = "mu")
  sd_mc<-mc_ci(param_matrix = mu_matrix, conf_level_mc = conf_level_mc, param_name = "sd")
  es95_mc<-mc_ci(param_matrix = es95_matrix, conf_level_mc = conf_level_mc, param_name = "ES95")
  es99_mc<-mc_ci(param_matrix = es99_matrix, conf_level_mc = conf_level_mc, param_name = "ES99")
  as.data.frame(do.call(cbind,(list(mu_mc, sd_mc, es95_mc, es99_mc))))
}

#Test efficiency frontier and plotting with n=10 sets
scenario.set.list<-scen.generator.list(vinecop = vine_dtail, fits = fit_double_tail, n = 5, scenario.set = scenario.set)
start<-Sys.time()
test1_mc<-eff.frontier.mc(scenario_set_list = scenario.set.list, mu_sigma_space = "ES", conf_level = 0.05)
end<-Sys.time()
end-start

#Parallelize efficient frontier
library(parallel)
eff.frontier.mc.par<-function(scenario_set_list, mu_sigma_space, conf_level, conf_level_mc=0.05){
  n<-length(scenario_set_list)
  frontier_size<-ncol(scenario_set_list[[1]])+5
  mu_matrix<-matrix(data = NA, nrow = frontier_size, ncol = n)
  sd_matrix<-matrix(data = NA, nrow = frontier_size, ncol = n)
  es95_matrix<-matrix(data = NA, nrow = frontier_size, ncol = n)
  es99_matrix<-matrix(data = NA, nrow = frontier_size, ncol = n)
  eff_frontier_list<-lapply(1:n, matrix, data=NA, nrow=frontier_size, ncol=4)
  #Parallelize
  cl<-makeCluster(detectCores())
  clusterEvalQ(cl, source("C:/Users/Pedro/Documents/Studium/WU Dokummente/SoSe 23/Thesis/Scripts/package.R"))
  eff_frontier_list<-parLapply(cl = cl, X = scenario_set_list, function(i){
    eff.frontier.b(scenario.set = i*(-1), mu_sigma_space = mu_sigma_space, conf_level = conf_level)
  })
  stopCluster(cl)
  mu_matrix<-sapply(eff_frontier_list, function(i) i[,1])
  sd_matrix<-sapply(eff_frontier_list, function(i) i[,2])
  es95_matrix<-sapply(eff_frontier_list, function(i) i[,3])
  es99_matrix<-sapply(eff_frontier_list, function(i) i[,4])
  mu_mc<-mc_ci(param_matrix = mu_matrix, conf_level_mc = conf_level_mc, param_name = "mu")
  sd_mc<-mc_ci(param_matrix = mu_matrix, conf_level_mc = conf_level_mc, param_name = "sd")
  es95_mc<-mc_ci(param_matrix = es95_matrix, conf_level_mc = conf_level_mc, param_name = "ES95")
  es99_mc<-mc_ci(param_matrix = es99_matrix, conf_level_mc = conf_level_mc, param_name = "ES99")
  as.data.frame(do.call(cbind,(list(mu_mc, sd_mc, es95_mc, es99_mc))))
}

start<-Sys.time()
test1_mc_par<-eff.frontier.mc.par(scenario_set_list = scenario.set.list, mu_sigma_space = "ES", conf_level = 0.05)
end<-Sys.time()
end-start

plot(test1_mc$ES95_mc, test1_mc$mu_mc, type="l")
lines(test1_mc$ES95_lb, test1_mc$mu_lb, type="l", lty=2)
lines(test1_mc$ES95_ub, test1_mc$mu_ub, type="l", lty=2)


#Plotting with CI's
set.seed(1234)
df <- data.frame(x =1:10,
                 F =runif(10,1,2),
                 L =runif(10,0,1),
                 U =runif(10,2,3))


plot(df$x, df$F, ylim = c(0,4), type = "l")
#make polygon where coordinates start with lower limit and 
# then upper limit in reverse order
polygon(c(df$x,rev(df$x)), c(df$L,rev(df$U)),col = "grey75", border = FALSE)
lines(df$x, df$F, lwd = 2)
#add red lines on borders of polygon
lines(df$x, df$U, col="red",lty=2)
lines(df$x, df$L, col="red",lty=2)

plot(test1_mc_par$ES95_mc, test1_mc_par$mu_mc, type="l", 
     xlim=c(min(test1_mc_par$ES95_mc+0.0001), max(test1_mc_par$ES95_mc+0.0001)),
     ylim=c(min(test1_mc_par$mu_mc+0.00001), max(test1_mc_par$mu_mc+0.00001)), xlab="ES95-MC", 
     ylab=expression(paste(mu, "-MC")))
polygon(c(test1_mc_par$ES95_mc, rev(test1_mc_par$ES95_mc)), c(test1_mc_par$mu_lb, rev(test1_mc_par$mu_ub)), 
        col="grey90", border=F)
lines(test1_mc_par$ES95_mc, test1_mc_par$mu_mc, col="grey15")
lines(test1_mc_par$ES95_mc, test1_mc_par$mu_lb, col="red", lty=2)
lines(test1_mc_par$ES95_mc, test1_mc_par$mu_ub, col="red", lty=2)

plot(eff_boot_es95$ES95_mc, eff_boot_es95$mu_mc, type="l", 
     xlim=c(min(eff_boot_es95$ES95_mc+0.0001), max(eff_boot_es95$ES95_mc+0.0001)),
     ylim=c(0, max(eff_boot_es95$mu_ub+0.00003)), xlab="ES95-MC", 
     ylab=expression(paste(mu, "-MC")))
polygon(c(eff_boot_es95$ES95_mc, rev(eff_boot_es95$ES95_mc)), c(eff_boot_es95$mu_lb, rev(eff_boot_es95$mu_ub)), 
        col="grey90", border=F)
lines(eff_boot_es95$ES95_mc, eff_boot_es95$mu_mc, col="red")
lines(eff_boot_es95$ES95_mc, eff_boot_es95$mu_lb, col="red3", lty=2)
lines(eff_boot_es95$ES95_mc, eff_boot_es95$mu_ub, col="red3", lty=2)

#Function to plot eff. frontiers with CI
plot_eff_frontier<-function(eff_frontier_list, hist_ef_included=T, risk_measure, col_vector){
  if(hist_ef_included){
    if(risk_measure=="ES95"){
      min_x<-min(min(eff_frontier_list[[1]]$frontier_es95), 
                 sapply(eff_frontier_list[2:3], function(i) min(i$ES95_mc)))
      max_x<-max(max(eff_frontier_list[[1]]$frontier_es95), 
                 sapply(eff_frontier_list[2:3], function(i) max(i$ES95_mc)))
      max_y<-max(sapply(eff_frontier_list[2:3], function(i) max(i$mu_ub)))
      plot(eff_frontier_list[[1]]$frontier_es95, eff_frontier_list[[1]]$frontier_mean, 
           type="l", col=col_vector[1], xlim=c(min_x-0.0005, max_x+0.0001), ylim=c(0, max_y+0.00001), 
           xlab=risk_measure, 
           ylab=expression(paste(mu, "-MC")))
      for(i in 2:length(eff_frontier_list)){
        polygon(c(eff_frontier_list[[i]]$ES95_mc, rev(eff_frontier_list[[i]]$ES95_mc)), 
                c(eff_frontier_list[[i]]$mu_lb, rev(eff_frontier_list[[i]]$mu_ub)), 
                col="grey90", border=F)
      }
      lines(eff_frontier_list[[1]]$frontier_es95, eff_frontier_list[[1]]$frontier_mean, 
            type="l", col=col_vector[1])
      for(i in 2:length(eff_frontier_list)){
        lines(eff_frontier_list[[i]]$ES95_mc, eff_frontier_list[[i]]$mu_mc, col=paste0(col_vector[i], 3))
        lines(eff_frontier_list[[i]]$ES95_mc, eff_frontier_list[[i]]$mu_lb, col=col_vector[i], lty=2)
        lines(eff_frontier_list[[i]]$ES95_mc, eff_frontier_list[[i]]$mu_ub, col=col_vector[i], lty=2)
      }
    }
    else if(risk_measure=="ES99"){
      min_x<-min(min(eff_frontier_list[[1]]$frontier_es99), 
                 sapply(eff_frontier_list[2:3], function(i) min(i$ES99_mc)))
      max_x<-max(max(eff_frontier_list[[1]]$frontier_es95), 
                 sapply(eff_frontier_list[2:3], function(i) max(i$ES99_mc)))
      max_y<-max(sapply(eff_frontier_list[2:3], function(i) max(i$mu_ub)))
      plot(eff_frontier_list[[1]]$frontier_es99, eff_frontier_list[[1]]$frontier_mean, 
           type="l", col=col_vector[1], xlim=c(min_x-0.0005, max_x+0.0001), ylim=c(0, max_y+0.00001), 
           xlab=risk_measure, 
           ylab=expression(paste(mu, "-MC")))
      for(i in 2:length(eff_frontier_list)){
        polygon(c(eff_frontier_list[[i]]$ES99_mc, rev(eff_frontier_list[[i]]$ES99_mc)), 
                c(eff_frontier_list[[i]]$mu_lb, rev(eff_frontier_list[[i]]$mu_ub)), 
                col="grey90", border=F)
      }
      lines(eff_frontier_list[[1]]$frontier_es99, eff_frontier_list[[1]]$frontier_mean, 
            type="l", col=col_vector[1])
      for(i in 2:length(eff_frontier_list)){
        lines(eff_frontier_list[[i]]$ES99_mc, eff_frontier_list[[i]]$mu_mc, col=paste0(col_vector[i], 3))
        lines(eff_frontier_list[[i]]$ES99_mc, eff_frontier_list[[i]]$mu_lb, col=col_vector[i], lty=2)
        lines(eff_frontier_list[[i]]$ES99_mc, eff_frontier_list[[i]]$mu_ub, col=col_vector[i], lty=2)
      }
    }
    else if(risk_measure=="var"){
      min_x<-min(min(eff_frontier_list[[1]]$frontier_sd), 
                 sapply(eff_frontier_list[2:3], function(i) min(i$sd_mc)))
      max_x<-max(max(eff_frontier_list[[1]]$frontier_sd), 
                 sapply(eff_frontier_list[2:3], function(i) max(i$sd_mc)))
      max_y<-max(sapply(eff_frontier_list[2:3], function(i) max(i$mu_ub)))
      plot(eff_frontier_list[[1]]$frontier_sd, eff_frontier_list[[1]]$frontier_mean, 
           type="l", col=col_vector[1], xlim=c(min_x-0.0005, max_x+0.001), ylim=c(0, max_y+0.00001), 
           xlab=risk_measure, 
           ylab=expression(paste(mu, "-MC")))
      for(i in 2:length(eff_frontier_list)){
        polygon(c(eff_frontier_list[[i]]$sd_mc, rev(eff_frontier_list[[i]]$sd_mc)), 
                c(eff_frontier_list[[i]]$mu_lb, rev(eff_frontier_list[[i]]$mu_ub)), 
                col="grey90", border=F)
      }
      lines(eff_frontier_list[[1]]$frontier_sd, eff_frontier_list[[1]]$frontier_mean, 
            type="l", col=col_vector[1])
      for(i in 2:length(eff_frontier_list)){
        lines(eff_frontier_list[[i]]$sd_mc, eff_frontier_list[[i]]$mu_mc, col=paste0(col_vector[i], 3))
        lines(eff_frontier_list[[i]]$sd_mc, eff_frontier_list[[i]]$mu_lb, col=col_vector[i], lty=2)
        lines(eff_frontier_list[[i]]$sd_mc, eff_frontier_list[[i]]$mu_ub, col=col_vector[i], lty=2)
      }
    }
  }
  else{
    if(risk_measure=="ES95"){
      min_x<-min(sapply(eff_frontier_list, function(i) min(i$ES95_mc)))
      max_x<-max(sapply(eff_frontier_list, function(i) max(i$ES95_mc)))
      max_y<-max(sapply(eff_frontier_list, function(i) max(i$mu_ub)))
      plot(eff_frontier_list[[1]]$ES95_mc, eff_frontier_list[[1]]$mu_mc, 
           type="l", col=col_vector[1], xlim=c(min_x-0.0005, max_x+0.0001), ylim=c(0, max_y+0.00001), 
           xlab=risk_measure, 
           ylab=expression(paste(mu, "-MC")))
      for(i in 1:length(eff_frontier_list)){
        polygon(c(eff_frontier_list[[i]]$ES95_mc, rev(eff_frontier_list[[i]]$ES95_mc)), 
                c(eff_frontier_list[[i]]$mu_lb, rev(eff_frontier_list[[i]]$mu_ub)), 
                col="grey90", border=F)
      }
      for(i in 1:length(eff_frontier_list)){
        lines(eff_frontier_list[[i]]$ES95_mc, eff_frontier_list[[i]]$mu_mc, col=paste0(col_vector[i], 3))
        lines(eff_frontier_list[[i]]$ES95_mc, eff_frontier_list[[i]]$mu_lb, col=col_vector[i], lty=2)
        lines(eff_frontier_list[[i]]$ES95_mc, eff_frontier_list[[i]]$mu_ub, col=col_vector[i], lty=2)
      }
    }
    else if(risk_measure=="ES99"){
      min_x<-min(sapply(eff_frontier_list, function(i) min(i$ES99_mc)))
      max_x<-max(sapply(eff_frontier_list, function(i) max(i$ES99_mc)))
      max_y<-max(sapply(eff_frontier_list, function(i) max(i$mu_ub)))
      plot(eff_frontier_list[[1]]$ES99_mc, eff_frontier_list[[1]]$mu_mc, 
           type="l", col=col_vector[1], xlim=c(min_x-0.0005, max_x+0.0001), ylim=c(0, max_y+0.00001), 
           xlab=risk_measure, 
           ylab=expression(paste(mu, "-MC")))
      for(i in 1:length(eff_frontier_list)){
        polygon(c(eff_frontier_list[[i]]$ES99_mc, rev(eff_frontier_list[[i]]$ES99_mc)), 
                c(eff_frontier_list[[i]]$mu_lb, rev(eff_frontier_list[[i]]$mu_ub)), 
                col="grey90", border=F)
      }
      for(i in 1:length(eff_frontier_list)){
        lines(eff_frontier_list[[i]]$ES99_mc, eff_frontier_list[[i]]$mu_mc, col=paste0(col_vector[i], 3))
        lines(eff_frontier_list[[i]]$ES99_mc, eff_frontier_list[[i]]$mu_lb, col=col_vector[i], lty=2)
        lines(eff_frontier_list[[i]]$ES99_mc, eff_frontier_list[[i]]$mu_ub, col=col_vector[i], lty=2)
      }
    }
    else if(risk_measure=="var"){
      min_x<-min(sapply(eff_frontier_list, function(i) min(i$sd_mc)))
      max_x<-max(sapply(eff_frontier_list, function(i) max(i$sd_mc)))
      max_y<-max(sapply(eff_frontier_list, function(i) max(i$mu_ub)))
      plot(eff_frontier_list[[1]]$sd_mc, eff_frontier_list[[1]]$mu_mc, 
           type="l", col=col_vector[1], xlim=c(min_x-0.0005, max_x+0.0001), ylim=c(0, max_y+0.00001), 
           xlab=risk_measure, 
           ylab=expression(paste(mu, "-MC")))
      for(i in 1:length(eff_frontier_list)){
        polygon(c(eff_frontier_list[[i]]$sd_mc, rev(eff_frontier_list[[i]]$sd_mc)), 
                c(eff_frontier_list[[i]]$mu_lb, rev(eff_frontier_list[[i]]$mu_ub)), 
                col="grey90", border=F)
      }
      for(i in 1:length(eff_frontier_list)){
        lines(eff_frontier_list[[i]]$sd_mc, eff_frontier_list[[i]]$mu_mc, col=paste0(col_vector[i], 3))
        lines(eff_frontier_list[[i]]$sd_mc, eff_frontier_list[[i]]$mu_lb, col=col_vector[i], lty=2)
        lines(eff_frontier_list[[i]]$sd_mc, eff_frontier_list[[i]]$mu_ub, col=col_vector[i], lty=2)
      }
    }
  }
}
dev.off()
par(mfrow=c(1,2), mar = c(6.5, 4.1, 3, 2.1))
plot_eff_frontier(eff_frontier_list = list(eff_hist_es95, eff_boot2_es95, eff_vc2_es95), 
                  risk_measure = "ES95", col_vector = c("green2", "red", "blue"))
plot_eff_frontier(eff_frontier_list = list(eff_hist_es95, eff_boot2_es99, eff_vc2_es99), 
                  risk_measure = "ES99", col_vector = c("green2", "red", "blue"))
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom', inset = c(-0.2, 0), legend = c("Historical", "Bootstrap", "Bootstrap-CI", "Vine Copula", "Vine Copula-CI"), 
       col = c("green2","red", "red3", "blue", "blue2"), 
       lty=c(1, 1, 3, 1, 3), lwd = 2, xpd = TRUE, horiz = TRUE, cex = 0.9, seg.len=1, bty = 'n')



