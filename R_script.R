#### Script with alll the relevant content for the print version ####
#Load package running script package.R
#Load data
prices_ss<-read.csv("C:/Users/Pedro/Documents/Studium/WU Dokummente/SoSe 23/Thesis/Scripts/Data/prices_ss.csv")
prices_ss$Date<-as.Date(prices_ss$Date)
#Build empirical loss distribution (scenario set) with log rets and losses in the right tail
scenario.set<-loss.dist(price_df = prices_ss, names = colnames(prices_ss)[-1], 
                        ret_type = "log", neg_rets = T, no_zeroes = F)
#Scenario set with losses in the left tail
scenario.set.p<-loss.dist(price_df = prices_ss, names = colnames(prices_ss)[-1], 
                          ret_type = "log", neg_rets = F, no_zeroes = F)

#### Summary statistics ####
library(xtable)
summary<-summary_stats(scenario.set = scenario.set, neg_rets = T, log_rets = T, alpha = 0.95)
summary

#Kendall's tau correlation
cor_scenario_set<-cor(scenario.set, method="kendall")
cor_scenario_set
cor_scenario_set_lower<-round(get_triangular(cor_matrix = cor_scenario_set, which.tri = "lower"), 2)
rownames(cor_scenario_set_lower)<-colnames(cor_scenario_set_lower)

#Max and min correlations
cor_scenario_set[cor_scenario_set==1]<-0
which(cor_scenario_set == min(cor_scenario_set), arr.ind = TRUE)
max(cor_scenario_set)
min(cor_scenario_set)

# Fitting of marginals via SPGPD ------------------------------------------
fit_double_tail<-spgpd.fit.list(scenario.set = scenario.set, method = "MLE", 
                                double.tail = TRUE)
fit_double_tail_2<-spgpd.fit.list(scenario.set = scenario.set, lower = 0.1, 
                                  upper = 0.9, method = "MLE", double.tail = TRUE)
fit_one_tail<-spgpd.fit.list(scenario.set = scenario.set, upper = 0.9, method = "MLE")

fit_double_tail_kernel<-spd.fit.list(scenario.set = scenario.set, upper = 0.9, 
                                     lower = 0.1, tailfit = "GPD", kernelfit = "epanech", 
                                     type = "mle")

#QQ Plots
qqplots(scenario.set)
#Mean excess plots comparing baseline (n_u=150) and top/bottom 10% percentile
mean.excess.plot.grid(df = scenario.set, gpd.fits = fit_double_tail, 
                      gpd.fits2 = fit_double_tail_2, lower.excesses = T)
mean.excess.plot.grid(df = scenario.set, gpd.fits = fit_double_tail, 
                      gpd.fits2 = fit_double_tail_2, lower.excesses = F)

#Shape plots comparing baseline (n_u=150) and top/bottom 10% percentile
GPD.shape.plot.grid(df = scenario.set, gpd.fits = fit_double_tail, gpd.fits2 = fit_double_tail_2)
GPD.shape.plot.grid(df = scenario.set, gpd.fits = fit_double_tail, 
                    gpd.fits2 = fit_double_tail_2, lower.excesses = T)

#Table of parameters and standard errors
param_tab_1<-fitted.param.tab(gpd.fits = fit_double_tail, tickers = colnames(scenario.set))
param_tab_2<-fitted.param.tab(gpd.fits = fit_double_tail_2, tickers = colnames(scenario.set))
print(xtable(param_tab_1, digits = 3))
print(xtable(param_tab_2, digits = 3))

# Multivariate estimation and simulation -----------------------------------------------------
#Empirical copulae
emp_cop_dtail<-copula.emp(scenario.set = scenario.set, fits = fit_double_tail, method = "SPGPD")
emp_cop_dtail2<-copula.emp(scenario.set = scenario.set, fits = fit_double_tail_2, method = "SPGPD")
emp_cop_onetail<-copula.emp(scenario.set = scenario.set, fits = fit_one_tail, method = "SPGPD")
emp_cop_kernel<-copula.emp(scenario.set = scenario.set, fits = fit_double_tail_kernel, method = "spd")
emp_cop_norm<-apply(scenario.set, MARGIN = 2, function(i) pnorm(q = i, mean = mean(i), sd = sd(i)))

#Vine copulae
start<-Sys.time()
vine_dtail<-vinecop(data = emp_cop_dtail, selcrit = "bic", cores = 12)
vine_dtail2<-vinecop(data = emp_cop_dtail2, selcrit = "bic", cores = 12)
vine_onetail<-vinecop(data = emp_cop_onetail, selcrit = "bic", cores = 12)
vine_kernel<-vinecop(data = emp_cop_kernel, selcrit = "bic", cores = 12)
vine_normal<-vinecop(data = emp_cop_norm, selcrit = "bic", cores=12)
end<-Sys.time()
end-start

#Amount of pairs belonging to a certain family
print(xtable(cbind(pair_cop_families(vine_dtail), pair_cop_families(vine_dtail2))))
pair_cop_families(vine_normal)

#Plots of vine trees
for(i in 1:(ncol(scenario.set)-1)){
  pdf(file = paste0("C:/Users/Pedro/Documents/Studium/WU Dokummente/SoSe 23/Thesis/Thesis draft/tree_", i, ".pdf"),
      width = 5, height = 5)
  plot(vine_dtail, tree=i)
  dev.off()
}

for(i in 1:(ncol(scenario.set)-1)){
  pdf(file = paste0("C:/Users/Pedro/Documents/Studium/WU Dokummente/SoSe 23/Thesis/Thesis draft/tree_9010_", i, ".pdf"),
      width = 5, height = 5)
  plot(vine_dtail2, tree=i)
  dev.off()
}

#Scenario set generation
start<-Sys.time()
set.seed(888)
#Bootstrap
ss_list_boot<-scen.generator.list(copula_df = emp_cop_dtail, 
                                  scenario.set = scenario.set, fits = fit_double_tail, n = 200)
ss_list_boot2<-scen.generator.list(copula_df = emp_cop_dtail2, 
                                  scenario.set = scenario.set, fits = fit_double_tail_2, n = 200)
#Vine Copula
ss_list_vine<-scen.generator.list(vinecop = vine_dtail, 
                                  scenario.set = scenario.set, fits = fit_double_tail, n = 200)
ss_list_vine2<-scen.generator.list(vinecop = vine_dtail2, 
                                 scenario.set = scenario.set, fits = fit_double_tail_2, n = 200)
ss_list_vine_onetail<-scen.generator.list(vinecop = vine_onetail, 
                                 scenario.set = scenario.set, fits = fit_one_tail, n = 200)
ss_list_vine_kernel<-scen.generator.list(vinecop = vine_onetail, 
                                 scenario.set = scenario.set, fits = fit_double_tail_kernel, 
                                 n = 200, method = "spd")
end<-Sys.time()
end-start


#For sensitivity tests: generate samples from multivariate normal
mu_vector<-colMeans(scenario.set)
sigma_matrix<-cov(scenario.set)
ss_list_mvnorm<-lapply(1:200, function(i){
  mvrnorm(nrow(scenario.set), mu = mu_vector, Sigma = sigma_matrix)
  })
ss_list_vine_norm<-scen.generator.list(vinecop = vine_normal, scenario.set = scenario.set, 
                                       n = 200, method = "normal")

# Optimal portfolios and efficiency frontiers -----------------------------
#Optimal portfolios
#Historical scenario set
methods<-c("var", "ES", "ES")
conf_levels<-c(NA, 0.95, 0.99)
opt_ports_historical<-lapply(seq_along(methods), function(i){
  stoch.optim(scenario.set = scenario.set*(-1), method = methods[i], 
              conf_level = conf_levels[i])
})

#Bootstrap
start<-Sys.time()
opt_ports_boot<-lapply(seq_along(methods), function(i){
  stoch.optim.mc(scenario_set_list = ss_list_boot, method = methods[i], 
                 conf_level = conf_levels[i])
})
print(xtable(t(do.call(rbind, opt_ports_boot)), digits=3))

opt_ports_boot2<-lapply(seq_along(methods), function(i){
  stoch.optim.mc(scenario_set_list = ss_list_boot2, method = methods[i], 
                 conf_level = conf_levels[i])
})

#Vine Copula
opt_ports_vc<-lapply(seq_along(methods), function(i){
  stoch.optim.mc(scenario_set_list = ss_list_vine, method = methods[i], 
                 conf_level = conf_levels[i])
})

opt_ports_vc2<-lapply(seq_along(methods), function(i){
  stoch.optim.mc(scenario_set_list = ss_list_vine2, method = methods[i], 
                 conf_level = conf_levels[i])
})
end<-Sys.time()
end-start

opt_ports<-c(opt_ports_historical, opt_ports_boot, opt_ports_boot2, 
             opt_ports_vc, opt_ports_vc2)
print(xtable(t(do.call(rbind, opt_ports)), digits=3))


#Efficiency frontiers
#Historical
eff_hist_var<-eff.frontier(scenario.set = scenario.set*(-1), 
                           mu_sigma_space = "var", conf_level = 0.05)
eff_hist_es95<-eff.frontier(scenario.set = scenario.set*(-1), 
                            mu_sigma_space = "ES", conf_level = 0.05)
eff_hist_es99<-eff.frontier(scenario.set = scenario.set*(-1), 
                            mu_sigma_space = "ES", conf_level = 0.01)

#Bootstrap
eff_boot_var<-eff.frontier.mc.par(scenario_set_list = ss_list_boot, 
                                  mu_sigma_space = "var", conf_level = 0.05)
eff_boot_es95<-eff.frontier.mc.par(scenario_set_list = ss_list_boot, 
                                   mu_sigma_space = "ES", conf_level = 0.05)
eff_boot_es99<-eff.frontier.mc.par(scenario_set_list = ss_list_boot, 
                                   mu_sigma_space = "ES", conf_level = 0.01)

eff_boot2_var<-eff.frontier.mc.par(scenario_set_list = ss_list_boot2, 
                                   mu_sigma_space = "var", conf_level = 0.05)
eff_boot2_es95<-eff.frontier.mc.par(scenario_set_list = ss_list_boot2, 
                                    mu_sigma_space = "ES", conf_level = 0.05)
eff_boot2_es99<-eff.frontier.mc.par(scenario_set_list = ss_list_boot2, 
                                    mu_sigma_space = "ES", conf_level = 0.01)

#Vine Copula
eff_vc_var<-eff.frontier.mc.par(scenario_set_list = ss_list_vine, 
                                mu_sigma_space = "var", conf_level = 0.05)
eff_vc_es95<-eff.frontier.mc.par(scenario_set_list = ss_list_vine,
                                 mu_sigma_space = "ES", conf_level = 0.05)
eff_vc_es99<-eff.frontier.mc.par(scenario_set_list = ss_list_vine, 
                                 mu_sigma_space = "ES", conf_level = 0.01)

eff_vc2_var<-eff.frontier.mc.par(scenario_set_list = ss_list_vine2, 
                                 mu_sigma_space = "var", conf_level = 0.05)
eff_vc2_es95<-eff.frontier.mc.par(scenario_set_list = ss_list_vine2, 
                                  mu_sigma_space = "ES", conf_level = 0.05)
eff_vc2_es99<-eff.frontier.mc.par(scenario_set_list = ss_list_vine2, 
                                  mu_sigma_space = "ES", conf_level = 0.01)

#For sensitivity tests: only ES95 needed
eff_vc_onetail<-eff.frontier.mc.par(scenario_set_list = ss_list_vine_onetail, 
                                    mu_sigma_space = "ES", conf_level = 0.05)
eff_vc_kernel<-eff.frontier.mc.par(scenario_set_list = ss_list_vine_kernel, 
                                    mu_sigma_space = "ES", conf_level = 0.05)
eff_mvnorm<-eff.frontier.mc.par(scenario_set_list = ss_list_mvnorm, 
                                mu_sigma_space = "ES", conf_level = 0.05)
start<-Sys.time()
eff_vc_normal<-eff.frontier.mc.par(scenario_set_list = ss_list_vine_norm, 
                                   mu_sigma_space = "ES", conf_level = 0.05)
end<-Sys.time()
end-start


#Plots of efficiency frontiers
#Baseline plots
dev.off()
par(mfrow=c(1,2), mar = c(6.5, 4.1, 3, 2.1))
plot_eff_frontier(eff_frontier_list = list(eff_hist_es95, eff_boot2_es95, eff_vc2_es95), 
                  risk_measure = "ES95", col_vector = c("green2", "red", "blue"), hist_ef_included = T)
plot_eff_frontier(eff_frontier_list = list(eff_hist_es95, eff_boot2_es99, eff_vc2_es99), 
                  risk_measure = "ES99", col_vector = c("green2", "red", "blue"))
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom', inset = c(-0.2, 0), legend = c("Historical", "Bootstrap", "Bootstrap CI", "Vine Copula", "Vine Copula CI"), 
       col = c("green2","red", "red3", "blue", "blue2"), 
       lty=c(1, 1, 3, 1, 3), lwd = 2, xpd = TRUE, horiz = TRUE, cex = 0.9, seg.len=1, bty = 'n')


#Different optimization methods
dev.off()
min_x<-min(c(eff_hist_var$frontier_es95, eff_hist_es95$frontier_es95))
max_x<-max(c(eff_hist_var$frontier_es95, eff_hist_es95$frontier_es95))
max_y<-max(c(eff_hist_var$frontier_mean, eff_hist_es95$frontier_mean))
par(mfrow=c(2,3), mar = c(6.5, 4.1, 3, 2.1))
plot(eff_hist_var$frontier_es95, eff_hist_var$frontier_mean, 
     type="l", col="green", xlim=c(min_x-0.0005, max_x+0.0001), ylim=c(0, max_y+0.00001), 
     xlab="ES95", main="Historical",
     ylab=expression(paste(mu, "-MC")))
lines(eff_hist_es95$frontier_es95, eff_hist_es95$frontier_mean, type="l", col="red")
plot_eff_frontier(eff_frontier_list = list(eff_boot2_var, eff_boot2_es95), 
                  risk_measure = "ES95", col_vector = c("green", "red"), 
                  hist_ef_included = F, title = "Bootstrap")
plot_eff_frontier(eff_frontier_list = list(eff_vc2_var, eff_vc2_es95), 
                  risk_measure = "ES95", col_vector = c("green", "red"),
                  hist_ef_included = F, title = "Vine Copula")
min_x<-min(c(eff_hist_var$frontier_es99, eff_hist_es95$frontier_es99))
max_x<-max(c(eff_hist_var$frontier_es99, eff_hist_es95$frontier_es99))
plot(eff_hist_var$frontier_es99, eff_hist_var$frontier_mean, 
     type="l", col="green", xlim=c(min_x-0.0005, max_x+0.0001), ylim=c(0, max_y+0.00001), 
     xlab="ES99", 
     ylab=expression(paste(mu, "-MC")))
lines(eff_hist_es95$frontier_es99, eff_hist_es95$frontier_mean, type="l", col="red")
plot_eff_frontier(eff_frontier_list = list(eff_boot2_var, eff_boot2_es99), 
                  risk_measure = "ES99", col_vector = c("green", "red"), hist_ef_included = F)
plot_eff_frontier(eff_frontier_list = list(eff_vc2_var, eff_vc2_es99), 
                  risk_measure = "ES99", col_vector = c("green", "red"), hist_ef_included = F)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom', legend = c("Mean-Variance", "Mean-Variance-CI", "ES", "ES-CI"), 
       col = c("green", "green3", "red", "red3"), 
       lty=c(1, 3, 1, 3), lwd = 2, xpd = TRUE, horiz = TRUE, cex = 0.9, seg.len=3, x.intersp = 0.25, text.width = 0.2, bty = 'n')

#Difference in threshold selection
dev.off()
par(mfrow=c(1,2), mar = c(6.5, 4.1, 3, 2.1))
plot_eff_frontier(eff_frontier_list = list(eff_hist_es95, eff_boot_es95, eff_boot2_es95), 
                  risk_measure = "ES95", col_vector = c("green", "red", "blue"), 
                  hist_ef_included = T, title = "Bootstrap")
plot_eff_frontier(eff_frontier_list = list(eff_hist_es95, eff_vc_es95, eff_vc2_es95), 
                  risk_measure = "ES95", col_vector = c("green", "red", "blue"), 
                  hist_ef_included = T, title = "Vine Copula")
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom', legend = c("Historical", expression(paste(n[u], "=150")), expression(paste(n[u], "=150 CI")),
                                                "10th/90th quantile", "10th/90th quantile CI"), 
       col = c("green","red", "red3", "blue", "blue2"), 
       lty=c(1, 1, 3, 1, 3), lwd = 2, xpd = TRUE, horiz = TRUE, cex = 0.9, seg.len=2, bty = 'n')

#Difference in selection of marginals (only using VC)
dev.off()
par(mfrow=c(1,2), mar = c(6.5, 4.1, 3, 2.1))
plot_eff_frontier(eff_frontier_list = list(eff_hist_es95, eff_vc2_es95, eff_vc_kernel, eff_vc_onetail, eff_mvnorm), 
                  risk_measure = "ES95", col_vector = c("orange", "green", "red", "blue", "magenta"), 
                  hist_ef_included = T, ci_background = F, title = "Different marginals")
plot_eff_frontier(eff_frontier_list = list(eff_vc2_es95, eff_mvnorm, eff_vc_onetail), 
                  risk_measure = "ES95", col_vector = c("green", "magenta", "blue"), 
                  hist_ef_included = F, ci_background = T, title = "Baseline vs. One tailed vs. multivariate normal")
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0.5, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom', legend = c("Historical", "SPD\nEmpirical Center", "SPD\nKernel Center", 
                            "SPD\nOnly loss tail", "Multivar. normal"), 
       col = c("orange", "green", "red", "blue", "magenta"), 
       lty=1, lwd = 2, xpd = TRUE, horiz = TRUE, cex = 0.9, seg.len=3, x.intersp = 2, bty = 'n')
