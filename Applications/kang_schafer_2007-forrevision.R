set.seed(12398123)
expit <- function(a) {exp(a) / (1 + exp(a))}
wtd_mean <- function(x, z, w){ sum(x * z * w) / sum(z * w) }
n_obs <- 200
n_sim <- 100



sim_res <- matrix(NA, nrow = n_sim, ncol = 6)
mean_bal_res1 <- matrix(NA, nrow = n_sim, ncol = 6)
mean_bal_res2 <- matrix(NA, nrow = n_sim, ncol = 6)
mean_bal_res3 <- matrix(NA, nrow = n_sim, ncol = 6)
mean_bal_res4 <- matrix(NA, nrow = n_sim, ncol = 6)

var_bal_res1 <- matrix(NA, nrow = n_sim, ncol = 6)
var_bal_res2 <- matrix(NA, nrow = n_sim, ncol = 6)
var_bal_res3 <- matrix(NA, nrow = n_sim, ncol = 6)
var_bal_res4 <- matrix(NA, nrow = n_sim, ncol = 6)

for(s in 1:n_sim){
  z1 <- rnorm(n_obs)
  z2 <- rnorm(n_obs)
  z3 <- rnorm(n_obs)
  z4 <- rnorm(n_obs)
  
  Z <- cbind(z1,z2,z3,z4)
  
  errs <- rnorm(n_obs)
  
  Yi <- 210 + 27.4 * z1 + 13.7 * z2 + 13.7 * z3 + 13.7 * z4 + errs
  prob_i <-  expit(-z1 + 0.5 * z2 - 0.25 * z3 - 0.1 * z4)
  
  ti <- rbinom(n_obs, 1, prob_i)
  
  X1 <- exp(z1 / 2)
  X2 <- z2 / (1 + exp(z1)) + 10
  X3 <- ((z1 * z3 /25) + 0.6)^3
  X4 <- (z2 + z4 + 20)^2
  
  X <- cbind(X1,X2,X3,X4)
  
  Xs <- cbind(scale(X1),
              scale(X2),
              scale(X3),
              scale(X4))
  
  dat <- data.frame('Y' = Yi,  'A' = ti, X)
  
  gpps_snp <- gpbalancer::gpbal(Xs, ti,
                                cov_function = gpbalancer::sqexp_poly,
                                init_theta = c(1,1))
  gpps_ard <- gpbalancer::gpbal(Xs, ti,
                                cov_function = gpbalancer::sqexp,
                                init_theta = c(1,1,1,1,1))
  cbps <- CBPS::CBPS('A ~ X1 + X2 + X3 + X4', data = dat, ATT = 0)
  bart <- BART::mc.pbart(X, ti, seed = s+2019)
  gbm <- twang::ps(A ~ X1 + X2 + X3 + X4, data = dat, verbose = F, estimand = 'ATE')
  
  true_wts <- ifelse(ti==1, 1/prob_i, 1/(1-prob_i))
  gpbal1_wts <- ifelse(ti == 1, 1/gpps_snp$ps, 1/(1-gpps_snp$ps))
  gpbal2_wts <- ifelse(ti == 1, 1/gpps_ard$ps, 1/(1-gpps_ard$ps))
  cbps_wts <- ifelse(ti == 1, 1/cbps$fitted.values, 1/(1-cbps$fitted.values))
  bart_wts <- ifelse(ti == 1, 1/bart$prob.train.mean, 1/(1-bart$prob.train.mean))
  gbm_wts <- ifelse(ti == 1, 1/gbm$ps[,2], 1/(1-gbm$ps[,2]))
  
  sim_res[s,1] <- wtd_mean(Yi, ti == 1, true_wts) - wtd_mean(Yi, ti == 0, true_wts)
  sim_res[s,2] <- wtd_mean(Yi, ti == 1, gpbal1_wts) - wtd_mean(Yi, ti == 0, gpbal1_wts)
  sim_res[s,3] <- wtd_mean(Yi, ti == 1, gpbal2_wts) - wtd_mean(Yi, ti == 0, gpbal2_wts)
  sim_res[s,4] <- wtd_mean(Yi, ti == 1, cbps_wts) - wtd_mean(Yi, ti == 0, cbps_wts)
  sim_res[s,5] <- wtd_mean(Yi, ti == 1, bart_wts) - wtd_mean(Yi, ti == 0, bart_wts)
  sim_res[s,6] <- wtd_mean(Yi, ti == 1, gbm_wts) - wtd_mean(Yi, ti == 0, gbm_wts)

  bal_true <- gpbalancer::bal_table(dat, 3:6, ti, wts = true_wts)
  bal_gpb1 <- gpbalancer::bal_table(dat, 3:6, ti, wts = gpbal1_wts)
  bal_gpb2 <- gpbalancer::bal_table(dat, 3:6, ti, wts = gpbal2_wts)
  bal_cbps <- gpbalancer::bal_table(dat, 3:6, ti, wts = cbps_wts)
  bal_bart <- gpbalancer::bal_table(dat, 3:6, ti, wts = bart_wts)
  bal_gbmt <- gpbalancer::bal_table(dat, 3:6, ti, wts = gbm_wts)
  
  mean_bal_res1[s, ] <- c(bal_true[1,7], bal_gpb1[1,7], bal_gpb2[1,7], bal_cbps[1,7], bal_bart[1,7], bal_gbmt[1,7])
  mean_bal_res2[s, ] <- c(bal_true[2,7], bal_gpb1[2,7], bal_gpb2[2,7], bal_cbps[2,7], bal_bart[2,7], bal_gbmt[2,7])
  mean_bal_res3[s, ] <- c(bal_true[3,7], bal_gpb1[3,7], bal_gpb2[3,7], bal_cbps[3,7], bal_bart[3,7], bal_gbmt[3,7])
  mean_bal_res4[s, ] <- c(bal_true[4,7], bal_gpb1[4,7], bal_gpb2[4,7], bal_cbps[4,7], bal_bart[4,7], bal_gbmt[4,7])
  
  var_bal_res1[s, ] <- c(bal_true[1,8], bal_gpb1[1,8], bal_gpb2[1,8], bal_cbps[1,8], bal_bart[1,8], bal_gbmt[1,8])
  var_bal_res2[s, ] <- c(bal_true[2,8], bal_gpb1[2,8], bal_gpb2[2,8], bal_cbps[2,8], bal_bart[2,8], bal_gbmt[2,8])
  var_bal_res3[s, ] <- c(bal_true[3,8], bal_gpb1[3,8], bal_gpb2[3,8], bal_cbps[3,8], bal_bart[3,8], bal_gbmt[3,8])
  var_bal_res4[s, ] <- c(bal_true[4,8], bal_gpb1[4,8], bal_gpb2[4,8], bal_cbps[4,8], bal_bart[4,8], bal_gbmt[4,8])
  
  message(paste(s,paste(rep('.',80),collapse = ''), sep=""), appendLF = T)
  print(round(apply(sim_res, 2, mean, na.rm=T), 3))
}



method_names <- c('True Propensity Score',
                  'GPPS - SE+NP',
                  'GPPS - SE-ARD',
                  'CBPS',
                  'BART', 
                  'GBM - twang')

colnames(sim_res) <- method_names


xtable::xtable(cbind(apply(sim_res - 0, 2, mean, na.rm=T),
                     apply(abs(sim_res - 0), 2, mean, na.rm=T),
                     sqrt(apply((sim_res - 0)^2, 2, mean, na.rm=T))), digits=3)



xtable::xtable(cbind(apply(abs(mean_bal_res1), 2, mean),
                     apply(abs(mean_bal_res2), 2, mean),
                     apply(abs(mean_bal_res3), 2, mean),
                     apply(abs(mean_bal_res4), 2, mean),
                     apply(abs(var_bal_res1), 2, mean),
                     apply(abs(var_bal_res2), 2, mean),
                     apply(abs(var_bal_res3), 2, mean),
                     apply(abs(var_bal_res4), 2, mean)), digits = 3)

pdf('~/Documents/manuscripts/2019-GPPS-revision/ks2007sim.pdf', height = 3, width = 8)
par(mfrow=c(1,1), mar = c(4,10,2,2))
plot(0, xlim = range(sim_res), ylim = c(0.5,6.5), axes = F,
     pch = 19, col = rgb(0,0,0,0),
     xlab = 'Estimated ATE',
     ylab = '')
abline(v = 0, col = rgb(1,0,0,0.5), lwd = 5)
boxplot(sim_res[,6:1], horizontal=T, axes = F, add = T, col = 'gray', pch = 19)
axis(1); axis(2, las = 2, at = 6:1, labels =  method_names)
dev.off()


par(mfrow=c(4,1))
plot(0, xlim = c(0, max(abs(mean_bal_res1))), ylim = c(0.5,6.5), axes = F,
     pch = 19, col = rgb(0,0,0,0),
     xlab = 'Absolute Standardized Mean Difference',
     ylab = '', main = expression(paste('Covariate ', X[1])))
abline(v = 0, col = rgb(0,0,0,1), lwd = 0)
abline(v = 0.2, col = rgb(1,0,0,0.5), lwd = 5)
boxplot(abs(mean_bal_res1[,6:1]), horizontal=T, axes = F, add = T, col = 'gray', pch = 19)
axis(1); axis(2, las = 2, at = 6:1, labels =  method_names)

plot(0, xlim = c(0, max(abs(mean_bal_res1))), ylim = c(0.5,6.5), axes = F,
     pch = 19, col = rgb(0,0,0,0),
     xlab = 'Absolute Standardized Mean Difference',
     ylab = '', main = expression(paste('Covariate ', X[2])))
abline(v = 0, col = rgb(0,0,0,1), lwd = 0)
abline(v = 0.2, col = rgb(1,0,0,0.5), lwd = 5)
boxplot(abs(mean_bal_res2[,6:1]), horizontal=T, axes = F, add = T, col = 'gray', pch = 19)
axis(1); axis(2, las = 2, at = 6:1, labels =  method_names)





# pdf('Documents/manuscripts/2019-GPPS-revision/gppaper-revision/bal_plot.pdf',
#     height = 5, width= 8)
# par(mfrow=c(1,1), mar = c(5,12,4,2)+0.1)
# apply(bal_res, 2, mean, na.rm = T)
# boxplot(bal_res[,9:1], horizontal = T, axes = F,
#         main = 'Average abs(Standardized Difference in Means)')
# axis(1)
# axis(2, at = 1:9, method_names[9:1], las = 2)
# abline(v = 0.1, lty = 3, lwd = 3, col = rgb(0, 0,0.5, 0.5))
# abline(v = 0.2, lty = 1, lwd = 3, col = rgb(0.5, 0,0, 0.5))
# # dev.off()
# 
# IPW_res <- cbind(apply(sim_res - 210, 2, mean, na.rm=T),
#                  sqrt(apply((sim_res - 210)^2, 2, mean, na.rm=T)))
# rownames(IPW_res) <- method_names
# 
# WLS_res <- cbind(apply(sim_res2 - 210, 2, mean, na.rm=T),
#                  sqrt(apply((sim_res2 - 210)^2, 2, mean, na.rm=T)))
# rownames(WLS_res) <- method_names
# 
# DR_res <- cbind(apply(sim_res3 - 210, 2, mean, na.rm=T),
#                 sqrt(apply((sim_res3 - 210)^2, 2, mean, na.rm=T)))
# rownames(DR_res) <- method_names
# 
# xtable::xtable(cbind(IPW_res, WLS_res, DR_res), digits = c(0, rep(3, 6)))
