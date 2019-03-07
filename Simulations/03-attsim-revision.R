#-------------------------------------------------------------------------------
# Number of sims to run
set.seed(89274)
# n_sims <- 1000
n_sims <- 1000

setwd('~/Documents/GitHub/paperOptBalGPPS/Simulations/att-simresults/')

lm_ps <- function(Y, X, wts, true_val = NULL){
  W <- diag(wts)
  invXtWX <- solve(t(X) %*% W %*% X)
  betas <- invXtWX %*% t(X) %*% W %*% Y

  Yhat <- X %*% betas
  resids <- Y - Yhat
  sighat <- as.double(sum(resids^2) / (length(Y) - 1))

  # varmat <- sighat * invXtWX %*% t(X) %*% W %*% W %*% X %*% invXtWX
  varmat <- invXtWX %*% t(X) %*% W %*% diag(as.vector(resids)^2) %*% W %*% X %*% invXtWX
  std_errs <- sqrt(diag(varmat))

  low_int <- betas - 1.96 * std_errs
  upp_int <- betas + 1.96 * std_errs

  res <- cbind(betas, std_errs, low_int, upp_int)
  colnames(res) <- c('coef', 'stderrs', 'low95', 'upp95')

  if(!is.null(true_val)){
    cover_test <- res[2,3] < true_val & res[2,4] > true_val
    return(list('ests' = data.frame(res),
                'covers' = cover_test))
  } else{
    return(list('ests' = res))
  }
}


res_mat_em <- matrix(NA, nrow=n_sims, ncol=12)
res_mat_lin <- matrix(NA, nrow=n_sims, ncol=12)
res_mat_non <- matrix(NA, nrow=n_sims, ncol=12)
bal_mats1 <- matrix(NA, nrow=n_sims, ncol=11)
bal_mats2 <- matrix(NA, nrow=n_sims, ncol=11)
bal_mats3 <- matrix(NA, nrow=n_sims, ncol=11)
bal_mats4 <- matrix(NA, nrow=n_sims, ncol=11)
bal_mats5 <- matrix(NA, nrow=n_sims, ncol=11)

n_obs <- 500

for(i in 1:n_sims){
  X1 <- rnorm(n_obs)
  X2 <- rnorm(n_obs)
  X3 <- rbinom(n_obs, 1, 0.3)
  X4 <- rbinom(n_obs, 1, prob = pnorm(X1))
  X5 <- rbinom(n_obs, 1, prob = pnorm(X2))
  X <- data.frame(X1, X2, X3, X4, X5)

  p <- pnorm(0.5 * X1 + 0.25 * X2 + 0.1 * X1 * X2 * X3 + 0.05 * X2 * X5 + 0.025 * X4)
  TA <- rbinom(n_obs, 1, p)

  dset <- data.frame(TA = TA,
                     X1 = X1,
                     X2 = X2,
                     X3 = X3,
                     X4 = X4,
                     X5 = X5)

  YT_em <- 5 * X1^2 + X1 * X3 - 4 * X2 + 50 * X5 + 10 * X1 - 3 * X2^3 + rnorm(n_obs)
  YC_em <- 5 * X1^2 + X1 * X3 - 4 * X2 + 50 * X5 + rnorm(n_obs)
  YO_em <- ifelse(TA==1, YT_em, YC_em)

  YT_lin <- 5 * X1^2 + X1 * X3 - 4 * X2 + 50 * X5 + 10 + rnorm(n_obs)
  YC_lin <- 5 * X1^2 + X1 * X3 - 4 * X2 + 50 * X5 + rnorm(n_obs)
  YO_lin <- ifelse(TA==1, YT_lin, YC_lin)

  # YT_non <- exp(scale(X[,1])) + 4 * scale(X[,1]) + 3 + rnorm(n_obs, sd = 0.25)
  # YC_non <- - scale(X[,1])^2 - exp(scale(X[,1])) + rnorm(n_obs, sd = 0.25)
  # YO_non <- TA * YT_non + (1-TA)*YC_non

  Xscaled <- cbind(scale(X[,1]),
                   scale(X[,2]),
                   ifelse(X[,3], 1, -1),
                   ifelse(X[,4], 1, -1),
                   ifelse(X[,5], 1, -1))

  est_ps_glm1 <- fitted(glm('TA ~ X1 + X2 + X3 + X4 + X5', data=dset, family = 'binomial'))
  est_ps_glm2 <- fitted(glm('TA ~ poly(X1, 2) + poly(X2, 2) + X3 + X4 + X5', data=dset, family = 'binomial'))
  est_ps_gp1 <- paperOptBalGPPS::gpbal(Xscaled, TA,
                                  paperOptBalGPPS::sqexp_poly,
                                  c(1,1),
                                  wts_vers = 'ATT')
  est_ps_gp2 <- paperOptBalGPPS::gpbal(Xscaled, TA,
                                  paperOptBalGPPS::sqexp_poly2,
                                  c(1,1),
                                  wts_vers = 'ATT')
  est_ps_gp3 <- paperOptBalGPPS::gpbal(Xscaled, TA,
                                  paperOptBalGPPS::sqexp_ard,
                                  c(1, 1, 0.5, 0.5, 0.5),
                                  wts_vers = 'ATT')

  bart_train <- BART::mc.pbart(Xscaled, TA, seed = i+2018)
  est_ps_bart <- bart_train$prob.train.mean
  capture.output(est_ps_cbps <- as.vector(fitted(CBPS::CBPS('TA ~ poly(X1, 2) + poly(X2, 2) + X3 + X4 + X5',
                                                            data=dset, family = 'binomial'))))
  est_ps_gbm1 <- twang::ps(TA ~ X1 + X2 + X3 + X4 + X5, data=dset,
                           estimand = 'ATT', verbose=F,
                           stop.method = 'ks.mean')$ps[,1]
  est_ps_gbm2 <- twang::ps(TA ~ X1 + X2 + X3 + X4 + X5, data=dset,
                           estimand = 'ATT', verbose=F,
                           stop.method = 'es.mean')$ps[,1]
  est_ps_gbm3 <- twang::ps(TA ~ X1 + X2 + X3 + X4 + X5, data=dset,
                           estimand = 'ATT', verbose=F,
                           stop.method = 'es.max')$ps[,1]

  wts_true <- ifelse(TA==1, 1, p/(1-p))
  wts_glm1 <- ifelse(TA==1, 1, est_ps_glm1/(1-est_ps_glm1))
  wts_glm2 <- ifelse(TA==1, 1, est_ps_glm2/(1-est_ps_glm2))
  wts_gp1 <- ifelse(TA==1, 1, est_ps_gp1$ps/(1-est_ps_gp1$ps))
  wts_gp2 <- ifelse(TA==1, 1, est_ps_gp2$ps/(1-est_ps_gp2$ps))
  wts_gp3 <- ifelse(TA==1, 1, est_ps_gp3$ps/(1-est_ps_gp3$ps))
  wts_bart <- ifelse(TA==1, 1, est_ps_bart/(1-est_ps_bart))
  wts_cbps <- ifelse(TA==1, 1, est_ps_cbps/(1-est_ps_cbps))
  wts_gbm1 <- ifelse(TA==1, 1, est_ps_gbm1/(1-est_ps_gbm1))
  wts_gbm2 <- ifelse(TA==1, 1, est_ps_gbm2/(1-est_ps_gbm2))
  wts_gbm3 <- ifelse(TA==1, 1, est_ps_gbm3/(1-est_ps_gbm3))

  bal_mats1[i, 1] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_true)[1, 7]
  bal_mats1[i, 2] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_gp1)[1, 7]
  bal_mats1[i, 3] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_gp2)[1, 7]
  bal_mats1[i, 4] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_gp3)[1, 7]
  bal_mats1[i, 5] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_bart)[1, 7]
  bal_mats1[i, 6] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_gbm1)[1, 7]
  bal_mats1[i, 7] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_gbm2)[1, 7]
  bal_mats1[i, 8] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_gbm3)[1, 7]
  bal_mats1[i, 9] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_cbps)[1, 7]
  bal_mats1[i, 10] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_glm2)[1, 7]
  bal_mats1[i, 11] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_glm1)[1, 7]


  bal_mats2[i, 1] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_true)[2, 7]
  bal_mats2[i, 2] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_gp1)[2, 7]
  bal_mats2[i, 3] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_gp2)[2, 7]
  bal_mats2[i, 4] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_gp3)[2, 7]
  bal_mats2[i, 5] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_bart)[2, 7]
  bal_mats2[i, 6] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_gbm1)[2, 7]
  bal_mats2[i, 7] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_gbm2)[2, 7]
  bal_mats2[i, 8] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_gbm3)[2, 7]
  bal_mats2[i, 9] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_cbps)[2, 7]
  bal_mats2[i, 10] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_glm2)[2, 7]
  bal_mats2[i, 11] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_glm1)[2, 7]


  bal_mats3[i, 1] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_true)[5, 7]
  bal_mats3[i, 2] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_gp1)[5, 7]
  bal_mats3[i, 3] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_gp2)[5, 7]
  bal_mats3[i, 4] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_gp3)[5, 7]
  bal_mats3[i, 5] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_bart)[5, 7]
  bal_mats3[i, 6] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_gbm1)[5, 7]
  bal_mats3[i, 7] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_gbm2)[5, 7]
  bal_mats3[i, 8] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_gbm3)[5, 7]
  bal_mats3[i, 9] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_cbps)[5, 7]
  bal_mats3[i, 10] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_glm2)[5, 7]
  bal_mats3[i, 11] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_glm1)[5, 7]

  bal_mats4[i, 1] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_true)[8, 7]
  bal_mats4[i, 2] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_gp1)[8, 7]
  bal_mats4[i, 3] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_gp2)[8, 7]
  bal_mats4[i, 4] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_gp3)[8, 7]
  bal_mats4[i, 5] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_bart)[8, 7]
  bal_mats4[i, 6] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_gbm1)[8, 7]
  bal_mats4[i, 7] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_gbm2)[8, 7]
  bal_mats4[i, 8] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_gbm3)[8, 7]
  bal_mats4[i, 9] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_cbps)[8, 7]
  bal_mats4[i, 10] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_glm2)[8, 7]
  bal_mats4[i, 11] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_glm1)[8, 7]

  bal_mats5[i, 1] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_true)[11, 7]
  bal_mats5[i, 2] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_gp1)[11, 7]
  bal_mats5[i, 3] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_gp2)[11, 7]
  bal_mats5[i, 4] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_gp3)[11, 7]
  bal_mats5[i, 5] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_bart)[11, 7]
  bal_mats5[i, 6] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_gbm1)[11, 7]
  bal_mats5[i, 7] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_gbm2)[11, 7]
  bal_mats5[i, 8] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_gbm3)[11, 7]
  bal_mats5[i, 9] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_cbps)[11, 7]
  bal_mats5[i, 10] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_glm2)[11, 7]
  bal_mats5[i, 11] <- paperOptBalGPPS::bal_table(X, 1:5, TA==1, wts=wts_glm1)[11, 7]

  res_mat_em[i, 1] <- mean((YT_em - YC_em)[TA==1])
  res_mat_em[i, 2] <- lm_ps(YO_em, cbind(1, TA), wts = wts_true)$ests[2,1]
  res_mat_em[i, 3] <- lm_ps(YO_em, cbind(1, TA), wts = wts_gp1)$ests[2,1]
  res_mat_em[i, 4] <- lm_ps(YO_em, cbind(1, TA), wts = wts_gp2)$ests[2,1]
  res_mat_em[i, 5] <- lm_ps(YO_em, cbind(1, TA), wts = wts_gp3)$ests[2,1]
  res_mat_em[i, 6] <- lm_ps(YO_em, cbind(1, TA), wts = wts_bart)$ests[2,1]
  res_mat_em[i, 7] <- lm_ps(YO_em, cbind(1, TA), wts = wts_gbm1)$ests[2,1]
  res_mat_em[i, 8] <- lm_ps(YO_em, cbind(1, TA), wts = wts_gbm2)$ests[2,1]
  res_mat_em[i, 9] <- lm_ps(YO_em, cbind(1, TA), wts = wts_gbm3)$ests[2,1]
  res_mat_em[i, 10] <- lm_ps(YO_em, cbind(1, TA), wts = wts_cbps)$ests[2,1]
  res_mat_em[i, 11] <- lm_ps(YO_em, cbind(1, TA), wts = wts_glm2)$ests[2,1]
  res_mat_em[i, 12] <- lm_ps(YO_em, cbind(1, TA), wts = wts_glm1)$ests[2,1]

  res_mat_lin[i, 1] <- mean((YT_lin - YC_lin)[TA==1])
  res_mat_lin[i, 2] <- lm_ps(YO_lin, cbind(1, TA), wts = wts_true)$ests[2,1]
  res_mat_lin[i, 3] <- lm_ps(YO_lin, cbind(1, TA), wts = wts_gp1)$ests[2,1]
  res_mat_lin[i, 4] <- lm_ps(YO_lin, cbind(1, TA), wts = wts_gp2)$ests[2,1]
  res_mat_lin[i, 5] <- lm_ps(YO_lin, cbind(1, TA), wts = wts_gp3)$ests[2,1]
  res_mat_lin[i, 6] <- lm_ps(YO_lin, cbind(1, TA), wts = wts_bart)$ests[2,1]
  res_mat_lin[i, 7] <- lm_ps(YO_lin, cbind(1, TA), wts = wts_gbm1)$ests[2,1]
  res_mat_lin[i, 8] <- lm_ps(YO_lin, cbind(1, TA), wts = wts_gbm2)$ests[2,1]
  res_mat_lin[i, 9] <- lm_ps(YO_lin, cbind(1, TA), wts = wts_gbm3)$ests[2,1]
  res_mat_lin[i, 10] <- lm_ps(YO_lin, cbind(1, TA), wts = wts_cbps)$ests[2,1]
  res_mat_lin[i, 11] <- lm_ps(YO_lin, cbind(1, TA), wts = wts_glm2)$ests[2,1]
  res_mat_lin[i, 12] <- lm_ps(YO_lin, cbind(1, TA), wts = wts_glm1)$ests[2,1]

  message(paste(i,':', sep=''), appendLF = 'F')
}

results_data <- list('Cov1_Balance' = bal_mats1,
                     'Cov2_Balance' = bal_mats2,
                     'Cov3_Balance' = bal_mats3,
                     'Cov4_Balance' = bal_mats4,
                     'Cov5_Balance' = bal_mats5,
                     'LinearResults' = res_mat_lin,
                     'EffModResults' = res_mat_em)

bias_em_mat <- res_mat_em - mean(res_mat_em[,1])
bias_lin_mat <- res_mat_lin - 10
# bias_non_mat <- res_mat_non - true_non

thresh <- 0.2
mean_bal1 <- apply(abs(bal_mats1) < thresh, 2, mean, na.rm=T)
mean_bal2 <- apply(abs(bal_mats2) < thresh, 2, mean, na.rm=T)
mean_balb <- apply(abs(bal_mats1) < thresh & abs(bal_mats2) < thresh & abs(bal_mats3) < thresh & abs(bal_mats4) < thresh & abs(bal_mats5) < thresh, 2, mean, na.rm=T)

bias_lin_mean <- apply(bias_lin_mat, 2, mean, na.rm=T)
bias_lin_abs <- apply(abs(bias_lin_mat), 2, mean, na.rm=T)
bias_lin_sd <- apply(bias_lin_mat, 2, sd, na.rm=T)
bias_lin_mse <- apply(bias_lin_mat^2, 2, mean, na.rm=T)

bias_em_mean <- apply(bias_em_mat, 2, mean, na.rm=T)
bias_em_abs <- apply(abs(bias_em_mat), 2, mean, na.rm=T)
bias_em_sd <- apply(bias_em_mat, 2, sd, na.rm=T)
bias_em_mse <- apply(bias_em_mat^2, 2, mean, na.rm=T)

# bias_non_mean <- apply(bias_non_mat, 2, mean, na.rm=T)
# bias_non_abs <- apply(abs(bias_non_mat), 2, mean, na.rm=T)
# bias_non_sd <- apply(bias_non_mat, 2, sd, na.rm=T)
# bias_non_mse <- apply(bias_non_mat^2, 2, mean, na.rm=T)

outro_lin <- rbind(c(0,mean_bal1), c(0,mean_bal2), c(0,mean_balb),
                   bias_lin_mean, bias_lin_abs, bias_lin_sd, bias_lin_mse)

outro_em <- rbind(c(0,mean_bal1), c(0,mean_bal2), c(0,mean_balb),
                  bias_em_mean, bias_em_abs, bias_em_sd, bias_em_mse)
# outro_non <- rbind(mean_bal1, mean_bal2, mean_balb,
#                   bias_non_mean, bias_non_abs, bias_non_sd, bias_non_mse)

colnames(outro_lin) <- c('TRUTH', "TRUEPS",
                         'OBGPPS:NPSE1:EP', 'OBGPPS:NPSE2:EP', 'OBGPPS:SE:EP',
                         'BART', 'GBM:KS.MEAN', 'GBM:ES.MEAN', 'GBM:ES.MAX',
                         'CBPS', 'GLM:POLY2', 'GLM:POLY1')
colnames(outro_em) <- c('TRUTH', "TRUEPS",
                        'OBGPPS:NPSE1:EP', 'OBGPPS:NPSE2:EP', 'OBGPPS:SE:EP',
                        'BART', 'GBM:KS.MEAN', 'GBM:ES.MEAN', 'GBM:ES.MAX',
                        'CBPS', 'GLM:POLY2', 'GLM:POLY1')
print(t(outro_lin)[2:12,])
print(t(outro_em)[2:12,])

saveRDS(results_data, paste(Sys.Date(), '-attsim-multivariate-results-revision.rds', sep=''))

