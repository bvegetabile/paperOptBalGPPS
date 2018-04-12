setwd('~/git/paperOptBalGPPS/Simulations/att-simresults/')

#-------------------------------------------------------------------------------
# Number of sims to run
set.seed(89274)
# n_sims <- 1000
n_sims <- 1
#-------------------------------------------------------------------------------

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

offset <- 1
m1_t <- -0.5 + offset
s1_t <- 1.5
m2_t <- 2.5 + offset
s2_t <- 0.5
p1_t <- 0.25

m1_c <- -1
s1_c <- 0.5
m2_c <- 2
s2_c <- 3
p1_c <- 0.5

prob_t <- 0.2

n_t <- 100
n_c <- n_t * (1-prob_t) / prob_t

m_t <- p1_t * m1_t + (1-p1_t) * m2_t
v_t <- p1_t * ((m1_t - m_t)^2 + s1_t^2) + (1 - p1_t) * ((m2_t - m_t)^2 + s2_t^2)

true_att_em <- v_t + m_t^2 + 2 - m_t
true_att_lin <- m_t + 3 - m_t

testexes <- seq(-10, 10, length.out = 1000)
densXt <- p1_t * dnorm(testexes, m1_t, s1_t) + (1 - p1_t) * dnorm(testexes, m2_t, s2_t)
densXc <- p1_c * dnorm(testexes, m1_c, s1_c) + (1 - p1_c) * dnorm(testexes, m2_c, s2_c)
p_tgivenx <- prob_t * densXt / (prob_t * densXt + (1-prob_t) * densXc)
par(mfrow=c(1,3))
plot(testexes, densXt, ylim=c(0, max(c(densXt, densXc))))
abline(v=m_t, lty=3, col='red')
plot(testexes, densXc, ylim=c(0, max(c(densXt, densXc))))
plot(testexes, p_tgivenx)


res_mat_em <- matrix(NA, nrow=n_sims, ncol=10)
res_mat_lin <- matrix(NA, nrow=n_sims, ncol=10)
bal_mats1 <- matrix(NA, nrow=n_sims, ncol=10)
bal_mats2 <- matrix(NA, nrow=n_sims, ncol=10)
opt_thetas <- matrix(NA, nrow=n_sims, ncol=2)
for(i in 1:n_sims){

  n_obs <- n_t + n_c
  Xt <- ifelse(rbinom(n_t, 1, p1_t), rnorm(n_t, m1_t, s1_t), rnorm(n_t, m2_t, s2_t))
  Xc <- ifelse(rbinom(n_c, 1, p1_c), rnorm(n_c, m1_c, s1_c), rnorm(n_c, m2_c, s2_c))
  X2 <- rbinom(n_obs, 1, prob = 0.4)

  X <- data.frame(X1 = c(Xt, Xc), X2 = X2)
  TA <- c(rep(1, n_t), rep(0, n_c))

  densXt <- p1_t * dnorm(X$X1, m1_t, s1_t) + (1 - p1_t) * dnorm(X$X1, m2_t, s2_t)
  densXc <- p1_c * dnorm(X$X1, m1_c, s1_c) + (1 - p1_c) * dnorm(X$X1, m2_c, s2_c)
  true_ps <- prob_t * densXt / (prob_t * densXt + (1-prob_t) * densXc)

  dset <- X
  dset$TA <- TA

  YT_em <- X[,1]^2 + 2 + rnorm(n_obs, sd=0.25)
  YC_em <- X[,1] + rnorm(n_obs, sd=0.25)
  YO_em <- ifelse(TA==1, YT_em, YC_em)

  YT_lin <- X[,1] + 3 + rnorm(n_obs, sd=0.25)
  YC_lin <- X[,1] + rnorm(n_obs, sd=0.25)
  YO_lin <- ifelse(TA==1, YT_lin, YC_lin)

  Xscaled <- cbind(scale(X[,1]), ifelse(X2==1, 1, -1))

  est_ps_glm1 <- fitted(glm('TA ~ poly(X1, 1) + X2', data=dset, family = 'binomial'))
  est_ps_glm2 <- fitted(glm('TA ~ poly(X1, 2) + X2', data=dset, family = 'binomial'))
  est_ps_gp1 <- paperOptBalGPPS::gpbal(Xscaled, TA,
                                 paperOptBalGPPS::sqexp_poly,
                                 # test_sqexp_poly,
                                 c(1,1),
                                 wts_vers = 'ATT')
  opt_thetas[i, ] <- est_ps_gp1$thetas

  est_ps_gp2 <- paperOptBalGPPS::gpbal(Xscaled, TA,
                                  paperOptBalGPPS::sqexp_ard,
                                  c(1, 0.5),
                                  wts_vers = 'ATT')
  bart_train <- BART::mc.pbart(Xscaled, TA, seed = i+2018)
  est_ps_bart <- pnorm(bart_train$yhat.train.mean)
  capture.output(est_ps_cbps <- as.vector(fitted(CBPS::CBPS('TA ~ poly(X1, 2) + X2',
                                                            data=dset, family = 'binomial'))))
  est_ps_gbm1 <- twang::ps(TA ~ X1 + X2, data=dset,
                           estimand = 'ATT', verbose=F,
                           stop.method = 'ks.mean')$ps[,1]
  est_ps_gbm2 <- twang::ps(TA ~ X1 + X2, data=dset,
                           estimand = 'ATT', verbose=F,
                           stop.method = 'es.mean')$ps[,1]
  est_ps_gbm3 <- twang::ps(TA ~ X1 + X2, data=dset,
                           estimand = 'ATT', verbose=F,
                           stop.method = 'es.max')$ps[,1]

  wts_true <- ifelse(TA==1, 1, true_ps/(1-true_ps))
  wts_glm1 <- ifelse(TA==1, 1, est_ps_glm1/(1-est_ps_glm1))
  wts_glm2 <- ifelse(TA==1, 1, est_ps_glm2/(1-est_ps_glm2))
  wts_gp1 <- ifelse(TA==1, 1, est_ps_gp1$ps/(1-est_ps_gp1$ps))
  wts_gp2 <- ifelse(TA==1, 1, est_ps_gp2$ps/(1-est_ps_gp2$ps))
  wts_bart <- ifelse(TA==1, 1, est_ps_bart/(1-est_ps_bart))
  wts_cbps <- ifelse(TA==1, 1, est_ps_cbps/(1-est_ps_cbps))
  wts_gbm1 <- ifelse(TA==1, 1, est_ps_gbm1/(1-est_ps_gbm1))
  wts_gbm2 <- ifelse(TA==1, 1, est_ps_gbm2/(1-est_ps_gbm2))
  wts_gbm3 <- ifelse(TA==1, 1, est_ps_gbm3/(1-est_ps_gbm3))

  bal_mats1[i, 1] <- paperOptBalGPPS::bal_table(X, 1:2, TA==1, wts=wts_true)[1, 7]
  bal_mats1[i, 2] <- paperOptBalGPPS::bal_table(X, 1:2, TA==1, wts=wts_gp1)[1, 7]
  bal_mats1[i, 3] <- paperOptBalGPPS::bal_table(X, 1:2, TA==1, wts=wts_gp2)[1, 7]
  bal_mats1[i, 4] <- paperOptBalGPPS::bal_table(X, 1:2, TA==1, wts=wts_bart)[1, 7]
  bal_mats1[i, 5] <- paperOptBalGPPS::bal_table(X, 1:2, TA==1, wts=wts_gbm1)[1, 7]
  bal_mats1[i, 6] <- paperOptBalGPPS::bal_table(X, 1:2, TA==1, wts=wts_gbm2)[1, 7]
  bal_mats1[i, 7] <- paperOptBalGPPS::bal_table(X, 1:2, TA==1, wts=wts_gbm3)[1, 7]
  bal_mats1[i, 8] <- paperOptBalGPPS::bal_table(X, 1:2, TA==1, wts=wts_cbps)[1, 7]
  bal_mats1[i, 9] <- paperOptBalGPPS::bal_table(X, 1:2, TA==1, wts=wts_glm2)[1, 7]
  bal_mats1[i, 10] <- paperOptBalGPPS::bal_table(X, 1:2, TA==1, wts=wts_glm1)[1, 7]

  bal_mats2[i, 1] <- paperOptBalGPPS::bal_table(X, 1:2, TA==1, wts=wts_true)[4, 7]
  bal_mats2[i, 2] <- paperOptBalGPPS::bal_table(X, 1:2, TA==1, wts=wts_gp1)[4, 7]
  bal_mats2[i, 3] <- paperOptBalGPPS::bal_table(X, 1:2, TA==1, wts=wts_gp2)[4, 7]
  bal_mats2[i, 4] <- paperOptBalGPPS::bal_table(X, 1:2, TA==1, wts=wts_bart)[4, 7]
  bal_mats2[i, 5] <- paperOptBalGPPS::bal_table(X, 1:2, TA==1, wts=wts_gbm1)[4, 7]
  bal_mats2[i, 6] <- paperOptBalGPPS::bal_table(X, 1:2, TA==1, wts=wts_gbm2)[4, 7]
  bal_mats2[i, 7] <- paperOptBalGPPS::bal_table(X, 1:2, TA==1, wts=wts_gbm3)[4, 7]
  bal_mats2[i, 8] <- paperOptBalGPPS::bal_table(X, 1:2, TA==1, wts=wts_cbps)[4, 7]
  bal_mats2[i, 9] <- paperOptBalGPPS::bal_table(X, 1:2, TA==1, wts=wts_glm2)[4, 7]
  bal_mats2[i, 10] <- paperOptBalGPPS::bal_table(X, 1:2, TA==1, wts=wts_glm1)[4, 7]

  res_mat_em[i, 1] <- lm_ps(YO_em, cbind(1, TA), wts = wts_true)$ests[2,1]
  res_mat_em[i, 2] <- lm_ps(YO_em, cbind(1, TA), wts = wts_gp1)$ests[2,1]
  res_mat_em[i, 3] <- lm_ps(YO_em, cbind(1, TA), wts = wts_gp2)$ests[2,1]
  res_mat_em[i, 4] <- lm_ps(YO_em, cbind(1, TA), wts = wts_bart)$ests[2,1]
  res_mat_em[i, 5] <- lm_ps(YO_em, cbind(1, TA), wts = wts_gbm1)$ests[2,1]
  res_mat_em[i, 6] <- lm_ps(YO_em, cbind(1, TA), wts = wts_gbm2)$ests[2,1]
  res_mat_em[i, 7] <- lm_ps(YO_em, cbind(1, TA), wts = wts_gbm3)$ests[2,1]
  res_mat_em[i, 8] <- lm_ps(YO_em, cbind(1, TA), wts = wts_cbps)$ests[2,1]
  res_mat_em[i, 9] <- lm_ps(YO_em, cbind(1, TA), wts = wts_glm2)$ests[2,1]
  res_mat_em[i, 10] <- lm_ps(YO_em, cbind(1, TA), wts = wts_glm1)$ests[2,1]

  res_mat_lin[i, 1] <- lm_ps(YO_lin, cbind(1, TA), wts = wts_true)$ests[2,1]
  res_mat_lin[i, 2] <- lm_ps(YO_lin, cbind(1, TA), wts = wts_gp1)$ests[2,1]
  res_mat_lin[i, 3] <- lm_ps(YO_lin, cbind(1, TA), wts = wts_gp2)$ests[2,1]
  res_mat_lin[i, 4] <- lm_ps(YO_lin, cbind(1, TA), wts = wts_bart)$ests[2,1]
  res_mat_lin[i, 5] <- lm_ps(YO_lin, cbind(1, TA), wts = wts_gbm1)$ests[2,1]
  res_mat_lin[i, 6] <- lm_ps(YO_lin, cbind(1, TA), wts = wts_gbm2)$ests[2,1]
  res_mat_lin[i, 7] <- lm_ps(YO_lin, cbind(1, TA), wts = wts_gbm3)$ests[2,1]
  res_mat_lin[i, 8] <- lm_ps(YO_lin, cbind(1, TA), wts = wts_cbps)$ests[2,1]
  res_mat_lin[i, 9] <- lm_ps(YO_lin, cbind(1, TA), wts = wts_glm2)$ests[2,1]
  res_mat_lin[i, 10] <- lm_ps(YO_lin, cbind(1, TA), wts = wts_glm1)$ests[2,1]

  message(paste(i,':', sep=''), appendLF = 'F')
}

results_data <- list('Cov1_Balance' = bal_mats1,
                     'Cov2_Balance' = bal_mats2,
                     'LinearResults' = res_mat_lin,
                     'EffModResults' = res_mat_em,
                     'NPSE_OptThetas' = opt_thetas)



bias_em_mat <- res_mat_em - true_att_em
bias_lin_mat <- res_mat_lin - true_att_lin

mean_bal1 <- apply(abs(bal_mats1) < 0.2, 2, mean, na.rm=T)
mean_bal2 <- apply(abs(bal_mats2) < 0.2, 2, mean, na.rm=T)
mean_balb <- apply(abs(bal_mats1) < 0.2 & abs(bal_mats2) < 0.2, 2, mean, na.rm=T)

bias_lin_mean <- apply(bias_lin_mat, 2, mean, na.rm=T)
bias_lin_abs <- apply(abs(bias_lin_mat), 2, mean, na.rm=T)
bias_lin_sd <- apply(bias_lin_mat, 2, sd, na.rm=T)
bias_lin_mse <- apply(bias_lin_mat^2, 2, mean, na.rm=T)

bias_em_mean <- apply(bias_em_mat, 2, mean, na.rm=T)
bias_em_abs <- apply(abs(bias_em_mat), 2, mean, na.rm=T)
bias_em_sd <- apply(bias_em_mat, 2, sd, na.rm=T)
bias_em_mse <- apply(bias_em_mat^2, 2, mean, na.rm=T)

outro_lin <- rbind(mean_bal1, mean_bal2, mean_balb,
                   bias_lin_mean, bias_lin_abs, bias_lin_sd, bias_lin_mse)

outro_em <- rbind(mean_bal1, mean_bal2, mean_balb,
                  bias_em_mean, bias_em_abs, bias_em_sd, bias_em_mse)
colnames(outro_lin) <- c('TRUTH', 'OBGPPS:NPSE', 'OBGPPS:SE', 'BART', 'GBM:KS.MEAN', 'GBM:ES.MEAN', 'GBM:ES.MAX',
                         'CBPS', 'GLM:POLY2', 'GLM:POLY1')
colnames(outro_em) <- c('TRUTH', 'OBGPPS', 'OBGPPS:SE', 'BART', 'GBM:KS.MEAN', 'GBM:ES.MEAN', 'GBM:ES.MAX',
                        'CBPS', 'GLM:POLY2', 'GLM:POLY1')
print(t(outro_lin))
print(t(outro_em))


# saveRDS(results_data, paste(Sys.Date(), '-attsim-mixtures-results.rds', sep=''))
