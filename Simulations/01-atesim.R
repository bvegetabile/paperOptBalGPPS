setwd('~/git/paperOptBalGPPS/Simulations/')

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

odd_ps_func <- function(x1, x2, inter,
                        beta1, beta2, beta3,
                        alph1, alph2){
  lin <- beta1 *x1 + beta2 * x2 + beta3 * x1 * x2 + inter
  ps <- alph1 * pnorm(lin) + alph2
  return(ps)
}

even_ps_func <- function(x1, x2, inter,
                         beta1, beta2, beta3,
                         alph1, alph2){
  lin <- inter + beta1 * x1 + beta2 * x1^2 + beta3 * x2
  ps <- alph1 * pnorm(lin) + alph2
}

make_wts <- function(ta, ps){
  wts <- data.frame(t=(ta/ps),
                    c=((1-ta)/(1-ps)))
  wts <- ifelse(ta==1, wts$t, wts$c)
  return(wts)
}

sim_settings <- list('nonparametric_odd' = list('type' = 'odd',
                                                'n_obs' = 500,
                                                'beta1' = 4,
                                                'beta2' = -0.5,
                                                'beta3' = -3,
                                                'inter' = 0.5,
                                                'alph1' = 0.7,
                                                'alph2' = 0.15,
                                                'simtitle'='Non-Parametric Function - Odd'),
                     'nonparametric_even' = list('type' = 'even',
                                                 'n_obs' = 500,
                                                 'beta1' = 3,
                                                 'beta2' = -4,
                                                 'beta3' = -2,
                                                 'inter' = 2.5,
                                                 'alph1' = 0.75,
                                                 'alph2' = 0.125,
                                                 'simtitle'='Non-Parametric Function - Even'))

true_ate <- 3

for(ss in 1:length(sim_settings)){
  # for(ss in 1:1){

  sim_vals <- sim_settings[[ss]]

  sim_name <- names(sim_settings)[ss]
  fn_type <- sim_vals$type
  beta1 <- sim_vals$beta1
  beta2 <- sim_vals$beta2
  beta3 <- sim_vals$beta3
  inter <- sim_vals$inter
  alph1 <- sim_vals$alph1
  alph2 <- sim_vals$alph2
  n_obs <- sim_vals$n_obs

  # Plotting Simulation Example
  X1_exp <- rnorm(n_obs)
  X2_exp <- rbinom(n_obs, 1, prob = 0.4)
  if(fn_type=='odd'){
    ps_exp <- odd_ps_func(X1_exp, X2_exp, inter, beta1, beta2, beta3, alph1, alph2)
  } else{
    ps_exp <- even_ps_func(X1_exp, X2_exp, inter, beta1, beta2, beta3, alph1, alph2)
  }

  # pdf(paste('./plot_', sim_name, '.pdf', sep=''), height=5, width = 10)
  plot(X1_exp[X2_exp ==1], ps_exp[X2_exp ==1],
       xlim = range(X1_exp), ylim=c(0,1),
       pch=19, col=rgb(0.75,0,0,0.5),
       xlab=expression(X[1]),
       ylab='Probability of Treatment',
       main=sim_vals$simtitle)
  points(X1_exp[X2_exp ==0], ps_exp[X2_exp ==0],
         pch=4, col=rgb(0,0,0.75,0.5), lwd=3)
  legend('bottomright', c(expression(X[2] == 0),
                          expression(X[2] == 1)),
         pch=c(4, 19), col=c(rgb(0,0,1,0.5),
                             rgb(1,0,0,0.5)),
         lty=c(0, NA), lwd=c(3,NA))
  # dev.off()

  res_mat_em <- matrix(NA, nrow=n_sims, ncol=12)
  res_mat_lin <- matrix(NA, nrow=n_sims, ncol=12)
  bal_mats1 <- matrix(NA, nrow=n_sims, ncol=12)
  bal_mats2 <- matrix(NA, nrow=n_sims, ncol=12)
  for(s in 1:n_sims){
    X1 <- rnorm(n_obs)
    X2 <- rbinom(n_obs, 1, prob = 0.4)
    Xscaled <- as.matrix(cbind(scale(X1, center = T, scale = T),
                               ifelse(X2==1, 1, -1)))

    if(fn_type=='odd'){
      true_ps <- odd_ps_func(X1, X2, inter, beta1, beta2, beta3, alph1, alph2)
    } else{
      true_ps <- even_ps_func(X1, X2, inter, beta1, beta2, beta3, alph1, alph2)
    }


    TA <- rbinom(n_obs, 1, true_ps)
    dset <- data.frame(X1 = X1, X2 = X2, TA=TA)
    Xdes <- cbind(1, TA)

    # Constant Treatment Effects - TE Related to X1
    Yt_lin <- X1 + 3 + rnorm(n_obs, sd = 0.25)
    Yc_lin <- X1 + rnorm(n_obs, sd = 0.25)
    Yo_lin <- TA * Yt_lin + (1-TA)*Yc_lin

    # Effect Modification
    Yt_em <- X1^2 + 2 + rnorm(n_obs, sd = 0.25)
    Yc_em <- X1 + rnorm(n_obs, sd = 0.25)
    Yo_em <- TA * Yt_em + (1-TA)*Yc_em

    # Naive Estimates
    naive_ate_lin <- lm(Yo_lin ~ TA)
    naive_ate_em <- lm(Yo_em ~ TA)
    naive_X1bal <- paperOptBalGPPS::bal_stats(X1, TA, 'continuous')[7]
    naive_X2bal <- paperOptBalGPPS::bal_stats(X2, TA, 'binary')[7]

    # Estimation
    if(fn_type=='odd'){
      est_ps_glm1 <- as.vector(fitted(glm(TA ~ X1 + X2 + X1*X2, family = 'binomial')))
    } else {
      est_ps_glm1 <- as.vector(fitted(glm(TA ~ X1 + I(X1^2) + X2, family = 'binomial')))
    }
    if(fn_type=='odd'){
      est_ps_glm2 <- as.vector(fitted(glm(TA ~ X1, family = 'binomial')))
    } else {
      est_ps_glm2 <- as.vector(fitted(glm(TA ~ X1, family = 'binomial')))
    }
    est_ps_gp1 <- paperOptBalGPPS::gpbal(Xscaled, TA,
                                    paperOptBalGPPS::sqexp_poly,
                                    c(1,1),
                                    wts_vers = 'ATE')
    est_ps_gp2 <- paperOptBalGPPS::gpbal(Xscaled, TA,
                                    paperOptBalGPPS::sqexp_ard,
                                    c(1, 0.5),
                                    wts_vers = 'ATE')
    bart_train <- BART::mc.pbart(Xscaled, TA, seed = s+2018)
    est_ps_bart <- pnorm(bart_train$yhat.train.mean)
    if(fn_type=='odd'){
      capture.output(est_ps_cbps1 <- as.vector(fitted(CBPS::CBPS(TA ~ X1 + X2 + X1*X2,
                                                                 ATT = 0, family = 'binomial'))))
    } else {
      capture.output(est_ps_cbps1 <- as.vector(fitted(CBPS::CBPS(TA ~ X1 + I(X1^2) + X2,
                                                                 ATT = 0, family = 'binomial'))))
    }
    if(fn_type=='odd'){
      capture.output(est_ps_cbps2 <- as.vector(fitted(CBPS::CBPS(TA ~ X1 + X2,
                                                                 ATT = 0, family = 'binomial'))))
    } else {
      capture.output(est_ps_cbps2 <- as.vector(fitted(CBPS::CBPS(TA ~ X1 + X2,
                                                                 ATT = 0, family = 'binomial'))))
    }
    est_ps_gbm1 <- twang::ps(TA ~ X1 + X2, data=dset,
                             estimand = 'ATE', verbose=F,
                             stop.method = 'ks.mean')$ps[,1]
    est_ps_gbm2 <- twang::ps(TA ~ X1 + X2, data=dset,
                             estimand = 'ATE', verbose=F,
                             stop.method = 'es.mean')$ps[,1]
    est_ps_gbm3 <- twang::ps(TA ~ X1 + X2, data=dset,
                             estimand = 'ATE', verbose=F,
                             stop.method = 'es.max')$ps[,1]

    wts_true <- make_wts(TA, true_ps)
    wts_glm1 <- make_wts(TA, est_ps_glm1)
    wts_glm2 <- make_wts(TA, est_ps_glm2)
    wts_gp1 <- make_wts(TA, est_ps_gp1$ps)
    wts_gp2 <- make_wts(TA, est_ps_gp2$ps)
    wts_bart <- make_wts(TA, est_ps_bart)
    wts_cbps1 <- make_wts(TA, est_ps_cbps1)
    wts_cbps2 <- make_wts(TA, est_ps_cbps2)
    wts_gbm1 <- make_wts(TA, est_ps_gbm1)
    wts_gbm2 <- make_wts(TA, est_ps_gbm2)
    wts_gbm3 <- make_wts(TA, est_ps_gbm3)

    bal_mats1[s, 1] <- naive_X1bal
    bal_mats1[s, 2] <- paperOptBalGPPS::bal_stats(X1, TA, 'continuous', wts_true)[7]
    bal_mats1[s, 3] <- paperOptBalGPPS::bal_stats(X1, TA, 'continuous', wts_gp1)[7]
    bal_mats1[s, 4] <- paperOptBalGPPS::bal_stats(X1, TA, 'continuous', wts_gp2)[7]
    bal_mats1[s, 5] <- paperOptBalGPPS::bal_stats(X1, TA, 'continuous', wts_bart)[7]
    bal_mats1[s, 6] <- paperOptBalGPPS::bal_stats(X1, TA, 'continuous', wts_gbm1)[7]
    bal_mats1[s, 7] <- paperOptBalGPPS::bal_stats(X1, TA, 'continuous', wts_gbm2)[7]
    bal_mats1[s, 8] <- paperOptBalGPPS::bal_stats(X1, TA, 'continuous', wts_gbm3)[7]
    bal_mats1[s, 9] <- paperOptBalGPPS::bal_stats(X1, TA, 'continuous', wts_glm1)[7]
    bal_mats1[s, 10] <- paperOptBalGPPS::bal_stats(X1, TA, 'continuous', wts_cbps1)[7]
    bal_mats1[s, 11] <- paperOptBalGPPS::bal_stats(X1, TA, 'continuous', wts_glm2)[7]
    bal_mats1[s, 12] <- paperOptBalGPPS::bal_stats(X1, TA, 'continuous', wts_cbps2)[7]


    bal_mats2[s, 1] <- naive_X2bal
    bal_mats2[s, 2] <- paperOptBalGPPS::bal_stats(X2, TA, 'binary', wts_true)[7]
    bal_mats2[s, 3] <- paperOptBalGPPS::bal_stats(X2, TA, 'binary', wts_gp1)[7]
    bal_mats2[s, 4] <- paperOptBalGPPS::bal_stats(X2, TA, 'binary', wts_gp2)[7]
    bal_mats2[s, 5] <- paperOptBalGPPS::bal_stats(X2, TA, 'binary', wts_bart)[7]
    bal_mats2[s, 6] <- paperOptBalGPPS::bal_stats(X2, TA, 'binary', wts_gbm1)[7]
    bal_mats2[s, 7] <- paperOptBalGPPS::bal_stats(X2, TA, 'binary', wts_gbm2)[7]
    bal_mats2[s, 8] <- paperOptBalGPPS::bal_stats(X2, TA, 'binary', wts_gbm3)[7]
    bal_mats2[s, 9] <- paperOptBalGPPS::bal_stats(X2, TA, 'binary', wts_glm1)[7]
    bal_mats2[s, 10] <- paperOptBalGPPS::bal_stats(X2, TA, 'binary', wts_cbps1)[7]
    bal_mats2[s, 11] <- paperOptBalGPPS::bal_stats(X2, TA, 'binary', wts_glm2)[7]
    bal_mats2[s, 12] <- paperOptBalGPPS::bal_stats(X2, TA, 'binary', wts_cbps2)[7]

    res_mat_lin[s, 1] <- naive_ate_lin$coefficients[2]
    res_mat_lin[s, 2] <- lm_ps(Yo_lin, cbind(1, TA), wts = wts_true)$ests[2,1]
    res_mat_lin[s, 3] <- lm_ps(Yo_lin, cbind(1, TA), wts = wts_gp1)$ests[2,1]
    res_mat_lin[s, 4] <- lm_ps(Yo_lin, cbind(1, TA), wts = wts_gp2)$ests[2,1]
    res_mat_lin[s, 5] <- lm_ps(Yo_lin, cbind(1, TA), wts = wts_bart)$ests[2,1]
    res_mat_lin[s, 6] <- lm_ps(Yo_lin, cbind(1, TA), wts = wts_gbm1)$ests[2,1]
    res_mat_lin[s, 7] <- lm_ps(Yo_lin, cbind(1, TA), wts = wts_gbm2)$ests[2,1]
    res_mat_lin[s, 8] <- lm_ps(Yo_lin, cbind(1, TA), wts = wts_gbm3)$ests[2,1]
    res_mat_lin[s, 9] <- lm_ps(Yo_lin, cbind(1, TA), wts = wts_glm1)$ests[2,1]
    res_mat_lin[s, 10] <- lm_ps(Yo_lin, cbind(1, TA), wts = wts_cbps1)$ests[2,1]
    res_mat_lin[s, 11] <- lm_ps(Yo_lin, cbind(1, TA), wts = wts_glm2)$ests[2,1]
    res_mat_lin[s, 12] <- lm_ps(Yo_lin, cbind(1, TA), wts = wts_cbps2)$ests[2,1]

    res_mat_em[s, 1] <- naive_ate_em$coefficients[2]
    res_mat_em[s, 2] <- lm_ps(Yo_em, cbind(1, TA), wts = wts_true)$ests[2,1]
    res_mat_em[s, 3] <- lm_ps(Yo_em, cbind(1, TA), wts = wts_gp1)$ests[2,1]
    res_mat_em[s, 4] <- lm_ps(Yo_em, cbind(1, TA), wts = wts_gp2)$ests[2,1]
    res_mat_em[s, 5] <- lm_ps(Yo_em, cbind(1, TA), wts = wts_bart)$ests[2,1]
    res_mat_em[s, 6] <- lm_ps(Yo_em, cbind(1, TA), wts = wts_gbm1)$ests[2,1]
    res_mat_em[s, 7] <- lm_ps(Yo_em, cbind(1, TA), wts = wts_gbm2)$ests[2,1]
    res_mat_em[s, 8] <- lm_ps(Yo_em, cbind(1, TA), wts = wts_gbm3)$ests[2,1]
    res_mat_em[s, 9] <- lm_ps(Yo_em, cbind(1, TA), wts = wts_glm1)$ests[2,1]
    res_mat_em[s, 10] <- lm_ps(Yo_em, cbind(1, TA), wts = wts_cbps1)$ests[2,1]
    res_mat_em[s, 11] <- lm_ps(Yo_em, cbind(1, TA), wts = wts_glm2)$ests[2,1]
    res_mat_em[s, 12] <- lm_ps(Yo_em, cbind(1, TA), wts = wts_cbps2)$ests[2,1]

    message(paste(s,':', sep=''), appendLF = 'F')
  }

  results_data <- list('Cov1_Balance' = bal_mats1,
                       'Cov2_Balance' = bal_mats2,
                       'LinearResults' = res_mat_lin,
                       'EffModResults' = res_mat_em)



  bias_em_mat <- res_mat_em - true_ate
  bias_lin_mat <- res_mat_lin - true_ate

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

  colnames(outro_lin) <- c('NAIVE', 'TRUEPS', 'OBGPPS:NPSE', 'OBGPPS:SE',
                           'BART', 'GBM:KS.MEAN', 'GBM:ES.MEAN', 'GBM:ES.MAX',
                           'GLM:CORRECT', 'CBPS:CORRECT', 'GLM:MISSPECIFIED', 'CBPS:MISSPECIFIED')
  colnames(outro_em) <- c('NAIVE', 'TRUEPS', 'OBGPPS', 'OBGPPS:SE',
                          'BART', 'GBM:KS.MEAN', 'GBM:ES.MEAN', 'GBM:ES.MAX',
                          'GLM:CORRECT', 'CBPS:CORRECT', 'GLM:MISSPECIFIED', 'CBPS:MISSPECIFIED')
  print(t(outro_lin))
  print(t(outro_em))

  # saveRDS(results_data, paste(Sys.Date(), '-', sim_name, '-atesim-results.rds', sep=''))
}
