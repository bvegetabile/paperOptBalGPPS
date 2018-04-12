setwd('~/git/paperOptBalGPPS/')

#-------------------------------------------------------------------------------
# Number of sims to run
set.seed(89274)
n_sims <- 5
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


dset_sizes <- c(100, 500, 1000)
time_results <- list()

for(ss in 1:1){
  # for(ss in 1:1){

  for(ds in 1:length(dset_sizes)){

    time_mat <- matrix(NA, nrow=n_sims, ncol=6)

    sim_vals <- sim_settings[[ss]]

    sim_name <- names(sim_settings)[ss]
    fn_type <- sim_vals$type
    beta1 <- sim_vals$beta1
    beta2 <- sim_vals$beta2
    beta3 <- sim_vals$beta3
    inter <- sim_vals$inter
    alph1 <- sim_vals$alph1
    alph2 <- sim_vals$alph2
    n_obs <- dset_sizes[ds]

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
      glm_start <- Sys.time()
      if(fn_type=='odd'){
        est_ps_glm1 <- as.vector(fitted(glm(TA ~ X1 + X2 + X1*X2, family = 'binomial')))
      } else {
        est_ps_glm1 <- as.vector(fitted(glm(TA ~ X1 + I(X1^2) + X2, family = 'binomial')))
      }
      glm_end <- Sys.time()
      glm_diff <- difftime(glm_end, glm_start, units = 'secs')

      gp1_start <- Sys.time()
      est_ps_gp1 <- paperOptBalGPPS::gpbal(Xscaled, TA,
                                      paperOptBalGPPS::sqexp_poly,
                                      c(1,1),
                                      wts_vers = 'ATE')
      gp1_end <- Sys.time()
      gp1_diff <- difftime(gp1_end, gp1_start, units = 'secs')

      gp2_start <- Sys.time()
      est_ps_gp2 <- paperOptBalGPPS::gpbal(Xscaled, TA,
                                      paperOptBalGPPS::sqexp_ard,
                                      c(1, 0.5),
                                      wts_vers = 'ATE')
      gp2_end <- Sys.time()
      gp2_diff <- difftime(gp2_end, gp2_start, units = 'secs')

      bart_start <- Sys.time()
      bart_train <- BART::mc.pbart(Xscaled, TA, seed = s+2018)
      est_ps_bart <- pnorm(bart_train$yhat.train.mean)
      bart_end <- Sys.time()
      bart_diff <- difftime(bart_end, bart_start, units = 'secs')


      if(fn_type=='odd'){
        cbps_start <- Sys.time()
        capture.output(est_ps_cbps1 <- as.vector(fitted(CBPS::CBPS(TA ~ X1 + X2 + X1*X2,
                                                                   ATT = 0, family = 'binomial'))))
      } else {
        cbps_start <- Sys.time()
        capture.output(est_ps_cbps1 <- as.vector(fitted(CBPS::CBPS(TA ~ X1 + I(X1^2) + X2,
                                                                   ATT = 0, family = 'binomial'))))
      }
      cbps_end <- Sys.time()
      cbps_diff <- difftime(cbps_end, cbps_start, units = 'secs')

      gbm_start <- Sys.time()
      est_ps_gbm3 <- twang::ps(TA ~ X1 + X2, data=dset,
                               estimand = 'ATE', verbose=F,
                               stop.method = 'es.max')$ps[,1]
      gbm_end <- Sys.time()
      gbm_diff <- difftime(gbm_end, gbm_start, units = 'secs')

      time_mat[s, 1] <- glm_diff[[1]]
      time_mat[s, 2] <- gp1_diff[[1]]
      time_mat[s, 3] <- gp2_diff[[1]]
      time_mat[s, 4] <- bart_diff[[1]]
      time_mat[s, 5] <- cbps_diff[[1]]
      time_mat[s, 6] <- gbm_diff[[1]]

      message(paste(s,':', sep=''), appendLF = 'F')
  }
  time_results[[ds]] <- time_mat
  }
}

outro <- matrix(NA, nrow=length(dset_sizes), ncol=6)
for(i in 1:length(dset_sizes)){
  outro[i,] <- apply(time_results[[i]], 2, mean)
}
colnames(outro) <- c('GLM', 'OBGPPS:NPSE', 'OBGPPS:SE', 'BART', 'CBPS', 'GBM')
rownames(outro) <- dset_sizes
print(t(outro))
