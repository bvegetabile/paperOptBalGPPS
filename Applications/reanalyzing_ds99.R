library(paperOptBalGPPS)
library(xtable)
setwd('~/Documents/GitHub/paperOptBalGPPS/Applications/')

# Total Covariate Balance ------------
mom_sq_bal <- function(dataset,
                       col_ind,
                       treat_ind,
                       wts = rep(1, length(treat_ind)),
                       max_uniq=5){

  treat_ind <- as.logical(treat_ind)
  outtable <- c()
  outbal <- 0.0
  counter <- 1
  # Iteration based on the order of column numbers provided
  for(i in 1:length(col_ind)){
    c <- col_ind[i]
    col_data <- dataset[, c]
    # Conversion of variables to categorical if it has less than the max_uniq
    # categories.
    if(length(unique(col_data)) <= max_uniq){
      col_data <- as.factor(col_data)
    }
    if(is.factor(col_data)){
      obs_lvls <- levels(col_data)
      outtable <- c()
      for(lvl in obs_lvls){
        stddiff <- bal_stats(col_data==lvl, treat_ind, 'binary', wts)
        outtable <- rbind(outtable, stddiff)
      }
      outbal <- outbal + sum(abs(outtable[2:nrow(outtable),7])^2) + sum(abs(outtable[2:nrow(outtable),8])^2)
    } else {
      stddiff <- bal_stats(col_data, treat_ind, 'continuous', wts)
      outbal <- outbal + abs(stddiff[7]^2) + abs(stddiff[8]^2)
    }
  }
  return(outbal)
}

# Merging Experimental and Observational Data -----

make_dataset <- function(treated, control){
  treated$ta <- 1
  control$ta <- 0
  rbind(treated,control)
}

# Creating a Scaled Design Matrix For Propensity Score Estimation -------------

make_Xmat <- function(X){
  for(c in 1:ncol(X)){
    if(length(unique(X[,c])) > 2){
      X[,c] <- scale(X[,c])
    } else{
      X[,c] <- ifelse(X[,c] == 1, 1, -1)
    }
  }
  X
}

# Creating ATT Weights --------

make_wts2 <- function(ta, ps){
  wts <- ta + (1-ta)*ps / (1-ps)
  as.vector(wts)
}

# Weighted Least Squares Function with Robust Variance --------

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

# Creating Design Matrix for ATT Estimation cbind(1 , TA) -----------

make_design <- function(dset){
  Xdes <- cbind(1, dset$ta)
  Xdes
}

#

# Function for statistics -------------------

mean_var <- function(dset, col){
  mean_val <- round(dset[,col],2)
  sd_val <- paste('(',round(sqrt(dset[,col+1]),2),')',sep='')
  cols <- cbind(mean_val, sd_val)
  stat_vals <- apply(cols, 1, paste, collapse=" ")
  stat_vals[c(1,2,5,8,11,14,15,16)]
}

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#
# Data Analysis - Some notes
# --- It is not advised to run ALL of the analyses -> CPS-1 in particular takes
#     a long time to finish running
# --- Github does not like how large the files are that the function outputs
#     --> This is why they are not included
#     --> A Dropbox with these files may be provided
#
#-------------------------------------------------------------------------------

# Loading all data ------------

lalonde_header <- scan('./data/nswre74_headers.txt', what=" ")
lalonde_treated <- read.fwf('./data/nswre74_treated.txt', widths=rep(16,10), col.names = lalonde_header)
lalonde_control <- read.fwf('./data/nswre74_control.txt', widths=rep(16,10), col.names = lalonde_header)

psid1_control <- read.fwf('./data/psid_controls.txt', widths=rep(16,10), col.names = lalonde_header)
psid2_control <- read.fwf('./data/psid2_controls.txt', widths=rep(16,10), col.names = lalonde_header)
psid3_control <- read.fwf('./data/psid3_controls.txt', widths=rep(16,10), col.names = lalonde_header)

cps_header <- scan('./data/cps_headers.txt', what=" ")
cps1_control <- read.fwf('./data/cps_controls.txt', widths=rep(16,10), col.names = cps_header)
cps2_control <- read.fwf('./data/cps2_controls.txt', widths=rep(16,10), col.names = cps_header)
cps3_control <- read.fwf('./data/cps3_controls.txt', widths=rep(16,10), col.names = cps_header)

# Making datasets
data_obs_psid1 <- make_dataset(lalonde_treated, psid1_control)
data_obs_psid2 <- make_dataset(lalonde_treated, psid2_control)
data_obs_psid3 <- make_dataset(lalonde_treated, psid3_control)
data_obs_cps1 <- make_dataset(lalonde_treated, cps1_control)
data_obs_cps2 <- make_dataset(lalonde_treated, cps2_control)
data_obs_cps3 <- make_dataset(lalonde_treated, cps3_control)


data_obs_psid1$u74 <- ifelse(data_obs_psid1$RE74 == 0, 1, 0)
data_obs_psid1$u75 <- ifelse(data_obs_psid1$RE75 == 0, 1, 0)
data_obs_psid2$u74 <- ifelse(data_obs_psid2$RE74 == 0, 1, 0)
data_obs_psid2$u75 <- ifelse(data_obs_psid2$RE75 == 0, 1, 0)
data_obs_psid3$u74 <- ifelse(data_obs_psid3$RE74 == 0, 1, 0)
data_obs_psid3$u75 <- ifelse(data_obs_psid3$RE75 == 0, 1, 0)

data_obs_cps1$u74 <- ifelse(data_obs_cps1$RE74 == 0, 1, 0)
data_obs_cps1$u75 <- ifelse(data_obs_cps1$RE75 == 0, 1, 0)
data_obs_cps2$u74 <- ifelse(data_obs_cps2$RE74 == 0, 1, 0)
data_obs_cps2$u75 <- ifelse(data_obs_cps2$RE75 == 0, 1, 0)
data_obs_cps3$u74 <- ifelse(data_obs_cps3$RE74 == 0, 1, 0)
data_obs_cps3$u75 <- ifelse(data_obs_cps3$RE75 == 0, 1, 0)


# PSID-1 ANALYSIS
dsw99_psid1mod <- glm(ta ~ age + I(age^2) + education + I(education^2) + married + nodegree + black + hispanic + RE74 + RE75 + I(RE74^2) + I(RE75^2) + I(u74 * black),
                      data = data_obs_psid1, family = 'binomial')
data_obs_psid1$ps <- predict(dsw99_psid1mod, type = 'response')
data_obs_psid1$wts <- make_wts2(data_obs_psid1$ta, data_obs_psid1$ps)
psid1_X <- cbind(1, data_obs_psid1$ta)
min_t_ps <- min(data_obs_psid1$ps[data_obs_psid1$ta == 1])
psid1_mask <- data_obs_psid1$ps > min_t_ps
psid1_est1 <- lm_ps(data_obs_psid1$RE78, psid1_X, data_obs_psid1$wts, true_val = 0)
psid1_est2 <- lm_ps(data_obs_psid1$RE78[psid1_mask],
              psid1_X[psid1_mask,],
              data_obs_psid1$wts[psid1_mask], true_val = 0)
paperOptBalGPPS::bal_table(data_obs_psid1, 2:9, as.logical(data_obs_psid1$ta), wts = data_obs_psid1$wts)
bal_total_psid1 <- mom_sq_bal(data_obs_psid1, 2:9, as.logical(data_obs_psid1$ta), wts =  data_obs_psid1$wts)

# PSID-2 ANALYSIS
dsw99_psid2mod <- glm(ta ~ age + I(age^2) + education + I(education^2) + married + nodegree + black + hispanic + RE74 + RE75 + I(RE74^2) + I(RE75^2) + u74 + u75,
                      data = data_obs_psid2, family = 'binomial')
data_obs_psid2$ps <- predict(dsw99_psid2mod, type = 'response')
data_obs_psid2$wts <- make_wts2(data_obs_psid2$ta, data_obs_psid2$ps)
psid2_X <- cbind(1, data_obs_psid2$ta)
min_t_ps <- min(data_obs_psid2$ps[data_obs_psid2$ta == 1])
psid2_mask <- data_obs_psid2$ps > min_t_ps
psid2_est1 <- lm_ps(data_obs_psid2$RE78, psid2_X, data_obs_psid2$wts, true_val = 0)
psid2_est2 <- lm_ps(data_obs_psid2$RE78[psid2_mask],
                    psid2_X[psid2_mask,],
                    data_obs_psid2$wts[psid2_mask], true_val = 0)
paperOptBalGPPS::bal_table(data_obs_psid2, 2:9, as.logical(data_obs_psid2$ta), wts = data_obs_psid2$wts)
bal_total_psid2 <- mom_sq_bal(data_obs_psid2, 2:9, as.logical(data_obs_psid2$ta), wts =  data_obs_psid2$wts)

# PSID-2 ANALYSIS
dsw99_psid3mod <- glm(ta ~ age + I(age^2) + education + I(education^2) + married + nodegree + black + hispanic + RE74 + RE75 + I(RE74^2) + I(RE75^2) + u74 + u75,
                      data = data_obs_psid3, family = 'binomial')
data_obs_psid3$ps <- predict(dsw99_psid3mod, type = 'response')
data_obs_psid3$wts <- make_wts2(data_obs_psid3$ta, data_obs_psid3$ps)
psid3_X <- cbind(1, data_obs_psid3$ta)
min_t_ps <- min(data_obs_psid3$ps[data_obs_psid3$ta == 1])
psid3_mask <- data_obs_psid3$ps > min_t_ps
psid3_est1 <- lm_ps(data_obs_psid3$RE78, psid3_X, data_obs_psid3$wts, true_val = 0)
psid3_est2 <- lm_ps(data_obs_psid3$RE78[psid3_mask],
                    psid3_X[psid3_mask,],
                    data_obs_psid3$wts[psid3_mask], true_val = 0)
paperOptBalGPPS::bal_table(data_obs_psid3, 2:9, as.logical(data_obs_psid3$ta), wts = data_obs_psid3$wts)
bal_total_psid3 <- mom_sq_bal(data_obs_psid3, 2:9, as.logical(data_obs_psid3$ta), wts =  data_obs_psid3$wts)
print(rbind(c(bal_total_psid1, psid1_est1$ests[2,1],psid1_est1$ests[2,2]),
            c(bal_total_psid2, psid2_est1$ests[2,1],psid2_est1$ests[2,2]),
            c(bal_total_psid3, psid3_est1$ests[2,1],psid3_est1$ests[2,2])))


# cps-1 ANALYSIS
dsw99_cps1mod <- glm(ta ~ age + I(age^2) + I(age^3) + education + I(education^2) + married + nodegree + black + hispanic + RE74 + RE75 + u74 + u75 + I(education * black),
                      data = data_obs_cps1, family = 'binomial')
data_obs_cps1$ps <- predict(dsw99_cps1mod, type = 'response')
data_obs_cps1$wts <- make_wts2(data_obs_cps1$ta, data_obs_cps1$ps)
cps1_X <- cbind(1, data_obs_cps1$ta)
min_t_ps <- min(data_obs_cps1$ps[data_obs_cps1$ta == 1])
cps1_mask <- data_obs_cps1$ps > min_t_ps
cps1_est1 <- lm_ps(data_obs_cps1$RE78, cps1_X, data_obs_cps1$wts, true_val = 0)
cps1_est2 <- lm_ps(data_obs_cps1$RE78[cps1_mask],
                    cps1_X[cps1_mask,],
                    data_obs_cps1$wts[cps1_mask], true_val = 0)
paperOptBalGPPS::bal_table(data_obs_cps1, 2:9, as.logical(data_obs_cps1$ta), wts = data_obs_cps1$wts)
bal_total_cps1 <- mom_sq_bal(data_obs_cps1, 2:9, as.logical(data_obs_cps1$ta), wts =  data_obs_cps1$wts)

# cps-2 ANALYSIS
dsw99_cps2mod <- glm(ta ~ age + I(age^2) + I(age^3) + education + I(education^2) + married + nodegree + black + hispanic + RE74 + RE75 + u74 + u75 + I(education * black),
                      data = data_obs_cps2, family = 'binomial')
data_obs_cps2$ps <- predict(dsw99_cps2mod, type = 'response')
data_obs_cps2$wts <- make_wts2(data_obs_cps2$ta, data_obs_cps2$ps)
cps2_X <- cbind(1, data_obs_cps2$ta)
min_t_ps <- min(data_obs_cps2$ps[data_obs_cps2$ta == 1])
cps2_mask <- data_obs_cps2$ps > min_t_ps
cps2_est1 <- lm_ps(data_obs_cps2$RE78, cps2_X, data_obs_cps2$wts, true_val = 0)
cps2_est2 <- lm_ps(data_obs_cps2$RE78[cps2_mask],
                    cps2_X[cps2_mask,],
                    data_obs_cps2$wts[cps2_mask], true_val = 0)
paperOptBalGPPS::bal_table(data_obs_cps2, 2:9, as.logical(data_obs_cps2$ta), wts = data_obs_cps2$wts)
bal_total_cps2 <- mom_sq_bal(data_obs_cps2, 2:9, as.logical(data_obs_cps2$ta), wts =  data_obs_cps2$wts)

# cps-2 ANALYSIS
dsw99_cps3mod <- glm(ta ~ age + I(age^2) + I(age^3) + education + I(education^2) + married + nodegree + black + hispanic + RE74 + RE75 + u74 + u75 + I(education * black),
                      data = data_obs_cps3, family = 'binomial')
data_obs_cps3$ps <- predict(dsw99_cps3mod, type = 'response')
data_obs_cps3$wts <- make_wts2(data_obs_cps3$ta, data_obs_cps3$ps)
cps3_X <- cbind(1, data_obs_cps3$ta)
min_t_ps <- min(data_obs_cps3$ps[data_obs_cps3$ta == 1])
cps3_mask <- data_obs_cps3$ps > min_t_ps
cps3_est1 <- lm_ps(data_obs_cps3$RE78, cps3_X, data_obs_cps3$wts, true_val = 0)
cps3_est2 <- lm_ps(data_obs_cps3$RE78[cps3_mask],
                    cps3_X[cps3_mask,],
                    data_obs_cps3$wts[cps3_mask], true_val = 0)
paperOptBalGPPS::bal_table(data_obs_cps3, 2:9, as.logical(data_obs_cps3$ta), wts = data_obs_cps3$wts)
bal_total_cps3 <- mom_sq_bal(data_obs_cps3, 2:9, as.logical(data_obs_cps3$ta), wts =  data_obs_cps3$wts)
print(rbind(c(bal_total_cps1, cps1_est1$ests[2,1],cps1_est1$ests[2,2]),
            c(bal_total_cps2, cps2_est1$ests[2,1],cps2_est1$ests[2,2]),
            c(bal_total_cps3, cps3_est1$ests[2,1],cps3_est1$ests[2,2])))

mom_sq_bal(data_obs_cps3,
           2:9, as.logical(data_obs_cps3$ta),
           wts =  ifelse(data_obs_cps3$ta == 1, 1/ sum(data_obs_cps3$ta == 1), 1/ sum(data_obs_cps3$ta == 0)))
