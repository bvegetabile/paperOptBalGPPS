library(gpbalancer)
library(xtable)
setwd('~/Dropbox/ucirvine/research/lalondeData/')

mom_bal <- function(dataset,
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
            outbal <- outbal + sum(abs(outtable[2:nrow(outtable),7])) + sum(abs(outtable[2:nrow(outtable),8]))
        } else {
            stddiff <- bal_stats(col_data, treat_ind, 'continuous', wts)
            outbal <- outbal + abs(stddiff[7]) + abs(stddiff[8])
        }
    }
    return(outbal)
}

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


lalonde_header <- scan('./nswre74_headers.txt', what=" ")
lalonde_treated <- read.fwf('./nswre74_treated.txt', widths=rep(16,10), col.names = lalonde_header)
lalonde_control <- read.fwf('./nswre74_control.txt', widths=rep(16,10), col.names = lalonde_header)

psid1_control <- read.fwf('./psid_controls.txt', widths=rep(16,10), col.names = lalonde_header)
psid2_control <- read.fwf('./psid2_controls.txt', widths=rep(16,10), col.names = lalonde_header)
psid3_control <- read.fwf('./psid3_controls.txt', widths=rep(16,10), col.names = lalonde_header)

cps_header <- scan('./cps_headers.txt', what=" ")
cps1_control <- read.fwf('./cps_controls.txt', widths=rep(16,10), col.names = cps_header)
cps2_control <- read.fwf('./cps2_controls.txt', widths=rep(16,10), col.names = cps_header)
cps3_control <- read.fwf('./cps3_controls.txt', widths=rep(16,10), col.names = cps_header)

setwd('~/Dropbox/ucirvine/research/papers/2017_obgpps/2018-03-DW99Replication/')


make_dataset <- function(treated, control){
    treated$ta <- 1
    control$ta <- 0
    rbind(treated,control)
}

data_experimental <- make_dataset(lalonde_treated, lalonde_control)
data_obs_psid1 <- make_dataset(lalonde_treated, psid1_control)
data_obs_psid2 <- make_dataset(lalonde_treated, psid2_control)
data_obs_psid3 <- make_dataset(lalonde_treated, psid3_control)
data_obs_cps1 <- make_dataset(lalonde_treated, cps1_control)
data_obs_cps2 <- make_dataset(lalonde_treated, cps2_control)
data_obs_cps3 <- make_dataset(lalonde_treated, cps3_control)

bal_before_exp <- gpbalancer::bal_table(data_experimental, 2:9, as.logical(data_experimental$ta))
bal_before_psid1 <- gpbalancer::bal_table(data_obs_psid1, 2:9, as.logical(data_obs_psid1$ta))
bal_before_psid2 <- gpbalancer::bal_table(data_obs_psid2, 2:9, as.logical(data_obs_psid2$ta))
bal_before_psid3 <- gpbalancer::bal_table(data_obs_psid3, 2:9, as.logical(data_obs_psid3$ta))
bal_before_cps1 <- gpbalancer::bal_table(data_obs_cps1, 2:9, as.logical(data_obs_cps1$ta))
bal_before_cps2 <- gpbalancer::bal_table(data_obs_cps2, 2:9, as.logical(data_obs_cps2$ta))
bal_before_cps3 <- gpbalancer::bal_table(data_obs_cps3, 2:9, as.logical(data_obs_cps3$ta))

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

Xmat_psid1 <- make_Xmat(data_obs_psid1[,2:9])
Xmat_psid2 <- make_Xmat(data_obs_psid2[,2:9])
Xmat_psid3 <- make_Xmat(data_obs_psid3[,2:9])
Xmat_cps1 <- make_Xmat(data_obs_cps1[,2:9])
Xmat_cps2 <- make_Xmat(data_obs_cps2[,2:9])
Xmat_cps3 <- make_Xmat(data_obs_cps3[,2:9])

# inits <- c(1,1,0.5,0.5,0.5,0.5,1,1)
#
est_ps_psid1 <- gpbalancer::gpbal(Xmat_psid1,
                                  data_obs_psid1$ta,
                                  gpbalancer::sqexp_poly,
                                  c(1,1),
                                  verbose=T,
                                  wts_vers = 'att')
est_ps_psid2 <- gpbalancer::gpbal(Xmat_psid2,
                                  data_obs_psid2$ta,
                                  gpbalancer::sqexp_poly,
                                  c(1,1),
                                  verbose=T,
                                  wts_vers = 'att')
est_ps_psid3 <- gpbalancer::gpbal(Xmat_psid3,
                                  data_obs_psid3$ta,
                                  gpbalancer::sqexp_poly,
                                  c(1,1),
                                  verbose=T,
                                  wts_vers = 'att')
est_ps_cps1 <- gpbalancer::gpbal(Xmat_cps1,
                                 data_obs_cps1$ta,
                                 gpbalancer::sqexp_poly,
                                 c(1,1),
                                 verbose=T,
                                 wts_vers = 'att')
est_ps_cps2 <- gpbalancer::gpbal(Xmat_cps2,
                                 data_obs_cps2$ta,
                                 gpbalancer::sqexp_poly,
                                 c(1,1),
                                 verbose=T,
                                 wts_vers = 'att')
est_ps_cps3 <- gpbalancer::gpbal(Xmat_cps3,
                                 data_obs_cps3$ta,
                                 gpbalancer::sqexp_poly,
                                 c(1,1),
                                 verbose=T,
                                 wts_vers = 'att')
#

saveRDS(est_ps_psid1, file = 'est_ps_psid1_sepoly.RData')
saveRDS(est_ps_psid2, file = 'est_ps_psid2_sepoly.RData')
saveRDS(est_ps_psid3, file = 'est_ps_psid3_sepoly.RData')
saveRDS(est_ps_cps1, file = 'est_ps_cps1_sepoly.RData')
saveRDS(est_ps_cps2, file = 'est_ps_cps2_sepoly.RData')
saveRDS(est_ps_cps3, file = 'est_ps_cps3_sepoly.RData')

est_ps_psid1 <- readRDS('est_ps_psid1_sepoly.RData')
est_ps_psid2 <- readRDS('est_ps_psid2_sepoly.RData')
est_ps_psid3 <- readRDS('est_ps_psid3_sepoly.RData')
est_ps_cps1 <- readRDS('est_ps_cps1_sepoly.RData')
est_ps_cps2 <- readRDS('est_ps_cps2_sepoly.RData')
est_ps_cps3 <- readRDS('est_ps_cps3_sepoly.RData')


make_wts2 <- function(ta, ps){
    wts <- ta + (1-ta)*ps / (1-ps)
    as.vector(wts)
}

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

make_design <- function(dset){
    Xdes <- cbind(1, dset$ta)
    Xdes
}

wts_psid1 <- make_wts2(data_obs_psid1$ta, est_ps_psid1$ps)
wts_psid2 <- make_wts2(data_obs_psid2$ta, est_ps_psid2$ps)
wts_psid3 <- make_wts2(data_obs_psid3$ta, est_ps_psid3$ps)
wts_cps1 <- make_wts2(data_obs_cps1$ta, est_ps_cps1$ps)
wts_cps2 <- make_wts2(data_obs_cps2$ta, est_ps_cps2$ps)
wts_cps3 <- make_wts2(data_obs_cps3$ta, est_ps_cps3$ps)

est_experiment <- lm_ps(data_experimental$RE78, X = make_design(data_experimental), wts = rep(1, nrow(data_experimental)))
est_ate_psid1 <- lm_ps(data_obs_psid1$RE78, X = make_design(data_obs_psid1), wts = wts_psid1)
est_ate_psid2 <- lm_ps(data_obs_psid2$RE78, X = make_design(data_obs_psid2), wts = wts_psid2)
est_ate_psid3 <- lm_ps(data_obs_psid3$RE78, X = make_design(data_obs_psid3), wts = wts_psid3)
est_ate_cps1 <- lm_ps(data_obs_cps1$RE78, make_design(data_obs_cps1), wts_cps1)
est_ate_cps2 <- lm_ps(data_obs_cps2$RE78, make_design(data_obs_cps2), wts_cps2)
est_ate_cps3 <- lm_ps(data_obs_cps3$RE78, make_design(data_obs_cps3), wts_cps3)

bal_experimental <- gpbalancer::bal_table(data_experimental, 2:9, as.logical(data_experimental$ta))
bal_after_psid1 <- gpbalancer::bal_table(data_obs_psid1, 2:9, as.logical(data_obs_psid1$ta), wts = wts_psid1)
bal_after_psid2 <- gpbalancer::bal_table(data_obs_psid2, 2:9, as.logical(data_obs_psid2$ta), wts = wts_psid2)
bal_after_psid3 <- gpbalancer::bal_table(data_obs_psid3, 2:9, as.logical(data_obs_psid3$ta), wts = wts_psid3)
bal_after_cps1 <- gpbalancer::bal_table(data_obs_cps1, 2:9, as.logical(data_obs_cps1$ta), wts = wts_cps1)
bal_after_cps2 <- gpbalancer::bal_table(data_obs_cps2, 2:9, as.logical(data_obs_cps2$ta), wts = wts_cps2)
bal_after_cps3 <- gpbalancer::bal_table(data_obs_cps3, 2:9, as.logical(data_obs_cps3$ta), wts = wts_cps3)

bal_total_exp <- mom_bal(data_experimental, 2:9, as.logical(data_experimental$ta))
bal_total_psid1 <- mom_bal(data_obs_psid1, 2:9, as.logical(data_obs_psid1$ta), wts = wts_psid1)
bal_total_psid2 <- mom_bal(data_obs_psid2, 2:9, as.logical(data_obs_psid2$ta), wts = wts_psid2)
bal_total_psid3 <- mom_bal(data_obs_psid3, 2:9, as.logical(data_obs_psid3$ta), wts = wts_psid3)
bal_total_cps1 <- mom_bal(data_obs_cps1, 2:9, as.logical(data_obs_cps1$ta), wts = wts_cps1)
bal_total_cps2 <- mom_bal(data_obs_cps2, 2:9, as.logical(data_obs_cps2$ta), wts = wts_cps2)
bal_total_cps3 <- mom_bal(data_obs_cps3, 2:9, as.logical(data_obs_cps3$ta), wts = wts_cps3)




ests <- c(est_ate_psid1$ests[2,1], est_ate_psid2$ests[2,1], est_ate_psid3$ests[2,1],
          est_ate_cps1$ests[2,1], est_ate_cps2$ests[2,1], est_ate_cps3$ests[2,1])
bals <- c(bal_total_psid1, bal_total_psid2, bal_total_psid3, bal_total_cps1, bal_total_cps2, bal_total_cps3)
idnames <- c('psid1', 'psid2', 'psid3',
             'cps1', 'cps2', 'cps3')

par(mfrow=c(1,1))
plot(bals, ests, pch=19, col=rgb(0,0,0,0.5),
     xlim = c(0, 3.5), ylim=range(ests))
text(bals, ests, idnames, pos = 4)
points(bal_total_exp, 1700, pch=4)
abline(h=0, lty=3)
abline(v=2)

balances <- cbind(bal_experimental[,7:8],
                  bal_after_psid1[,7:8],
                  bal_after_psid2[,7:8],
                  bal_after_psid3[,7:8],
                  bal_after_cps1[,7:8],
                  bal_after_cps2[,7:8],
                  bal_after_cps3[,7:8])
balances <- rbind(c(nrow(lalonde_control),
                    nrow(psid1_control), nrow(psid2_control), nrow(psid3_control),
                    nrow(cps1_control), nrow(cps2_control), nrow(cps3_control)),
                  balances)

rows_i_want <- c(0,1,2,3,6,9,12,15,16)+1
bals_subset <- balances[c(0,1,2,5,8,11,14,15,16)+1,]
rownames(bals_subset) <- rownames(balances)[rows_i_want]
rownames(bals_subset)[1] <- "Nobs"
colnames(bals_subset) <- c('Experimental',
                           'PSID1', 'PSID2', 'PSID3',
                           'CPS1', 'CPS2', 'CPS3')

mean_var <- function(dset, col){
    mean_val <- round(dset[,col],2)
    sd_val <- paste('(',round(sqrt(dset[,col+1]),2),')',sep='')
    cols <- cbind(mean_val, sd_val)
    stat_vals <- apply(cols, 1, paste, collapse=" ")
    stat_vals[c(1,2,5,8,11,14,15,16)]
}

table_one <- cbind(mean_var(bal_before_exp, 2),
                   mean_var(bal_before_exp, 5),
                   mean_var(bal_before_psid1, 5),
                   mean_var(bal_before_psid2, 5),
                   mean_var(bal_before_psid3, 5),
                   mean_var(bal_before_cps1, 5),
                   mean_var(bal_before_cps2, 5),
                   mean_var(bal_before_cps3, 5))
table_one <- rbind(c(nrow(lalonde_treated), nrow(lalonde_control),
                     nrow(psid1_control), nrow(psid2_control), nrow(psid3_control),
                     nrow(cps1_control), nrow(cps2_control), nrow(cps3_control)),
                   table_one)
rownames(table_one) <- rownames(balances)[rows_i_want]
rownames(table_one)[1] <- "Nobs"
colnames(table_one) <- c('NSW_t', 'NSW_c',
                         'PSID1', 'PSID2', 'PSID3',
                         'CPS1', 'CPS2', 'CPS3')
xtable(t(table_one))


befores <- round(cbind(bal_experimental[,7:8],
                       bal_before_psid1[,7:8],
                       bal_before_psid2[,7:8],
                       bal_before_psid3[,7:8],
                       bal_before_cps1[,7:8],
                       bal_before_cps2[,7:8],
                       bal_before_cps3[,7:8]),3)
befores <- befores[c(0,1,2,5,8,11,14,15,16),]
rownames(befores) <- rownames(bals_subset)[2:9]
xtable(befores)



afters <- round(cbind(bal_experimental[,7:8],
                      bal_after_psid1[,7:8],
                      bal_after_psid2[,7:8],
                      bal_after_psid3[,7:8],
                      bal_after_cps1[,7:8],
                      bal_after_cps2[,7:8],
                      bal_after_cps3[,7:8]),3)
afters <- afters[c(0,1,2,5,8,11,14,15,16),]

rownames(afters) <- rownames(bals_subset)[2:9]
meanbals <- afters[,seq(1,13,2)]
colnames(meanbals) <- c('NSW', 'PSID-1', 'PSID-2', 'PSID-3',
                        'CPS-1', 'CPS-2', 'CPS-3')
knitr::kable(meanbals)

xtable(afters)


apply(abs(afters), 2, mean)[seq(1,13,2)]

round(c(bal_total_exp, bals),2)


bal_start_exp <- mom_sq_bal(data_experimental, 2:9, as.logical(data_experimental$ta))
bal_start_psid1 <- mom_sq_bal(data_obs_psid1, 2:9, as.logical(data_obs_psid1$ta))
bal_start_psid2 <- mom_sq_bal(data_obs_psid2, 2:9, as.logical(data_obs_psid2$ta))
bal_start_psid3 <- mom_sq_bal(data_obs_psid3, 2:9, as.logical(data_obs_psid3$ta))
bal_start_cps1 <- mom_sq_bal(data_obs_cps1, 2:9, as.logical(data_obs_cps1$ta))
bal_start_cps2 <- mom_sq_bal(data_obs_cps2, 2:9, as.logical(data_obs_cps2$ta))
bal_start_cps3 <- mom_sq_bal(data_obs_cps3, 2:9, as.logical(data_obs_cps3$ta))


bals2 <- c(bal_start_exp, bal_start_psid1, bal_start_psid2, bal_start_psid3, bal_start_cps1, bal_start_cps2, bal_start_cps3)


bal_end_exp <- mom_sq_bal(data_experimental, 2:9, as.logical(data_experimental$ta))
bal_end_psid1 <- mom_sq_bal(data_obs_psid1, 2:9, as.logical(data_obs_psid1$ta), wts = wts_psid1)
bal_end_psid2 <- mom_sq_bal(data_obs_psid2, 2:9, as.logical(data_obs_psid2$ta), wts = wts_psid2)
bal_end_psid3 <- mom_sq_bal(data_obs_psid3, 2:9, as.logical(data_obs_psid3$ta), wts = wts_psid3)
bal_end_cps1 <- mom_sq_bal(data_obs_cps1, 2:9, as.logical(data_obs_cps1$ta), wts = wts_cps1)
bal_end_cps2 <- mom_sq_bal(data_obs_cps2, 2:9, as.logical(data_obs_cps2$ta), wts = wts_cps2)
bal_end_cps3 <- mom_sq_bal(data_obs_cps3, 2:9, as.logical(data_obs_cps3$ta), wts = wts_cps3)


bals3 <- c(bal_end_exp, bal_end_psid1, bal_end_psid2, bal_end_psid3, bal_end_cps1, bal_end_cps2, bal_end_cps3)


uhh <- rbind(est_experiment$ests[2,],
             est_ate_cps2$ests[2,],
             est_ate_cps1$ests[2,],
             est_ate_cps3$ests[2,],
             est_ate_psid3$ests[2,],
             est_ate_psid2$ests[2,],
             est_ate_psid1$ests[2,])
cbind(uhh, uhh[,4] - uhh[,3])

uhh <- rbind(est_experiment$ests[2,],
             est_ate_psid1$ests[2,],
             est_ate_psid2$ests[2,],
             est_ate_psid3$ests[2,],
             est_ate_cps1$ests[2,],
             est_ate_cps2$ests[2,],
             est_ate_cps3$ests[2,])
uhh
rownames(uhh) <- c('NSW', 'PSID-1', 'PSID-2', 'PSID-3',
                   'CPS-1', 'CPS-2', 'CPS-3')
knitr::kable(uhh)
round(uhh[,1],0)
round(uhh[,2],0)

cbind(uhh, uhh[,4] - uhh[,3])

hist(est_ps_psid2$ps[data_obs_psid2$ta==0])
hist(est_ps_psid2$ps[data_obs_psid2$ta==1],add=T)

hist(est_ps_psid3$ps[data_obs_psid3$ta==0])
hist(est_ps_psid3$ps[data_obs_psid3$ta==1],add=T)

hist(est_ps_cps3$ps[data_obs_cps3$ta==0])
hist(est_ps_cps3$ps[data_obs_cps3$ta==1],add=T)

plot( uhh[,1], c(1794, 1473, 1480, 1549, 1616, 1563, 662))
abline(0,1)

