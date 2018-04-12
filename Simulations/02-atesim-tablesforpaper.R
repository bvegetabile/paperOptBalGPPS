setwd('~/Dropbox/ucirvine/research/papers/2017_obgpps/2018-03-Simulations/ate-simresults/')
ateresults <- readRDS('2018-03-22-nonparametric_odd-atesim-results.rds')

true_ate <- 3

meanbal1 <- apply(abs(ateresults$Cov1_Balance) < 0.1 & abs(ateresults$Cov2_Balance) < 0.1, 2, mean)
meanbal15 <- apply(abs(ateresults$Cov1_Balance) < 0.15 & abs(ateresults$Cov2_Balance) < 0.15, 2, mean)
meanbal2 <- apply(abs(ateresults$Cov1_Balance) < 0.2 & abs(ateresults$Cov2_Balance) < 0.2, 2, mean)

# Bias
biases_lin <- ateresults$LinearResults - true_ate
niave_bias <- matrix(rep(biases_lin[,1], 12), nrow=nrow(biases_lin), ncol=12)
lin_avg_bias <- apply(biases_lin, 2, mean)
lin_bias_red <- apply(1 - abs(biases_lin) / abs(niave_bias), 2, mean) * 100
lin_avg_absbias <- apply(abs(biases_lin), 2, mean)
lin_emp_std_err <- apply(ateresults$LinearResults, 2, sd)
lin_emp_mse <- apply(biases_lin^2, 2, mean)

biases_em <- ateresults$EffModResults - true_ate
niave_bias <- matrix(rep(biases_em[,1], 12), nrow=nrow(biases_em), ncol=12)
em_avg_bias <- apply(biases_em, 2, mean)
em_bias_red <- apply(1 - abs(biases_em) / abs(niave_bias), 2, mean) * 100
em_avg_absbias <- apply(abs(biases_em), 2, mean)
em_emp_std_err <- apply(ateresults$EffModResults, 2, sd)
em_emp_mse <- apply(biases_em^2, 2, mean)


outro <- rbind(meanbal1, meanbal15, meanbal2,
               lin_avg_bias, lin_avg_absbias, lin_bias_red, lin_emp_std_err, lin_emp_mse,
               em_avg_bias, em_avg_absbias, em_bias_red, em_emp_std_err, em_emp_mse)
colnames(outro) <- c('NAIVE', 'TRUEPS', 'OBGPPS:NPSE', 'OBGPPS:SE',
                     'BART', 'GBM:KS.MEAN', 'GBM:ES.MEAN', 'GBM:ES.MAX',
                     'GLM:CORRECT', 'CBPS:CORRECT', 'GLM:MISSPECIFIED', 'CBPS:MISSPECIFIED')
xtable::xtable(t(outro), digits=3)


t(outro)


ateresults <- readRDS('2018-03-23-nonparametric_even-atesim-results.rds')

true_ate <- 3

meanbal1 <- apply(abs(ateresults$Cov1_Balance) < 0.1 & abs(ateresults$Cov2_Balance) < 0.1, 2, mean)
meanbal15 <- apply(abs(ateresults$Cov1_Balance) < 0.15 & abs(ateresults$Cov2_Balance) < 0.15, 2, mean)
meanbal2 <- apply(abs(ateresults$Cov1_Balance) < 0.2 & abs(ateresults$Cov2_Balance) < 0.2, 2, mean)

# Bias
biases_lin <- ateresults$LinearResults - true_ate
niave_bias <- matrix(rep(biases_lin[,1], 12), nrow=nrow(biases_lin), ncol=12)
lin_avg_bias <- apply(biases_lin, 2, mean)
lin_bias_red <- apply(1 - abs(biases_lin) / abs(niave_bias), 2, mean) * 100
lin_avg_absbias <- apply(abs(biases_lin), 2, mean)
lin_emp_std_err <- apply(ateresults$LinearResults, 2, sd)
lin_emp_mse <- apply(biases_lin^2, 2, mean)

biases_em <- ateresults$EffModResults - true_ate
niave_bias <- matrix(rep(biases_em[,1], 12), nrow=nrow(biases_em), ncol=12)
em_avg_bias <- apply(biases_em, 2, mean)
em_bias_red <- apply(1 - abs(biases_em) / abs(niave_bias), 2, mean) * 100
em_avg_absbias <- apply(abs(biases_em), 2, mean)
em_emp_std_err <- apply(ateresults$EffModResults, 2, sd)
em_emp_mse <- apply(biases_em^2, 2, mean)


outro <- rbind(meanbal1, meanbal15, meanbal2,
               lin_avg_bias, lin_avg_absbias, lin_bias_red, lin_emp_std_err, lin_emp_mse,
               em_avg_bias, em_avg_absbias, em_bias_red, em_emp_std_err, em_emp_mse)
colnames(outro) <- c('NAIVE', 'TRUEPS', 'OBGPPS:NPSE', 'OBGPPS:SE',
                     'BART', 'GBM:KS.MEAN', 'GBM:ES.MEAN', 'GBM:ES.MAX',
                     'GLM:CORRECT', 'CBPS:CORRECT', 'GLM:MISSPECIFIED', 'CBPS:MISSPECIFIED')
xtable::xtable(t(outro), digits=3)


t(outro)
