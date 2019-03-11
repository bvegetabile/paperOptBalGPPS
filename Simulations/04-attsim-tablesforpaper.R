setwd('~/Documents/GitHub/paperOptBalGPPS/Simulations/att-simresults/')

attresults <- readRDS('2019-02-11-attsim-multivariate-results-revision.rds')

calcbal <- function(res, thresh = 0.1){
  bal1 <- abs(res$Cov1_Balance) < thresh
  bal2 <- abs(res$Cov2_Balance) < thresh
  bal3 <- abs(res$Cov3_Balance) < thresh
  bal4 <- abs(res$Cov4_Balance) < thresh
  bal5 <- abs(res$Cov5_Balance) < thresh
  
  apply(bal1 & bal2 & bal3 & bal4 & bal5, 2, mean)
}


true_att_em <- mean(attresults$EffModResults[,1])
true_att_lin <- 10

# Balance Summaries
meanbal1 <- calcbal(attresults, 0.1)
meanbal15 <- calcbal(attresults, 0.15)
meanbal2 <- calcbal(attresults, 0.2)

# Bias
biases_lin <- attresults$LinearResults - true_att_lin
lin_avg_bias <- apply(biases_lin, 2, mean)
lin_avg_absbias <- apply(abs(biases_lin), 2, mean)
lin_emp_std_err <- apply(attresults$LinearResults, 2, sd)
lin_emp_mse <- apply(biases_lin^2, 2, mean)

biases_em <- attresults$EffModResults - true_att_em
em_avg_bias <- apply(biases_em, 2, mean)
em_avg_absbias <- apply(abs(biases_em), 2, mean)
em_emp_std_err <- apply(attresults$EffModResults, 2, sd)
em_emp_mse <- apply(biases_em^2, 2, mean)


outro <- rbind(meanbal1, meanbal15, meanbal2,
               lin_avg_bias, lin_avg_absbias, lin_emp_std_err, lin_emp_mse,
               em_avg_bias, em_avg_absbias, em_emp_std_err, em_emp_mse)
colnames(outro) <- c('TRUTH', "TRUEPS",
                     'OBGPPS:NPSE1:EP', 'OBGPPS:NPSE2:EP', 'OBGPPS:SE:EP',
                     'BART', 'GBM:KS.MEAN', 'GBM:ES.MEAN', 'GBM:ES.MAX',
                     'CBPS', 'GLM:POLY2', 'GLM:POLY1')
t(outro)[2:12,]

xtable::xtable(t(outro)[2:12,], digits=3)

