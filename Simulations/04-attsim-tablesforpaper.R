setwd('~/Dropbox/ucirvine/research/papers/2017_obgpps/2018-03-Simulations/att-simresults/')
attresults <- readRDS('2018-03-21-attsim-mixtures-results.rds')

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

pdf('att-simulationsetup.pdf', height=4, width=12)
testexes <- seq(-10, 10, length.out = 1000)
densXt <- p1_t * dnorm(testexes, m1_t, s1_t) + (1 - p1_t) * dnorm(testexes, m2_t, s2_t)
densXc <- p1_c * dnorm(testexes, m1_c, s1_c) + (1 - p1_c) * dnorm(testexes, m2_c, s2_c)
p_tgivenx <- prob_t * densXt / (prob_t * densXt + (1-prob_t) * densXc)
par(mfrow=c(1,3))
plot(testexes, densXt, ylim=c(0, max(c(densXt, densXc))), type='l',
     lwd=3, col = rgb(0.5,0,0,0.5),
     xlab=expression(X[1]), ylab='Density',
     main=expression(paste('Density Function of ', X[1], '\n in the Treated Group',sep='')))
plot(testexes, densXc, ylim=c(0, max(c(densXt, densXc))), type='l',
     lwd=3, col = rgb(0,0,0.5,0.5),
     xlab=expression(X[1]), ylab='Density',
     main=expression(paste('Density Function of ', X[1], '\n in the Control Group',sep='')))
plot(testexes, p_tgivenx, type='l',
     ylim=c(0,1),
     lwd=3, col = rgb(0,0,0,0.5),
     xlab=expression(X[1]), ylab=expression("Pr" * (T ==1 ~ "|" ~ X)),
     main=expression(paste('True Probability of Treatment',sep='')))
dev.off()

true_att_em
true_att_lin

# Balance Summaries
meanbal1 <- apply(abs(attresults$Cov1_Balance) < 0.1 & abs(attresults$Cov2_Balance) < 0.1, 2, mean)
meanbal15 <- apply(abs(attresults$Cov1_Balance) < 0.15 & abs(attresults$Cov2_Balance) < 0.15, 2, mean)
meanbal2 <- apply(abs(attresults$Cov1_Balance) < 0.2 & abs(attresults$Cov2_Balance) < 0.2, 2, mean)

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
colnames(outro) <- c('TRUTH', 'OBGPPS:NPSE', 'OBGPPS:SE', 'BART', 'GBM:KS.MEAN', 'GBM:ES.MEAN', 'GBM:ES.MAX',
                     'CBPS', 'GLM:POLY2', 'GLM:POLY1')
t(outro)

xtable::xtable(t(outro), digits=3)

