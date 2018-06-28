# Testing SPM for multiple imputations #
# One-dimensional prediction #
library(stpm)
library(grid)
library(ggplot2)
library(openxlsx)



#source("Z:/data/work/iz12/spm/DataImputation/scripts/multiplot.R")
source("~/Projects/spm/DataImputation/scripts/multiplot.R")

#######################################################
############## One dimensional case ###################
#######################################################

nexp <- 100 # Number of experiments

# Store parameter estimates for true data
spm.results <- data.frame(matrix(ncol=7, nrow=nexp))
colnames(spm.results) <- c("a", "f1", "Q", "f", "b", "mu0", "theta")
# Store parameter estimates for the imputed data
spm.results.imputed <- data.frame(matrix(ncol=7, nrow=nexp))
spm.results.incomplete <- data.frame(matrix(ncol=7, nrow=nexp))
colnames(spm.results.imputed) <- c("a", "f1", "Q", "f", "b", "mu0", "theta")
colnames(spm.results.incomplete) <- c("a", "f1", "Q", "f", "b", "mu0", "theta")

diff.table <- data.frame(matrix(NA,nrow=1, ncol = 2))
colnames(diff.table) <- c("AD", "RMSE")
rownames(diff.table) <- c("stpm")



N <- 1000 # Number of individuals in a dataset
nd <- 1
n.imp <- 5
theta_range <- seq(0.01, 0.15, by=0.001)
#pdir <- "/Volumes/G/Projects/spm_mi/"
#pdir <- "Z:/data/work/iz12/spm/spm.impute/"
pdir <- "~/Projects/spm/spm.impute/"
dir.create(paste(pdir, "simulated/", sep=""))

for(jj in 1:nexp) {
    print(paste("Experiment:", jj))
    #--------------- stpm ---------------#
    ############ Data preparation ##############
    data <- simdata_discr(N=N, dt = 2, format="short")
    #write.table(x=data, file=paste(pdir, "simulated/", "data.", nd, ".", N, "_", jj, ".txt", sep=""), row.names = FALSE, col.names = TRUE)
    
    
    miss.id <- sample(x=dim(data)[1], size=round(dim(data)[1]/4)) # ~25% missing data
    incomplete.data <- data
    incomplete.data[miss.id,4] <- NA
    #write.table(x=incomplete.data, file=paste(pdir, "simulated/", "incomplete.data.", nd, ".", N, "_", jj, ".txt", sep=""), row.names = FALSE, col.names = TRUE)
    # End of data preparation #
  
    ##### Multiple imputation with SPM #####
    imp.data <- spm.impute(x=incomplete.data, id=1, case="xi", t1=3, covariates="y1", minp=n.imp, theta_range=theta_range)$imputed
    # Estimate parameters from the 'true' dataset #
    ddd <- prepare_data(x=data, col.id="id", col.status = "xi", col.age = "t", col.age.event = "t", covariates = "y1", interval=1)
    p <- spm_discrete(ddd[[1]], theta_range = theta_range)
    spm.results$a[jj] <- p$cmodel$a
    spm.results$f1[jj] <- p$cmodel$f1
    spm.results$Q[jj] <- p$cmodel$Q
    spm.results$f[jj] <- p$cmodel$f
    spm.results$b[jj] <- p$cmodel$b
    spm.results$mu0[jj] <- p$cmodel$mu0
    spm.results$theta[jj] <- p$cmodel$theta
    
     
    # Estimate SPM parameters from imputed data and compare them to the p:
    ddd.imp <- prepare_data(x=imp.data, col.id="id", col.status = "xi", col.age = "t", col.age.event = "t", covariates = "y1")
    pp.test <- spm_discrete(ddd.imp[[1]], theta_range = theta_range)
    spm.results.imputed$a[jj] <- pp.test$cmodel$a
    spm.results.imputed$f1[jj] <- pp.test$cmodel$f1
    spm.results.imputed$Q[jj] <- pp.test$cmodel$Q
    spm.results.imputed$f[jj] <- pp.test$cmodel$f
    spm.results.imputed$b[jj] <- pp.test$cmodel$b
    spm.results.imputed$mu0[jj] <- pp.test$cmodel$mu0
    spm.results.imputed$theta[jj] <- pp.test$cmodel$theta
    
    # -------------------------------------Writing data-------------------------------------#
    #write.table(x=imp.data, file=paste(pdir, "simulated/", "imp.data.spm.1d.", N, "_", jj, ".txt", sep=""),col.names=TRUE,row.names=FALSE)
    
    ### Calculating the differences: AD and RMSE ###
    # AD
    ad.spm <- mean(abs(imp.data[,4] - data[, 4]), na.rm=T) # Average absolute distance
    
    # RMSE
    rmse.spm <- sqrt(mean((imp.data[,4] - data[, 4])^2, na.rm = T))
    
    diff.table["stpm","AD"] <- ad.spm
    diff.table["stpm","RMSE"] <- rmse.spm
    
    #write.table(diff.table, file=paste(pdir, "simulated/", "results.1d.", N, "_", jj,".txt", sep=""),row.names = T, col.names = T)
    #write.xlsx(diff.table, file=paste(pdir, "simulated/", "results.1d.", N, "_", jj, ".xlsx", sep=""),row.names = T, col.names = T)
}

write.table(x=spm.results, file=paste(pdir, "simulated/", "spm.results.", N, ".tsv", sep=""), row.names = FALSE, col.names = TRUE)
write.table(x=spm.results.imputed, file=paste(pdir, "simulated/", "spm.results.imputed.", N, ".tsv", sep=""), row.names = FALSE, col.names = TRUE)


#-------------------------------- Results processing ---------------------------------#
#### SPM ####
# Combine your two dataframes into one.  First make a new column in each.
True <- spm.results
Imputed <- spm.results.imputed
True$test <- 'True'
Imputed$test <- 'Imputed'
# And then combine into your new data frame vegLengths
comb <- rbind(True, Imputed)

a.pval <- round(t.test(spm.results$a, spm.results.imputed$a)$p.value, 3)
f1.pval <- round(t.test(spm.results$f1, spm.results.imputed$f1)$p.value, 3)
Q.pval <- round(t.test(spm.results$Q, spm.results.imputed$Q)$p.value, 3)
f.pval <- round(t.test(spm.results$f, spm.results.imputed$f)$p.value, 3)
mu0.pval <- round(t.test(spm.results$mu0, spm.results.imputed$mu0)$p.value, 3)
b.pval <- round(t.test(spm.results$b, spm.results.imputed$b)$p.value, 3)
theta.pval <- round(t.test(spm.results$theta, spm.results.imputed$theta)$p.value, 3)

p1 <- ggplot(comb, aes(a, fill = test)) + geom_density(alpha=.2) + scale_fill_discrete("") + labs(x=(paste("a (p.value=", a.pval, ")")))
p2 <- ggplot(comb, aes(f1, fill = test)) + geom_density(alpha = 0.2) + scale_fill_discrete("") + labs(x=(paste("f1 (p.value=", f1.pval, ")")))
p3 <- ggplot(comb, aes(Q, fill = test)) + geom_density(alpha = 0.2) + scale_fill_discrete("") + labs(x=(paste("Q (p.value=", Q.pval, ")")))
p4 <- ggplot(comb, aes(f, fill = test)) + geom_density(alpha = 0.2) + scale_fill_discrete("") + labs(x=(paste("f (p.value=", f.pval, ")")))
p5 <- ggplot(comb, aes(b, fill = test)) + geom_density(alpha = 0.2) + scale_fill_discrete("") + labs(x=(paste("b (p.value=", b.pval, ")")))
p6 <- ggplot(comb, aes(mu0, fill = test)) + geom_density(alpha = 0.2) + scale_fill_discrete("") + labs(x=(paste("mu0 (p.value=", mu0.pval, ")")))
p7 <- ggplot(comb, aes(theta, fill = test)) + geom_density(alpha = 0.2) + scale_fill_discrete("") + labs(x=(paste("theta (p.value=", theta.pval, ")")))

pdf(file=paste(pdir, "/true_vs_imputed_spm.", nd, ".pdf", sep=""), width = 16, height = 12)
    multiplot(p1, p2, p3, p4, p5, p6, p7, cols = 3, rows=3, title = "Parameters estimated from true against imputed data with method SPM")
dev.off()
