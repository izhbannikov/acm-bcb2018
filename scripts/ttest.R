library(grid)
library(ggplot2)

source("~/Projects/spm/DataImputation/scripts/multiplot.R")

#-------------------------------- Results processing ---------------------------------#
#### SPM ####
# Combine your two dataframes into one.  First make a new column in each.
spm.results <- read.table("~/Projects/spm/spm.impute/simulated/spm.results.1000.tsv", header=T)

spm.results.imputed <- read.table("~/Projects/spm/spm.impute/simulated/spm.results.imputed.1000.tsv", header=T)

#spm.results.imputed <- read.table("~/Projects/spm/spm.impute/simulated/spm.results.imputed.mice.pmm.1000.100.tsv", header=T)

#spm.results.imputed <- read.table("~/Projects/spm/spm.impute/simulated/spm.results.imputed.missforest.1000.100.tsv", header=T)


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

p1 <- ggplot(comb, aes(a, fill = test)) + geom_density(alpha=.2) + scale_fill_discrete("") + labs(x=(paste("a (p.value=", a.pval, ")"))) + 
  theme(text = element_text(size=14)) 
p2 <- ggplot(comb, aes(f1, fill = test)) + geom_density(alpha = 0.2) + scale_fill_discrete("") + labs(x=(paste("f1 (p.value=", f1.pval, ")"))) +
  theme(text = element_text(size=14)) 
p3 <- ggplot(comb, aes(Q, fill = test)) + geom_density(alpha = 0.2) + scale_fill_discrete("") + labs(x=(paste("Q (p.value=", Q.pval, ")"))) +
  theme(text = element_text(size=14)) 
p4 <- ggplot(comb, aes(f, fill = test)) + geom_density(alpha = 0.2) + scale_fill_discrete("") + labs(x=(paste("f (p.value=", f.pval, ")"))) +
  theme(text = element_text(size=14))
p5 <- ggplot(comb, aes(b, fill = test)) + geom_density(alpha = 0.2) + scale_fill_discrete("") + labs(x=(paste("b (p.value=", b.pval, ")"))) +
  theme(text = element_text(size=14))
p6 <- ggplot(comb, aes(mu0, fill = test)) + geom_density(alpha = 0.2) + scale_fill_discrete("") + labs(x=(paste("mu0 (p.value=", mu0.pval, ")"))) +
  theme(text = element_text(size=14))
p7 <- ggplot(comb, aes(theta, fill = test)) + geom_density(alpha = 0.2) + scale_fill_discrete("") + labs(x=(paste("theta (p.value=", theta.pval, ")"))) +
  theme(text = element_text(size=14))

#pdf(file=paste(pdir, "/true_vs_imputed_spm.", nd, ".pdf", sep=""), width = 16, height = 12)
#multiplot(p1, p2, p3, p4, p5, p6, p7, cols = 3, rows=3, title = "Parameters estimated from true against imputed data with method SPM")
#dev.off()

multiplot(p1, p2, p3, p4, p5, p6, p7, cols = 2, rows=4)


