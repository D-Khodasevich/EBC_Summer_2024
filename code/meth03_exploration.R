.libPaths( c('.Rpackages',.libPaths() ) )

options(warn=-1)
library(data.table)
suppressMessages(library(limma))
suppressMessages(library(ENmix))
suppressMessages(library(sva))
suppressMessages(library(ChAMP))
options(warn=0)


#+ setdir01, echo = F
# knitr::opts_knit$set(root.dir = "../")
load("data/processed.rda")

#'# Exploring global DNA Methylation variability via PCs
# Let's look at the effect of sex
#+ fig.width=8, fig.height=6, dpi=300
plotMDS(beta, top=10000, gene.selection="common",
        pch=17,col=c("deeppink","blue")[factor(pheno$sex)],
        dim=c(1,2),cex=1.5)
legend("center", legend=levels(factor(pheno$sex)),bty='n',
       cex=1.5,pch=17,col=c("deeppink","blue"))

# Let's look at the effect of sex but let's remove the X and Y chromosomes
betas.clean = beta[manifest[probe_type=="cg" & !chr %in% c("X","Y")]$index,]
plotMDS(betas.clean, top=10000, gene.selection="common",
        pch=17,col=c("deeppink","blue")[factor(pheno$sex)],
        dim=c(1,2),cex=1.5)
legend("center", legend=levels(factor(pheno$sex)),bty='n',
       cex=1.5,pch=17,col=c("deeppink","blue"))
       
# Let's look at the effect of technical variables
# Row
cols <- rainbow(n = length(levels(factor(pheno$row))))
plotMDS(betas.clean, top=10000, gene.selection="common",
        pch=17,col=cols[factor(pheno$row)],
        dim=c(1,2),cex=1.5)
legend("top", legend=levels(factor(pheno$row)),bty='n',
       cex=1.5,pch=17,col=cols)
# If clustering by technical variables is present, batch effects can be adjusted for using ComBat or using covariates for batch in downstream analyses.
# Example of ComBat
# Impute missing Beta-values (ComBat will produce an error with missingness)
# sum(is.na(betas.clean))
# betas.impute = champ.impute(beta=betas.clean, pd=pheno, k=5, ProbeCutoff=0.2, SampleCutoff=0.1)
# betas.impute = betas.impute$beta
# Run ComBat
# batch <- factor(pheno$Sentrix_ID)
# modcombat <- model.matrix(~1, data = pheno)
# betaCombat <- ComBat(dat = as.matrix(betas.impute), batch = batch, mod = modcombat)


#' It would be useful to look at several traits with global variability
cov<-data.frame(pheno[,c('smoker','sex','CD4','CD8','NK','MO','GR','B')])
npc <- 20 # Top 20 PCs
svd <- prcomp(t(na.omit(betas.clean)))
screeplot(svd, npc, type = "barplot")
eigenvalue <- svd[["sdev"]]^2
prop <- (sum(eigenvalue[1:npc])/sum(eigenvalue)) * 100
cat("Top ", npc, " principal components can explain ", 
    prop, "% of data \n    variation", "\n")
screeplot(svd,npc,type="barplot")
pcrplot(na.omit(betas.clean), cov, npc=10) # Already saved in working directory


#'# Epigenetic Age
#'Load package for Age-Prediction
options(warn=-1)
suppressMessages(library(wateRmelon))
options(warn=0)
DNAmAge<-agep(beta)$horvath.age
hist(DNAmAge)
boxplot(DNAmAge);stripchart(DNAmAge, vertical = T,method = "jitter", add = T, pch = 20, col = 'red')

#'Agreement between Horvath's Epigenetic age and Hannum's clock
data(hannumCoef)
length(hannumCoef)
DNAmAge.Hannum<-agep(beta,coeff=hannumCoef,method = "hannum")$custom_age


#' Correlation; agreement
cor.test(DNAmAge,DNAmAge.Hannum)
plot(DNAmAge.Hannum,DNAmAge,pch=21,ylab="Horvath's DNAm Age",
     xlab="Hannum's DNAm Age",cex=1.2, bg=alpha("deepskyblue",0.45),main="Epigenetic Clocks")
legend("topleft",legend=paste0("r=",round(cor(DNAmAge.Hannum,DNAmAge),2)),bty="n")
abline(lm(DNAmAge~DNAmAge.Hannum),col="red",lw=2)


#' Age Acceleration Residuals
AgeAccelerationResidual<-residuals(lm(DNAmAge.Hannum~DNAmAge))
boxplot(AgeAccelerationResidual ~pheno$smoker, col=c("blue","red"))
wilcox.test(AgeAccelerationResidual ~ pheno$smoker)

#' Differences by Smoking status
boxplot(DNAmAge ~pheno$smoker, col=c("blue","red"))
wilcox.test(DNAmAge ~ pheno$smoker)
boxplot(DNAmAge.Hannum ~pheno$smoker, col=c("blue","red"))
wilcox.test(DNAmAge.Hannum ~ pheno$smoker)

#' Load Age Acceleration measures
#' See  [Horvath New Methylation Age Calculator](http://dnamage.genetics.ucla.edu/).  
Online<-read.csv("Clock/datout_New.output.csv")
Online<- Online[match(pheno$gsm, Online$SID),]
#+ fig.width=8, fig.height=6, dpi=300
group <- NA
group[pheno$smoker=="smoker"] <- 1
group[pheno$smoker=="non-smoker"] <- 2
pairs(~DNAmAge + DNAmPhenoAge + DNAmAgeSkinBloodClock + 
        DNAmGrimAgeBasedOnRealAge + DNAmTL + DNAmAgeHannum, 
    col = c("red","blue")[group],
    pch = c(8, 18)[group], 
    data = Online)


#' Look at GrimAge Acceleration
boxplot(Online$AgeAccelGrim ~pheno$smoker, col=c("blue","red"))
wilcox.test(Online$AgeAccelGrim ~ pheno$smoker)

#' Pack years predicted from DNA methylation
boxplot(Online$DNAmPACKYRS ~pheno$smoker, col=c("blue","red"))
wilcox.test(Online$DNAmPACKYRS ~ pheno$smoker)

#' DNA methylation estimate of TL
boxplot(Online$DNAmTLAdjAge ~pheno$smoker, col=c("blue","red"))
t.test(Online$DNAmTLAdjAge ~ pheno$smoker)

#' End of script 03