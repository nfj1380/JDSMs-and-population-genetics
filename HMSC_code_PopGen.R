#-----------------------------------------------------------------------------
#################Joint species distributionmodels (JSDM) by Nick Fountain-Jones#################
#-----------------------------------------------------------------------------

#prelims
rm(list = ls())
library(HMSC)
set.seed(29698)
library(extendedForest)
library(gradientForest) #needs extended forest. Both are available here:https://r-forge.r-project.org/R/?group_id=973



#-----------------------------------------------------------------------------
#################Read Data#################
#-----------------------------------------------------------------------------

# Community matrix
Y <- read.csv("GammSNPsWS.csv", head= T, row.names=1)
#make sure data is binary
Y <- ifelse(Y > 0, 1, 0)
# Covariates (no NAs  -must be numeric)
X <- read.csv("EnvGammaWS.csv", head= T, row.names=1);str(X)

if(colnames(X)==colnames(Y))
{
  print("identical")
}
#select variables to add to the model
library(dplyr)
SubX <- X %>% select( 'Year','Hunting', 'EuclideanDIstPCA1', 'EuclideanDIstPCA2')
SubX$Hunting <-as.factor(SubX$Hunting);str(SubX)

#check correlations
library(GGally)
ggpairs(SubX, aes(alpha = 0.4))

# Random effects
Pi <- read.csv("PiGammaWS.csv", head= T, row.names=1);str(Pi)

# Covert all columns of Pi to a factor
Pi <- data.frame(apply(Pi,2,as.factor));str(Pi)

#Turn into matrices for HMSC 
YMat<-as.matrix(Y)
XMat<-model.matrix(~.,data=SubX)



#-----------------------------------------------------------------------------
#################Gradient Forests################# 
#-----------------------------------------------------------------------------
#useful for larger datasets

nSites <- dim(Y)[1]
nSpecs <- dim(Y)[2]

#work out how many predictors should be used at each node (mtry)
lev <- floor(log2(nSites * 0.368/2))
lev
#run analysis
gf <- gradientForest(cbind(SubX, Y), predictor.vars = colnames(SubX), response.vars = colnames(Y), ntree = 1000, transform = NULL, compact = T,nbin = 201, maxLevel = lev, corr.threshold = 0.5)#gf doesn't handle factors
gf
#variable importance
most_important <- names(importance(gf))[1:25]
par(mgp = c(2, 2, 2))
#plots
par(mfrow = c(1,1))
plot(gf, plot.type = "O")
plot(gf, plot.type = "C", imp.vars = most_important, show.species = F, common.scale = T, cex.axis = 0.6, cex.lab = 0.7, line.ylab = 0.9, par.args = list(mgp = c(1.5,0.5, 0), mar = c(2.5, 1, 0.1, 0.5), omi = c(0, 0.3, 0, 0)))
plot(gf, plot.type = "C", imp.vars = most_important, show.overall = F, legend = T, leg.posn = "topleft",leg.nspecies = 10, cex.lab = 0.7, cex.legend = 0.9,cex.axis = 0.6, line.ylab = 0.9, par.args = list(mgp = c(1.5,0.5, 0), mar = c(2.5, 1, 0.1, 0.5), omi = c(0, 0.3, 0, 0)))
plot(gf, plot.type = "P", show.names = T, horizontal = F,cex.axis = 1, cex.labels = 0.7, line = 2.5)

dev.off()
#-----------------------------------------------------------------------------
#################Create HMSC objects and set priors#################
#-----------------------------------------------------------------------------


#get data into right format
formdata <- as.HMSCdata(Y = YMat, X = XMat, Random = Pi, interceptX = FALSE,
                         scaleX=TRUE)


formprior <- as.HMSCprior(formdata)
formparam <- as.HMSCparam(formdata, formprior)
#-----------------------------------------------------------------------------
#################Run the HMSC model#################
#-----------------------------------------------------------------------------
#run the model  
modelWSGamma<- hmsc(formdata , family = "probit", niter = 60000, nburn = 6000,
               thin = 250)
#save
save(modelWSGamma,  file = "modelWSGamma .RData")
load(".RData")

#check model
mixing4 <- as.mcmc(modelWSGamma , parameters = "paramX")

par(mfrow=c(1,1))
#check diagnostics (ESS - uses coda)
effectiveSize(mixing4)

### Convert the mixing object to a matrix
mixingDF <- as.data.frame(mixing4)
str(mixingDF)
### Draw beanplots
library(beanplot)
par(mar = c(7, 4, 4, 2))
beanplot(mixingDF, las = 2)

library(ggmcmc)
full <- ggs(mixing4)
p <- ggs_caterpillar(full)
p + geom_point(size = 0.01, stroke = 0, shape = 16)+
  theme(text = element_text(size=5),
        axis.text.x = element_text(angle=0, hjust=1))  
### Draw boxplot for each parameters
par(mar = c(7, 4, 4, 2))
boxplot(mixingDF, las = 2)
### Summary table
average <- apply(modelWSGamma $results$estimation$paramX, 1:2, mean)
CI.025 <- apply(modelWSGamma $results$estimation$paramX, 1:2, quantile, probs=0.025)
CI.975 <- apply(modelWSGamma $results$estimation$paramX, 1:2, quantile, probs=0.975)
CI <- cbind(as.vector(CI.025), as.vector(CI.975))

paramXCITable <- cbind(unlist(as.data.frame(average)),
                       unlist(as.data.frame(CI.025)),
                       unlist(as.data.frame(CI.975)))
colnames(paramXCITable) <- c("paramX", "lowerCI", "upperCI")
rownames(paramXCITable) <- paste(rep(colnames(average),
                                     each = nrow(average)), "_",
                                 rep(rownames(average),
                                     ncol(average)), sep="")
write.csv(paramXCITable, file='WSGamaaModel1.csv')
 
#-----------------------------------------------------------------------------
#################variation partioning#################
#-----------------------------------------------------------------------------
#group covariates
head(SubX)
groupCov <- c(rep ("Temporal", 2),rep("Spatial", 2))
variationPart <- variPart(modelWSGamma ,groupX= c("Intercept",groupCov))
par(mar=c(5,5,5,1))
barplot(t(variationPart), las=2, col= rainbow(5), cex.names=0.75, cex.axis=0.75,
        legend.text=paste(colnames(variationPart)," ",
                          signif(100*colMeans(variationPart),2),"%",sep=""),
        args.legend=list(y=1.2, xjust=1, horiz=F, bty="n",cex=0.75))
#R2 

# Prevalence
prevSp <- colSums(Y)
# Coefficient of multiple determination
R2 <- Rsquared(modelWSGamma , averageSp = FALSE)
R2comm <- Rsquared(modelWSGamma , averageSp = TRUE)
par(mar=c(5,6,0.5,0.5))
plot(prevSp, R2, xlab = "Prevalence",
     ylab = expression(R^2), pch=19, las=1,cex.lab = 2)
abline(h = R2comm, col = "blue", lwd = 2)

#-----------------------------------------------------------------------------
#################Testing for associations#################
#-----------------------------------------------------------------------------
assoMat <- corRandomEff(modelWSGamma)
library(corrplot)
library(circlize)
# Average
IndividualMean <- apply(assoMat[, , , 1], 1:2, mean)

#=======================
### Associations to draw
#=======================
#--------------------
### Individual level effect
#--------------------
# Build matrix of colours for chordDiagram
IndivDrawCol <- matrix(NA, nrow = nrow(IndividualMean),
                      ncol = ncol(IndividualMean))
IndivDrawCol[which(IndividualMean > 0.8, arr.ind=TRUE)]<-"red"
IndivDrawCol[which(IndividualMean < -0.8, arr.ind=TRUE)]<-"blue"
# Build matrix of "significance" for corrplot
IndivDraw <- IndivDrawCol
IndivDraw[which(!is.na(IndivDraw), arr.ind = TRUE)] <- 0
IndivDraw[which(is.na(IndivDraw), arr.ind = TRUE)] <- 1
IndivDraw <- matrix(as.numeric(IndivDraw), nrow = nrow(IndividualMean),
                   ncol = ncol(IndividualMean))
par(mfrow=c(1,1))
# Matrix plot
Colour <- colorRampPalette(c("blue", "white", "red"))(200)
corrplot(IndividualMean, method = "color", col = Colour, type = "lower",
         diag = FALSE, p.mat = IndivDraw, tl.srt = 45)
# Chord diagram
chordDiagram(IndividualMean, symmetric = TRUE,
             annotationTrack = c("name", "grid"),
             grid.col = "grey", col = IndivDrawCol)
write.csv(paramXCITable, file='ChannelTest1.csv')



#-------------------------
#Latent effectsModel - work in progress
#-------------------------
assoMat <- corRandomEff(model6)
LatentMean <- apply(assoMat[, , , 1], 1:2, mean)
# Build matrix of colours for chordDiagram
LatentDrawCol <- matrix(NA, nrow = nrow(LatentMean),
                      ncol = ncol(LatentMean))
LatentDrawCol[which(LatentMean > 0.4, arr.ind=TRUE)]<-"red"
LatentDrawCol[which(LatentMean < -0.4, arr.ind=TRUE)]<-"blue"
# Build matrix of "significance" for corrplot
LatentDraw <- LatentDrawCol
LatentDraw[which(!is.na(LatentDraw), arr.ind = TRUE)] <- 0
LatentDraw[which(is.na(LatentDraw), arr.ind = TRUE)] <- 1
LatentDraw <- matrix(as.numeric(LatentDraw), nrow = nrow(LatentMean),
                   ncol = ncol(LatentMean))
par(mfrow=c(1,2))
# Matrix plot
Colour <- colorRampPalette(c("blue", "white", "red"))(200)
corrplot(LatentMean, method = "color", col = Colour, type = "lower",
         diag = FALSE, p.mat = LatentDraw, tl.srt = 45)
# Chord diagram
chordDiagram(YearMean, symmetric = TRUE,
             annotationTrack = c("name", "grid"),
             grid.col = "grey", col = LatentDrawCol)