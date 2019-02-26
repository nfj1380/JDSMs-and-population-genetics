#-----------------------------------------------------------------------------
#################HMSC Analysis adapted by Nick Fountain-Jones#################
#-----------------------------------------------------------------------------

#prelims
rm(list = ls())
library(HMSC)
set.seed(29698)

#-----------------------------------------------------------------------------
#################Read Data#################
#-----------------------------------------------------------------------------

# Community matrix
Y <- read.csv("GammSNPsWS.csv", head= T, row.names=1)
#make sure data is binary
Y <- ifelse(Y > 0, 1, 0)
# Covariates (no NAs  -must be numeric)
X <- read.csv("EnvGammaWS.csv", head= T, row.names=1);str(X)
#select variables to add to the model
library(dplyr)
SubX <- X %>% select( 'Year', 'Hunting', 'EuclideanDIstPCA1', 'EuclideanDIstPCA2', 'rdsPCO1', 'rdsPCO2')


#check correlations
library(Hmisc)
res <-rcorr(XMat, type=c("spearman"))
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
flattenCorrMatrix(res$r, res$P)

# Random effects
Pi <- read.csv("PiGammaWS.csv", head= T, row.names=1);str(Pi)

# Covert all columns of Pi to a factor
Pi <- data.frame(apply(Pi,2,as.factor));str(Pi)



#Turn into matrices 
YMat<-as.matrix(Y)
XMat<-model.matrix(~.,data=SubX)

if(colnames(X)==colnames(Y))
{
  print("identical")
  }
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
modelWSGamma<- hmsc(formdata , family = "probit", niter = 30000, nburn = 3000,
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
average <- apply(model1$results$estimation$paramX, 1:2, mean)
CI.025 <- apply(model1$results$estimation$paramX, 1:2, quantile, probs=0.025)
CI.975 <- apply(model1$results$estimation$paramX, 1:2, quantile, probs=0.975)
CI <- cbind(as.vector(CI.025), as.vector(CI.975))

paramXCITable <- cbind(unlist(as.data.frame(average)),
                       unlist(as.data.frame(CI.025)),
                       unlist(as.data.frame(CI.975)))
colnames(paramXCITable) <- c("paramX", "lowerCI", "upperCI")
rownames(paramXCITable) <- paste(rep(colnames(average),
                                     each = nrow(average)), "_",
                                 rep(rownames(average),
                                     ncol(average)), sep="")
write.csv(paramXCITable, file='CoinfectionModel1.csv')
 
#-----------------------------------------------------------------------------
#################variation partioning#################
#-----------------------------------------------------------------------------
#group covariates
head(SubX)
groupCov <- c(rep ("Temporal", 2),rep("Spatial", 2), rep("Roads", 2))
variationPart <- variPart(modelWSGamma ,groupX= c("Intercept",groupCov))
par(mar=c(5,5,5,1))
barplot(t(variationPart), las=2, col= rainbow(9), cex.names=0.75, cex.axis=0.75,
        legend.text=paste(colnames(variationPart)," ",
                          signif(100*colMeans(variationPart),2),"%",sep=""),
        args.legend=list(y=1.2, xjust=1, horiz=F, bty="n",cex=0.75))
#R2 

# Prevalence
prevSp <- colSums(Y)
# Coefficient of multiple determination
R2 <- Rsquared(modelprev5BabFIV , averageSp = FALSE)
R2comm <- Rsquared(modelprev5BabFIV , averageSp = TRUE)
par(mar=c(5,6,0.5,0.5))
plot(prevSp, R2, xlab = "Prevalence",
     ylab = expression(R^2), pch=19, las=1,cex.lab = 2)
abline(h = R2comm, col = "blue", lwd = 2)

#-----------------------------------------------------------------------------
#################Testing for associations#################
#-----------------------------------------------------------------------------
assoMat <- corRandomEff(modelMay2018traits)
library(corrplot)
library(circlize)
# Average
IndividualMean <- apply(assoMat[, , , 1], 1:2, mean)
YearMean <- apply(assoMat[, , , 3], 1:2, mean)
PrideYearMean<- apply(assoMat[, , , 2], 1:2, mean)
#=======================
### Associations to draw
#=======================
#--------------------
### Individual level effect
#--------------------
# Build matrix of colours for chordDiagram
IndivDrawCol <- matrix(NA, nrow = nrow(IndividualMean),
                      ncol = ncol(IndividualMean))
IndivDrawCol[which(IndividualMean > 0.4, arr.ind=TRUE)]<-"red"
IndivDrawCol[which(IndividualMean < -0.4, arr.ind=TRUE)]<-"blue"
# Build matrix of "significance" for corrplot
IndivDraw <- IndivDrawCol
IndivDraw[which(!is.na(IndivDraw), arr.ind = TRUE)] <- 0
IndivDraw[which(is.na(IndivDraw), arr.ind = TRUE)] <- 1
IndivDraw <- matrix(as.numeric(IndivDraw), nrow = nrow(IndividualMean),
                   ncol = ncol(IndividualMean))
par(mfrow=c(1,2))
# Matrix plot
Colour <- colorRampPalette(c("blue", "white", "red"))(200)
corrplot(IndividualMean, method = "color", col = Colour, type = "lower",
         diag = FALSE, p.mat = IndivDraw, tl.srt = 45)
# Chord diagram
chordDiagram(IndividualMean, symmetric = TRUE,
             annotationTrack = c("name", "grid"),
             grid.col = "grey", col = IndivDrawCol)
write.csv(paramXCITable, file='ChannelTest1.csv')

#--------------------
### Pride-Year level effect
#--------------------
# Build matrix of colours for chordDiagram
PrideYearDrawCol <- matrix(NA, nrow = nrow(PrideYearMean),
                       ncol = ncol(PrideYearMean))
PrideYearDrawCol[which(PrideYearMean > 0.3, arr.ind=TRUE)]<-"red"
IndivDrawCol[which(PrideYearMean < -0.3, arr.ind=TRUE)]<-"blue"
# Build matrix of "significance" for corrplot
PrideYearDraw <- PrideYearDrawCol
PrideYearDraw[which(!is.na(PrideYearDraw), arr.ind = TRUE)] <- 0
PrideYearDraw[which(is.na(PrideYearDraw), arr.ind = TRUE)] <- 1
PrideYearDraw <- matrix(as.numeric(PrideYearDraw), nrow = nrow(PrideYearMean),
                    ncol = ncol(PrideYearMean))
par(mfrow=c(1,2))
# Matrix plot
Colour <- colorRampPalette(c("blue", "white", "red"))(200)
corrplot(PrideYearMean, method = "color", col = Colour, type = "lower",
         diag = FALSE, p.mat = PrideYearDraw, tl.srt = 45)
# Chord diagram
chordDiagram(PrideYearMean, symmetric = TRUE,
             annotationTrack = c("name", "grid"),
             grid.col = "grey", col = PrideYearDrawCol, link.lty=0.1)
#--------------------
### Year level effect
#--------------------
# Build matrix of colours for chordDiagram. Change to 0.5 to make focus on strong interactions
YearDrawCol <- matrix(NA, nrow = nrow(YearMean),
                      ncol = ncol(YearMean))
YearDrawCol[which(YearMean > 0.4, arr.ind=TRUE)]<-"red"
YearDrawCol[which(YearMean < -0.4, arr.ind=TRUE)]<-"blue"

# Build matrix of "significance" for corrplot

YearDraw <- YearDrawCol
YearDraw[which(!is.na(YearDraw), arr.ind = TRUE)] <- 0
YearDraw[which(is.na(YearDraw), arr.ind = TRUE)] <- 1
YearDraw <- matrix(as.numeric(YearDraw), nrow = nrow(YearMean),
                   ncol = ncol(YearMean))

#do the following code 
par(mfrow=c(1,2))
# Matrix plot
Colour <- colorRampPalette(c("blue", "white", "red"))(200)
corrplot(YearMean, method = "color", col = Colour, type = "lower",
         diag = FALSE, p.mat = YearDraw, tl.srt = 45)
# Chord diagram
chordDiagram(YearMean, symmetric = TRUE,
             annotationTrack = c("name", "grid"),
             grid.col = "grey", col = YearDrawCol)

#=======================
### Traits
#=======================
mixingTr <-as.mcmc(modelMay2018traitsHighRes, parameters = "paramTr")
effectiveSize(mixingTr)
full <- ggs(mixingTr)
p <- ggs_caterpillar(full)
p + geom_point(size = 0.01, stroke = 0, shape = 16) 
#trait table
average <- apply(modelMay2018traits$results$estimation$paramTr, 1:2, mean)
CI.025 <- apply(modelMay2018traits$results$estimation$paramTr, 1:2, quantile,
                probs = 0.025)
CI.975 <- apply(modelMay2018traits$results$estimation$paramTr, 1:2, quantile,
                probs = 0.975)
CI <- cbind(as.vector(CI.025), as.vector(CI.975))
plot(0, 0, xlim = c(1, nrow(CI)), ylim = range(CI), type = "n",
     xlab = "", ylab = "", main = "paramX")
abline(h = 0, col = "grey")

arrows(x0 = 1:nrow(CI), x1 = 1:nrow(CI), y0 = CI[, 1], y1 = CI[, 2],
       code = 3, angle = 90, length = 0.05)
points(1:nrow(CI), average, pch = 15, cex = 1.5)

paramXCITrTable <- cbind(unlist(as.data.frame(average)),
                         unlist(as.data.frame(CI.025)),
                         unlist(as.data.frame(CI.975)))
colnames(paramXCITable) <- c("paramX", "lowerCI", "upperCI")
rownames(paramXCITable) <- paste(rep(colnames(average),
                                     each = nrow(average)), "_",
                                 rep(rownames(average),
                                     ncol(average)), sep="")
write.csv(paramXCITable, file='Pathtraits.csv')

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