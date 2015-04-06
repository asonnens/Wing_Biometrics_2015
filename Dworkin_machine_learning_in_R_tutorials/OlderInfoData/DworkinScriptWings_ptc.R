# Analysis of wild-type VS patched using both the 15 landmarks and the full set of landmarks and semi-landmarks.
# Last Updated Dec15th


# Examining the wings for the patched (10514) vs wild-type (SAM) comparison. 

# Using lda and qda, how different are the results,
# Do we get similar results using the 15 LM (26 dimensions), then the high dimensional data with semi-landmarks (58 dimensions)


require(shapes) # library for Procrustres Superimposition
require(MASS)   # for lda and qda


setwd("/Users/ian/Dropbox/NSF-IIS Proposal 2011/wings")

wings_15LM <- read.csv("Wings_15LM_combined.csv", h=T)
str(wings_15LM)

wings_58Dim <- read.csv("Wings_Sam10514_HighDim.csv", h=T)
str(wings_58Dim)


# The 15 LM data has not yet undergone any registration.

range_for_plot <- range(wings_15LM[, 7:36])
plot(wings_15LM[,7], wings_15LM[,8], xlim=range_for_plot, ylim=range_for_plot, type="n", 
     xlab="", ylab="", main= "unregistered",
     frame.plot=F, axes=F, asp=1 )

for(i in seq(7,36, by=2)){
	points(wings_15LM[, i], wings_15LM[, i+1], col=c("red", "blue")[wings_15LM$Genotype1], cex=0.5 )
}


# register using Procrustes superimposition (GLS)
landmarks_15_unregistered <- t(as.matrix(wings_15LM[, 7:36])) # the program wants it as rows for LM, cols for observations
dim(landmarks_15_unregistered)

landmarks_15_unregistered <- array(landmarks_15_unregistered, dim = c(2, 15, ncol(landmarks_15_unregistered))) #dimensions, # of landmarks, sample number

landmarks_15_unregistered  <- aperm(landmarks_15_unregistered , c(2, 1, 3)) #gets the array into the right shape for Dryden's program

class(landmarks_15_unregistered)

landmarks_15_aligned <- procGPA(landmarks_15_unregistered)
plotshapes(landmarks_15_aligned$mshape)
plotshapes(landmarks_15_aligned$rotated)


landmarks_15_registered <-landmarks_15_aligned$rotated
dim(landmarks_15_registered)

landmarks_15_aligned$size
# Double check the plot....
range_for_plot.2 <- range(c(landmarks_15_registered[,1,],landmarks_15_registered[,2,] ))
plot(landmarks_15_registered[1,,], landmarks_15_registered[2,,],  
      xlim=range_for_plot.2, ylim=range_for_plot.2, type="n", xlab="", ylab="")

for(i in 1:15){
	points(landmarks_15_registered[i, 1,], landmarks_15_registered[i,2,],  cex=0.5 )
}




# Now I want to get it back into a matrix (I hate the arrays in R!!!)


crap1 <- t(as.matrix(landmarks_15_registered[,1,]))
colnames(crap1) <- c('LM_X1' ,  'LM_X2' ,  'LM_X3' ,  'LM_X4' ,  'LM_X5' ,  'LM_X6' ,  'LM_X7' ,  'LM_X8' ,  'LM_X9' ,  'LM_X10' ,  'LM_X11' ,  'LM_X12' ,  'LM_X13' ,  'LM_X14' ,  'LM_X15')

crap2 <- t(as.matrix(landmarks_15_registered[,2,]))
colnames(crap2) <- c('LM_Y1' ,  'LM_Y2' ,  'LM_Y3' ,  'LM_Y4' ,  'LM_Y5' ,  'LM_Y6' ,  'LM_Y7' ,  'LM_Y8' ,  'LM_Y9' ,  'LM_Y10' ,  'LM_Y11' ,  'LM_Y12' ,  'LM_Y13' ,  'LM_Y14' ,  'LM_Y15')

crap3 <- cbind(crap1, crap2)
colnames(crap3)


wings_15LM.new <- data.frame(wings_15LM, crap3, CSize = (landmarks_15_aligned$size) )
colnames(wings_15LM.new)



# Check one last time that everything looks fine!
range_for_plot.3 <- range(wings_15LM.new[, 37:66])
plot(x=wings_15LM.new[, 37], y=wings_15LM.new[, 37+15], 
      xlim=range(wings_15LM.new[, 37:51]), ylim=range(wings_15LM.new[, 52:66]), 
      type="n", xlab="", ylab="", main = "registered (Procrustes)", frame.plot=F, axes=F, asp=1)

for(i in 1:15){
	points(wings_15LM.new[, 36 + i], wings_15LM.new[, 36 + 15 + i], col=c("red", "blue")[wings_15LM.new$Genotype1], cex=0.5 )
}

#### Newbit as of Dec 15th 2011
# Address the finding by Dirk about the increase in variance for the mutant (Ptc or 10514)  as compared to the wild-type ( SAM)
# Use the total variance (which is the trace of the covariance matrix)
TotVarSam <- sum(diag(cov(wings_15LM.new[wings_15LM.new$Genotype1=="Sam", 37:66])))
TotVarPtc <- sum(diag(cov(wings_15LM.new[wings_15LM.new$Genotype1=="10514", 37:66])))

TotVarSam
TotVarPtc
# Nope, the finding for this mutant is in the opposite direction across all Landmarks (we could look at subset that were biologically most important)

#Confirm with eigenvalues (should be identical to what I computed above, just a computational double check)
EigenCovSam <- eigen(cov(wings_15LM.new[wings_15LM.new$Genotype1=="Sam", 37:66]))
EigenCovPtc <- eigen(cov(wings_15LM.new[wings_15LM.new$Genotype1=="10514", 37:66]))

sum(EigenCovSam$values) # Same estimate as TotVarSam
sum(EigenCovPtc$values)
#########



# Now do a plot for the 58D data set.

plot(x=wings_58Dim[, 7], y=wings_58Dim[, 8], 
      xlim=range(wings_58Dim[, seq(7, 101, by=2)]), ylim=range(wings_58Dim[, seq(8,102, by=2)]), 
      type="n", xlab="", ylab="", main = "registered (Procrustes with semi-landmarks)", 
      frame.plot=F, axes=F, asp=1)

for(i in seq(7,101, by=2)){
	points(wings_58Dim[,  i], wings_58Dim[, 1 + i], col=c("red", "blue")[wings_58Dim$Genotype1], cex=0.5 )
}


#### Newbit as of Dec 15th 2011
# Address the finding by Dirk about the increase in variance for the mutant (Ptc or 10514)  as compared to the wild-type ( SAM)
# Use the total variance (which is the trace of the covariance matrix).
# For all landmarks and semilandmarks
TotVarSam58D <- sum(diag(cov(wings_58Dim[wings_58Dim$Genotype1=="sam", 7:102])))
TotVarPtc58D <- sum(diag(cov(wings_58Dim[wings_58Dim$Genotype1=="10514", 7:102])))

TotVarSam58D
TotVarPtc58D
# Same as before.

#Confirm with eigenvalues (should be identical to what I computed above, just a computational double check)
EigenCovSam58D <- eigen(cov(wings_58Dim[wings_58Dim$Genotype1=="sam", 7:102]))
EigenCovPtc58D <- eigen(cov(wings_58Dim[wings_58Dim$Genotype1=="10514", 7:102]))

sum(EigenCovSam58D$values) # Same estimate as TotVarSam
sum(EigenCovPtc58D$values)
#########



# it is worth pointing out that centroid size has been removed, so it may need to be added back in for analysis.



# Let's perform LDA for the 15 LM, unregistered.

LDA_15_unregistered <- lda(Genotype1 ~ as.matrix(wings_15LM[, 7:36]), data=wings_15LM, CV=T)
LDA_15_unregistered$posterior
LDA_15_unregistered$class == wings_15LM$Genotype1

# 100% correct calls for the whole data set.
sum(LDA_15_unregistered$class == wings_15LM$Genotype1)/length(wings_15LM$Genotype1)
# still need to perform cross validation.....


# Principal components
PCA_15_unregistered <- prcomp(as.matrix(wings_15LM[, 7:36]))
biplot(PCA_15_unregistered)
plot(PCA_15_unregistered)


wings_15LM_PCA <- data.frame(wings_15LM, PCA_15_unregistered$x)
str(wings_15LM_PCA)

plot(PC1 ~ PC2, data=wings_15LM_PCA, col=c("red", "blue")[wings_15LM_PCA$Genotype1], asp=1)




# How about for the registered data. Column 67 is Csize (divided by 1000 so it is ~ same scale as shape)
LDA_15LM_registered <- lda(Genotype1 ~ as.matrix(wings_15LM.new[, 37:66]), data=wings_15LM.new, CV=T)
LDA_15LM_registered$posterior
LDA_15LM_registered$class == wings_15LM.new$Genotype1

# 100 correct calls for the whole data set.
sum(LDA_15LM_registered$class == wings_15LM.new$Genotype1)/length(wings_15LM.new$Genotype1) # it actually misses one!

# Principal components
PCA_15_registered <- prcomp(as.matrix(wings_15LM.new[, 37:66]))
biplot(PCA_15_registered)
plot(PCA_15_registered)


wings_15LM.new <- data.frame(wings_15LM.new, PCA_15_registered$x)
str(wings_15LM.new)

plot(PC1 ~ PC2, data=wings_15LM.new, col=c("red", "blue")[wings_15LM.new$Genotype1], asp=1)


# PC1 seems to be strongly influenced by size
plot(PC1 ~ CSize, data =wings_15LM.new, col=c("red", "blue")[wings_15LM.new$Genotype1])

# 