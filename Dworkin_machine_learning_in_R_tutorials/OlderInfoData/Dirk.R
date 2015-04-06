#!/usr/bin/env Rscript
# Analysis of wild-type VS patched using both the 15 landmarks and the full set of landmarks and semi-landmarks.
# Last Updated Dec15th


#Process inputarguments to get input file names
args <- commandArgs(trailingOnly = TRUE)
cat("Input Arguments =", args, "\n")

if(length(args) == 0) {
	inputfile <- "cv2_6342_landmarks/cv2_15landmarks.csv"
	cat("No input argument given for csv file, using default:",inputfile,"\n")
}
if(length(args) > 0) { 
	inputfile <- args[1]
}

wings_15LM <- read.csv(inputfile, head=TRUE, strip.white=TRUE)
str(wings_15LM)
length(wings_15LM)


if(length(args) > 1) {
	start_col <- as.integer(args[2])
} else {
	start_col <- 7
}
if(length(args) > 2) {
	end_col <- as.integer(args[3])
} else {
	end_col <- length(wings_15LM)
}

cat("start_col=",start_col,"end_col=",end_col,"\n")
range <- start_col:end_col

# Examining the wings for the cv-2 (6342) mutant vs wild-type (SAM) comparison. 

# Using lda and qda, how different are the results,
# Do we get similar results using the 15 LM (26 dimensions), then the high dimensional data with semi-landmarks (58 dimensions)

require(MASS)   # for lda and qda


LDA_15_registered <- lda(G1 ~ as.matrix(wings_15LM[, range]), data=wings_15LM, CV=T)
#LDA_15_registered$posterior
#LDA_15_registered$class == wings_15LM$G1

# 91% correct calls for the whole data set.
CorrectCalls <- sum(LDA_15_registered$class == wings_15LM$G1)/length(wings_15LM$G1)

#write.csv(CorrectCalls, "../Routput")
cat("RESULTS:",CorrectCalls,"\n")
cat("Program Completed Succssfully\n")

# still need to perform cross validation.....

# Let's look at the distribution of the linear combinations (disrcriminant scores)
# LDA_15_registered.1 <- lda(G1 ~ as.matrix(wings_15LM[, 7:36]), data=wings_15LM, CV=F, method="mle") # other methods are moment and "t"
# DiscrimVect <- as.numeric(LDA_15_registered.1$scaling) # vector of 'loadings'

# wings_15LM$DiscrimValues <- as.matrix(wings_15LM[, 7:36]) %*% DiscrimVect  # project the landmark data onto the discriminant vector
# par(mfrow=c(2,1))
# boxplot(wings_15LM$DiscrimValues ~ wings_15LM$G1)
# plot(density(wings_15LM$DiscrimValues[wings_15LM$G1=="Samar"]), main="kernel density estimate for discriminant function", 
  # ylim=c(0,0.5), xlim=c(-40, -30), xlab="arbitrary discriminant score")
# lines(density(wings_15LM$DiscrimValues[wings_15LM$G1=="6342"], bw=0.3435), col="red")
# # LiFeng is sampling along index using 90% for estimating and 10% for testing. I need to Implement this.


# # Principal components
# PCA_15_registered <- prcomp(as.matrix(wings_15LM[, 7:36]))
# biplot(PCA_15_registered)
# plot(PCA_15_registered)


# wings_15LM_PCA <- data.frame(wings_15LM, PCA_15_registered$x)
# str(wings_15LM_PCA)

# plot(PC1 ~ PC2, data=wings_15LM_PCA, col=c("red", "blue")[wings_15LM_PCA$G1], asp=1)




# # PC1 and PC2 seems to be somewhat influenced by size
# plot(PC2 ~ Log_Centroid_Size, data =wings_15LM_PCA, col=c("red", "blue")[wings_15LM_PCA$G1])
