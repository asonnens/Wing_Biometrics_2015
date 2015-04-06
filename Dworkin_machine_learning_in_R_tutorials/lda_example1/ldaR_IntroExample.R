
## ----libraries-----------------------------------------------------------
require(MASS) # contains lda() and qda() functions


## ----data----------------------------------------------------------------
wings <- read.csv("~/Dropbox/WingBiometrics2014/PredationWingsALL.csv", h=T)


## ------------------------------------------------------------------------
str(wings)


## ------------------------------------------------------------------------
table(wings$Sex)


## ----RHS-----------------------------------------------------------------
RHS <- paste("ProcCoord", 1:26, sep="", collapse= " + ")


## ----RHSoutput-----------------------------------------------------------
RHS


## ----LHS-----------------------------------------------------------------
DiscrimEqn <- as.formula(paste("Sex ~ ", RHS, sep="", collapse=""))


## ----LHSout--------------------------------------------------------------
DiscrimEqn


## ----lda-----------------------------------------------------------------
linDiscrim_0 <- lda(DiscrimEqn, data = wings, CV=F)


## ----lda_scaling---------------------------------------------------------
linDiscrim_0$scaling


## ----lda_1---------------------------------------------------------------
linDiscrim_1 <- lda(DiscrimEqn, data = wings, CV=T)


## ----lda_classification--------------------------------------------------
linDiscrim_1$class


## ----lda_posterior-------------------------------------------------------
head(linDiscrim_1$posterior, n=10)
tail(linDiscrim_1$posterior, n=10)


## ----lda_correct---------------------------------------------------------
linDiscrim_1$class == wings$Sex


## ----lda_prob_correct----------------------------------------------------
prob_correct <- sum(linDiscrim_1$class == wings$Sex)/length(wings$Sex)
prob_correct  


## ----vectors-------------------------------------------------------------
males <- wings[wings$Sex=="M",9:34]
females <- wings[wings$Sex=="F",9:34]
male.means <- colMeans(males)
female.means <- colMeans(females)


## ----VCV-----------------------------------------------------------------
pooled.cov <- cov(wings[, 9:34]) # This is S
S.inv <- solve(pooled.cov) # If you set the tolerance too low, all hell breaks loose, so for this example I am using a covariance matrix of full rank,


## ----disc_vec------------------------------------------------------------
a = S.inv %*% (male.means - female.means)

# rescaling it to compare to the original vector
a.prime = a*(34.488708/31.918649)
cbind(a.prime,linDiscrim_0$scaling)


