
## ----libraries-----------------------------------------------------------
require(MASS)     # contains lda() and qda() functions
require(sampling) # easier way to generate subsets of data
require(mda)      # flexible discriminant analysis and mixture discrimant analysis
require(class)    # k nearest neighbours (knn)
#require(rpart)    # Classification and regression trees (CART)
require(adabag)   # bagging()
require(ada)      # boosting function, ada()
require(nnet)     # neural net function, nnet()
require(e1071)    # svm()
require(randomForest) # randomForest() 


## ----data----------------------------------------------------------------
wings <- read.csv("~/Dropbox/WingBiometrics2014/PredationWingsALL.csv", h=T)


## ------------------------------------------------------------------------
str(wings)


## ----princomp------------------------------------------------------------
wings_pc <- prcomp(wings[,9:38])
summary(wings_pc)
head(wings_pc$x[,1:26]) # All of the variables we need
wings <- data.frame(wings, wings_pc$x[,1:26]) # appended the PCs


## ----table---------------------------------------------------------------
tab1 <- table(wings$Sex)
tab1


## ----strata--------------------------------------------------------------
wings_strata <- strata(wings, stratanames=c("Sex"), method="srswor",
    size=c( round(tab1[1]*(2/3)), round(tab1[2]*(2/3))))

names(wings_strata)    


## ----lda_1---------------------------------------------------------------
linDiscrim_1 <- lda(formula=Sex ~.,data = wings_training, tol = 1.0e-4, CV=FALSE)


## ----predict_table-------------------------------------------------------
lda_training_table <- table(actual = wings_training$Sex,
                          predicted = predict(linDiscrim_1, newdata=wings_training)$class)
lda_training_table                          


## ----predict_table_test--------------------------------------------------
lda_test_table <- table(actual = wings_test$Sex,
                          predicted = predict(linDiscrim_1, newdata=wings_test)$class)
lda_test_table                          


## ----qda_1---------------------------------------------------------------
qda_1 <- qda(formula=Sex ~.,data = wings_training, tol = 1.0e-4)
qda_training_table <- table(actual = wings_training$Sex,
                          predicted = predict(qda_1, newdata=wings_training)$class)
qda_training_table                          

qda_test_table <- table(actual = wings_test$Sex,
                          predicted = predict(qda_1, newdata=wings_test)$class)
qda_test_table                          


## ----fda_1---------------------------------------------------------------
fda_1 <- fda(formula=Sex ~.,data = wings_training)
fda_1
fda_1$confusion
fda_training_table <- table(actual = wings_training$Sex,
                          predicted = predict(fda_1, newdata=wings_training, type="class"))
fda_training_table                          

fda_test_table <- table(actual = wings_test$Sex,
                          predicted = predict(fda_1, newdata=wings_test, type="class"))
fda_test_table                          


## ----mda_1---------------------------------------------------------------
mda_1 <- mda(formula=Sex ~.,data = wings_training)
mda_1
mda_1$confusion
mda_training_table <- table(actual = wings_training$Sex,
                          predicted = predict(mda_1, newdata=wings_training))
mda_training_table                          

mda_test_table <- table(actual = wings_test$Sex,
                          predicted = predict(mda_1, newdata=wings_test))
mda_test_table  


## ----knn-----------------------------------------------------------------
wings_knn <- knn(train=wings_training[, -1], test=wings_test[,-1], cl=wings_training$Sex, k=1, prob=TRUE)
summary(wings_knn)

table_knn <- table(predicted=wings_knn, actual=wings_test$Sex)
table_knn


## ----knn_k5--------------------------------------------------------------
wings_knn_k5 <- knn(train=wings_training[, -1], test=wings_test[,-1], 
    cl=wings_training$Sex, k=5, prob=TRUE)
summary(wings_knn_k5)

table_knn_k5 <- table(predicted=wings_knn_k5, actual=wings_test$Sex)
table_knn_k5


## ----bag_1---------------------------------------------------------------
bag_1 <- bagging(formula=Sex ~.,data = wings_training)

bag_training_table <- table(actual = wings_training$Sex,
                          predicted = predict(bag_1, newdata=wings_training)$class)
bag_training_table                          

bag_test_table <- table(actual = wings_test$Sex,
                          predicted = predict(bag_1, newdata=wings_test)$class)
bag_test_table                          


## ----boo_1---------------------------------------------------------------
boo_1 <- ada(formula=Sex ~.,data = wings_training, loss="logistic")

boo_1

boo_training_table <- table(actual = wings_training$Sex,
                          predicted = predict(boo_1, newdata=wings_training))
boo_training_table                          

boo_test_table <- table(actual = wings_test$Sex,
                          predicted = predict(boo_1, newdata=wings_test))
boo_test_table                          


## ----nnet----------------------------------------------------------------
nnet_1 <- nnet(formula=Sex ~.,data = wings_training, size=10, decay=0)

nnet_training_table <- table(actual = wings_training$Sex,
                          predicted = predict(nnet_1, type="class"))
nnet_training_table                          

nnet_test_table <- table(actual = wings_test$Sex,
                          predicted = predict(nnet_1, newdata=wings_test, type="class"))
nnet_test_table                          



## ----svm-----------------------------------------------------------------
wings_svm <- svm(Sex ~., data = wings_training)

svm_training_table <- table(actual = wings_training$Sex,
                          predicted = predict(wings_svm, type="class"))
svm_training_table                          

svm_test_table <- table(actual = wings_test$Sex,
                          predicted = predict(wings_svm, newdata=wings_test, type="class"))
svm_test_table                          



## ----randomForest--------------------------------------------------------
wings_rf <- randomForest(Sex ~., data = wings_training)
wings_rf

rf_training_table <- table(actual = wings_training$Sex,
                          predicted = predict(wings_rf, type="class"))
rf_training_table                          

rf_test_table <- table(actual = wings_test$Sex,
                          predicted = predict(wings_rf, newdata=wings_test, type="class"))
rf_test_table                          



