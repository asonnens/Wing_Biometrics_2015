## Extending "classification" approaches with lda, qda, etc..

Author: Ian Dworkin. Created: April 14/2014

This continues on from the first tutoria into classifiying observations into classes for purposes of future predictions. In particular we will utilize so-called **supervised learning** approaches like _linear discriminants analysis_ (lda) _quadratic discriminants analysis_ (qda) and a number of other approaches like _support vector machines_ (svm) and non-parametric lda (mda). For all of this I am more or less using the approach from *R in a nutshell* by Joseph Adler (chapter 21). It is a useful `R` reference. 

In addition (and unlike the first tutorial), I will be subsetting the data into a training and testing partitions.

First we call the libraries we want to use (note I am not making any statements of which ones to use, just demonstrating performance on a very specific datasets).

```r
require(MASS)  # contains lda() and qda() functions
require(sampling)  # easier way to generate subsets of data
require(mda)  # flexible discriminant analysis and mixture discrimant analysis
require(class)  # k nearest neighbours (knn)
# require(rpart) # Classification and regression trees (CART)
require(adabag)  # bagging()
require(ada)  # boosting function, ada()
require(nnet)  # neural net function, nnet()
require(e1071)  # svm()
require(randomForest)  # randomForest() 
```


Now we read in data. Here we are using 15, 2-dimensional landmarks on *Drosophila* wings

```r
wings <- read.csv("~/Dropbox/WingBiometrics2014/PredationWingsALL.csv", h = T)
```


### Notes on this data.
1. The landmarks were acquired "manually" from images using the ImageJ plug-in from Chris Klingenberg. This is not the 58 dimensional data (landmarks and semi-landmarks) from the semi-automated [wing machine](http://bio.fsu.edu/~dhoule/Software/) software from the lab of Dr. David Houle that we normally use these days (2008-present).
2. These images are from an analysis of phenotypic selection of juvenile mantids on *Drosophila melanogaster*.
3. For the purposes of this overview, we only care about the shape data and sex of flies. We are going to ignore other variables (like survivorship)
4. Side refers to Left, Right or unknown (from wings that fell off when the flies were being eaten by predator). We should take care of this to avoid pseudo-replication, but for the moment we will ignore it for the moment, but will come back to it. Do not ignore it for a real analysis.
5. Since we are starting by examining sex, we should also include size as variable, or at least use size corrected shape variables. Again, for the purposes of getting this all started, we are ignoring it.
6.Importanly, since we are using procrustes superimposed data that has been registered (centered) and scaled to centroid size, we only have 26 dimensions, despite having thirty variables. Thus we need to extract the first 26 principal components from the shape data (which contains **ALL** the information)

Let's look at the structure of the data object.

```r
str(wings)
```

```
## 'data.frame':	988 obs. of  38 variables:
##  $ Unique_ID        : Factor w/ 971 levels "G4_F_C_1_U","G4_F_C_10_U",..: 81 82 83 84 85 86 87 88 89 90 ...
##  $ Population       : Factor w/ 3 levels "G4","G9","S4": 2 2 2 2 2 2 2 2 2 2 ...
##  $ Sex              : Factor w/ 2 levels "F","M": 1 1 1 1 1 1 1 1 1 1 ...
##  $ Treatment        : Factor w/ 3 levels "C","D","S": 1 1 1 1 1 1 1 1 1 1 ...
##  $ Individual       : int  1 1 2 2 3 3 4 4 5 5 ...
##  $ Side             : Factor w/ 3 levels "L","R","U": 1 2 1 2 1 2 1 2 1 2 ...
##  $ Centroid.Size    : num  2.28 2.28 2.6 2.59 2.78 ...
##  $ Log.Centroid.Size: num  0.823 0.823 0.954 0.952 1.021 ...
##  $ ProcCoord1       : num  0.25 0.249 0.249 0.248 0.25 ...
##  $ ProcCoord2       : num  0.0697 0.0693 0.0688 0.0691 0.0687 ...
##  $ ProcCoord3       : num  0.239 0.241 0.239 0.241 0.241 ...
##  $ ProcCoord4       : num  0.0227 0.0239 0.0225 0.0227 0.0231 ...
##  $ ProcCoord5       : num  0.267 0.266 0.269 0.268 0.267 ...
##  $ ProcCoord6       : num  -0.00989 -0.00954 -0.01114 -0.01197 -0.00808 ...
##  $ ProcCoord7       : num  0.206 0.206 0.207 0.208 0.208 ...
##  $ ProcCoord8       : num  -0.0543 -0.0529 -0.0541 -0.0534 -0.0551 ...
##  $ ProcCoord9       : num  0.198 0.198 0.198 0.199 0.197 ...
##  $ ProcCoord10      : num  -0.031 -0.0317 -0.0316 -0.0319 -0.0325 ...
##  $ ProcCoord11      : num  0.173 0.171 0.168 0.168 0.17 ...
##  $ ProcCoord12      : num  0.0234 0.0241 0.0252 0.0248 0.0251 ...
##  $ ProcCoord13      : num  0.138 0.138 0.141 0.141 0.135 ...
##  $ ProcCoord14      : num  0.0877 0.0887 0.0861 0.0872 0.0861 ...
##  $ ProcCoord15      : num  0.0522 0.0566 0.062 0.0652 0.0664 ...
##  $ ProcCoord16      : num  0.0214 0.0215 0.0208 0.0214 0.0199 ...
##  $ ProcCoord17      : num  0.0533 0.0552 0.0647 0.0658 0.0689 ...
##  $ ProcCoord18      : num  -0.00409 -0.00452 -0.00366 -0.00254 -0.00355 ...
##  $ ProcCoord19      : num  -0.0809 -0.084 -0.0891 -0.0915 -0.086 ...
##  $ ProcCoord20      : num  -0.0304 -0.0338 -0.0277 -0.0269 -0.0302 ...
##  $ ProcCoord21      : num  -0.075 -0.0745 -0.0862 -0.0886 -0.0877 ...
##  $ ProcCoord22      : num  -0.0957 -0.0976 -0.0951 -0.095 -0.0994 ...
##  $ ProcCoord23      : num  -0.309 -0.313 -0.319 -0.317 -0.329 ...
##  $ ProcCoord24      : num  0.169 0.168 0.161 0.16 0.157 ...
##  $ ProcCoord25      : num  -0.482 -0.478 -0.477 -0.476 -0.473 ...
##  $ ProcCoord26      : num  0.0515 0.0535 0.0499 0.0483 0.0493 ...
##  $ ProcCoord27      : num  -0.458 -0.456 -0.454 -0.456 -0.453 ...
##  $ ProcCoord28      : num  -0.0352 -0.0333 -0.0323 -0.031 -0.0243 ...
##  $ ProcCoord29      : num  -0.172 -0.174 -0.173 -0.174 -0.174 ...
##  $ ProcCoord30      : num  -0.184 -0.185 -0.179 -0.181 -0.176 ...
```


### Extracting the principal components 
We will use `prcomp()` as this uses singular value decomposition which is a bit more robust than eigendecomposition. 


```r
wings_pc <- prcomp(wings[, 9:38])
summary(wings_pc)
```

```
## Importance of components:
##                           PC1     PC2     PC3     PC4     PC5     PC6
## Standard deviation     0.0108 0.00842 0.00752 0.00648 0.00582 0.00507
## Proportion of Variance 0.2780 0.16731 0.13353 0.09910 0.07996 0.06077
## Cumulative Proportion  0.2780 0.44531 0.57885 0.67795 0.75791 0.81868
##                           PC7     PC8     PC9    PC10    PC11    PC12
## Standard deviation     0.0035 0.00338 0.00288 0.00267 0.00247 0.00225
## Proportion of Variance 0.0289 0.02704 0.01963 0.01688 0.01438 0.01200
## Cumulative Proportion  0.8476 0.87460 0.89423 0.91111 0.92548 0.93748
##                           PC13    PC14    PC15    PC16    PC17    PC18
## Standard deviation     0.00204 0.00199 0.00178 0.00172 0.00153 0.00141
## Proportion of Variance 0.00986 0.00940 0.00748 0.00700 0.00556 0.00467
## Cumulative Proportion  0.94734 0.95673 0.96421 0.97122 0.97678 0.98145
##                           PC19    PC20    PC21    PC22     PC23    PC24
## Standard deviation     0.00126 0.00114 0.00111 0.00102 0.000907 0.00088
## Proportion of Variance 0.00376 0.00306 0.00291 0.00248 0.001940 0.00183
## Cumulative Proportion  0.98521 0.98828 0.99119 0.99366 0.995610 0.99744
##                            PC25    PC26     PC27     PC28     PC29
## Standard deviation     0.000814 0.00065 5.39e-09 2.98e-10 2.81e-10
## Proportion of Variance 0.001560 0.00100 0.00e+00 0.00e+00 0.00e+00
## Cumulative Proportion  0.999000 1.00000 1.00e+00 1.00e+00 1.00e+00
##                            PC30
## Standard deviation     2.73e-10
## Proportion of Variance 0.00e+00
## Cumulative Proportion  1.00e+00
```

```r
head(wings_pc$x[, 1:26])  # All of the variables we need
```

```
##            PC1        PC2        PC3       PC4        PC5       PC6
## [1,]  0.011070 -0.0074440  0.0134332 -0.013453  0.0012779  0.008113
## [2,]  0.007639 -0.0067695  0.0124251 -0.006600  0.0020827  0.006221
## [3,] -0.003166 -0.0004478  0.0001809 -0.005185 -0.0007805 -0.002859
## [4,] -0.002419 -0.0022569 -0.0046148 -0.002972 -0.0016012 -0.003358
## [5,] -0.012588  0.0056047  0.0008588  0.005268 -0.0018362 -0.001347
## [6,] -0.015720  0.0051138 -0.0053862  0.006833 -0.0059565 -0.004863
##            PC7        PC8       PC9       PC10      PC11       PC12
## [1,]  0.002699  4.014e-03 -0.001773  0.0078021  0.003272  9.673e-04
## [2,]  0.002238  4.119e-03  0.001531  0.0046466  0.001639 -4.841e-05
## [3,]  0.001188  2.358e-03 -0.004367  0.0017926 -0.003163  1.408e-04
## [4,]  0.002299 -3.694e-05 -0.004084  0.0022248 -0.003433  1.116e-03
## [5,] -0.003676 -3.260e-03 -0.003028  0.0025903 -0.001859  2.942e-04
## [6,] -0.004221 -1.911e-03  0.001883 -0.0008215 -0.002413  7.454e-05
##            PC13       PC14       PC15       PC16       PC17       PC18
## [1,]  0.0005957 -0.0016598 -0.0005507  5.200e-04  0.0004744  1.117e-03
## [2,]  0.0024150 -0.0009685 -0.0008530  8.252e-04 -0.0012561 -5.657e-04
## [3,]  0.0008143  0.0003967  0.0011783 -1.489e-03  0.0022350  8.649e-04
## [4,]  0.0018286  0.0017987 -0.0010462  6.860e-04  0.0033499  1.965e-03
## [5,] -0.0012558 -0.0014640  0.0014623  7.427e-05 -0.0013443  1.852e-04
## [6,] -0.0019348 -0.0031586  0.0005488  1.059e-03 -0.0018435 -5.414e-06
##            PC19       PC20       PC21       PC22       PC23       PC24
## [1,] -0.0004079 -0.0008209 -0.0002909  0.0006315 -5.770e-04  0.0007675
## [2,]  0.0002205 -0.0029924  0.0002370 -0.0014635  8.130e-04  0.0011403
## [3,]  0.0001877 -0.0005082 -0.0004081  0.0006695 -2.964e-04 -0.0012019
## [4,] -0.0002751 -0.0017065 -0.0001227  0.0001086  1.145e-04 -0.0005314
## [5,] -0.0006763  0.0006450 -0.0016224 -0.0017930 -8.994e-05  0.0006803
## [6,]  0.0006898 -0.0006194 -0.0014567 -0.0001087  9.250e-04  0.0006271
##            PC25       PC26
## [1,]  7.969e-04 -0.0003434
## [2,]  9.106e-04 -0.0011449
## [3,] -1.647e-04 -0.0001115
## [4,] -3.311e-04 -0.0002310
## [5,] -5.659e-05  0.0001270
## [6,] -8.437e-04  0.0001708
```

```r
wings <- data.frame(wings, wings_pc$x[, 1:26])  # appended the PCs
```



### Back to the important stuff
In particular we are interested in how many observations we have for each sex.

```r
tab1 <- table(wings$Sex)
tab1
```

```
## 
##   F   M 
## 583 405
```

So we have 583 females and 405 males in this data set. For subsetting we want to use ~2/3 for the training set and ~1/3 for test set, so 389 for females and 270 for males. We will do this using the `strata()` in the sampling library.


```r
wings_strata <- strata(wings, stratanames = c("Sex"), method = "srswor", size = c(round(tab1[1] * 
    (2/3)), round(tab1[2] * (2/3))))

names(wings_strata)
```

```
## [1] "Sex"     "ID_unit" "Prob"    "Stratum"
```


`ID_unit` provides the row number. We could have also done this with the `sample(,replace=FALSE)` by specifying probabilities, and accounting for both sexes, but this should be easier. To create the training and test sets we will utilize the index for matching observations using `ID_unit`.
 
```{R subsetting}
wings_2 <- wings[, c(3, 39:64)] # Just extracting the variables needed for this tutorial
wings_training <- wings_2[rownames(wings_2) %in% wings_strata$ID_unit, ]
wings_test <- wings_2[!rownames(wings_2) %in% wings_strata$ID_unit, ] # notice the '!'
nrow(wings_training)
nrow(wings_test)
nrow(wings)
``` 

## linear discriminant analysis.

I am running the lda a bit different from how I performed it in the previous tutorial. Instead of writing out the left and right hand sides of the model, I am going to do it all internal to the call to `lda()`. 

We start by running this on our training set as follows. Note we had to set the tolerance pretty low as the matrix was difficult to invert.

```r
linDiscrim_1 <- lda(formula = Sex ~ ., data = wings_training, tol = 1e-04, CV = FALSE)
```


We can ask which ones did lda classify correctly, as compared to the true sex. We can check this using a table, which will make this easier to see and use. Note `CV=FALSE` needs to be set in `lda()` for this approach to work. Not sure why they set it up that way. Does not really matter since we do not want to use 1-Fold cross validation.


```r
lda_training_table <- table(actual = wings_training$Sex, predicted = predict(linDiscrim_1, 
    newdata = wings_training)$class)
lda_training_table
```

```
##       predicted
## actual   F   M
##      F 368  21
##      M  17 253
```

Which results in 94.2337 percent correct classification.

How about for the training set? We do it almost identically in terms of generating the table.


```r
lda_test_table <- table(actual = wings_test$Sex, predicted = predict(linDiscrim_1, 
    newdata = wings_test)$class)
lda_test_table
```

```
##       predicted
## actual   F   M
##      F 187   7
##      M   7 128
```

Which results in 95.7447 percent correct classification. So it seems like our linear discriminant classifier did pretty well to start with.

## quadratic discriminant analysis

Now we move on to quadratic discriminant analysis, with the `qda()` in the MASS library. This has almost the same function call. Indeed many of the functions we will use for classification do.


```r
qda_1 <- qda(formula = Sex ~ ., data = wings_training, tol = 1e-04)
qda_training_table <- table(actual = wings_training$Sex, predicted = predict(qda_1, 
    newdata = wings_training)$class)
qda_training_table
```

```
##       predicted
## actual   F   M
##      F 377  12
##      M  25 245
```

```r

qda_test_table <- table(actual = wings_test$Sex, predicted = predict(qda_1, 
    newdata = wings_test)$class)
qda_test_table
```

```
##       predicted
## actual   F   M
##      F 186   8
##      M  15 120
```

Which results in 94.3854 percent correct classification for the training data and 93.0091 percent correct classification for the test data.

## Flexible discriminant analysis (FDA)

This uses a non-parametric regression instead of the 'linear' regression type approach in lda. We will use the `fda()` function in the mda library. Otherwise the syntax is similar.


```r
fda_1 <- fda(formula = Sex ~ ., data = wings_training)
fda_1
```

```
## Call:
## fda(formula = Sex ~ ., data = wings_training)
## 
## Dimension: 1 
## 
## Percent Between-Group Variance Explained:
##  v1 
## 100 
## 
## Degrees of Freedom (per dimension): 27 
## 
## Training Misclassification Error: 0.05766 ( N = 659 )
```

```r
fda_1$confusion
```

```
##          true
## predicted   F   M
##         F 368  17
##         M  21 253
## attr(,"error")
## [1] 0.05766
```

```r
fda_training_table <- table(actual = wings_training$Sex, predicted = predict(fda_1, 
    newdata = wings_training, type = "class"))
fda_training_table
```

```
##       predicted
## actual   F   M
##      F 368  21
##      M  17 253
```

```r

fda_test_table <- table(actual = wings_test$Sex, predicted = predict(fda_1, 
    newdata = wings_test, type = "class"))
fda_test_table
```

```
##       predicted
## actual   F   M
##      F 187   7
##      M   7 128
```

Which results in 94.2337 percent correct classification for the training data and 95.7447 percent correct classification for the test data.

## mixture discrimant analysis
`mda()` is used.



```r
mda_1 <- mda(formula = Sex ~ ., data = wings_training)
mda_1
```

```
## Call:
## mda(formula = Sex ~ ., data = wings_training)
## 
## Dimension: 5 
## 
## Percent Between-Group Variance Explained:
##     v1     v2     v3     v4     v5 
##  61.72  82.47  93.60  98.69 100.00 
## 
## Degrees of Freedom (per dimension): 27 
## 
## Training Misclassification Error: 0.05463 ( N = 659 )
## 
## Deviance: 176.1
```

```r
mda_1$confusion
```

```
##          true
## predicted   F   M
##         F 371  18
##         M  18 252
## attr(,"error")
## [1] 0.05463
```

```r
mda_training_table <- table(actual = wings_training$Sex, predicted = predict(mda_1, 
    newdata = wings_training))
mda_training_table
```

```
##       predicted
## actual   F   M
##      F 371  18
##      M  18 252
```

```r

mda_test_table <- table(actual = wings_test$Sex, predicted = predict(mda_1, 
    newdata = wings_test))
mda_test_table
```

```
##       predicted
## actual   F   M
##      F 187   7
##      M   7 128
```


Which results in 94.5372 percent correct classification for the training data and 95.7447 percent correct classification for the test data.

## k nearest neighbours

Using the `knn()`. The `k=1` specifies the number of nearest neighbours to use (defaults to one). Important to tune this. Looks like you need to make sure that only the numeric parts enter.


```r
wings_knn <- knn(train = wings_training[, -1], test = wings_test[, -1], cl = wings_training$Sex, 
    k = 1, prob = TRUE)
summary(wings_knn)
```

```
##   F   M 
## 197 132
```

```r

table_knn <- table(predicted = wings_knn, actual = wings_test$Sex)
table_knn
```

```
##          actual
## predicted   F   M
##         F 182  15
##         M  12 120
```

Which means we get 91.7933, which is much lower than with other methods.

### Trying a different k

Let us try a different value of 'k' such as k=5. I am not sure (yet) whether there is an automated way to find the best k, but we could always write a for loop if need be.


```r
wings_knn_k5 <- knn(train = wings_training[, -1], test = wings_test[, -1], cl = wings_training$Sex, 
    k = 5, prob = TRUE)
summary(wings_knn_k5)
```

```
##   F   M 
## 203 126
```

```r

table_knn_k5 <- table(predicted = wings_knn_k5, actual = wings_test$Sex)
table_knn_k5
```

```
##          actual
## predicted   F   M
##         F 183  20
##         M  11 115
```

Which means we get 90.5775 percent. Still not great.

## Bagging
Using the `bagging()` in the `adabag` library. This is a bit slow.


```r
bag_1 <- bagging(formula = Sex ~ ., data = wings_training)

bag_training_table <- table(actual = wings_training$Sex, predicted = predict(bag_1, 
    newdata = wings_training)$class)
bag_training_table
```

```
##       predicted
## actual   F   M
##      F 361  28
##      M  22 248
```

```r

bag_test_table <- table(actual = wings_test$Sex, predicted = predict(bag_1, 
    newdata = wings_test)$class)
bag_test_table
```

```
##       predicted
## actual   F   M
##      F 171  23
##      M  22 113
```


Which results in 92.4127 percent correct classification for the training data and 86.3222 percent correct classification for the test data.

## boosting
Using the `ada()` in the `ada` library

```r
boo_1 <- ada(formula = Sex ~ ., data = wings_training, loss = "logistic")

boo_1
```

```
## Call:
## ada(Sex ~ ., data = wings_training, loss = "logistic")
## 
## Loss: logistic Method: discrete   Iteration: 50 
## 
## Final Confusion Matrix for Data:
##           Final Prediction
## True value   F   M
##          F 387   2
##          M   4 266
## 
## Train Error: 0.009 
## 
## Out-Of-Bag Error:  0.047  iteration= 49 
## 
## Additional Estimates of number of iterations:
## 
## train.err1 train.kap1 
##         50         50
```

```r

boo_training_table <- table(actual = wings_training$Sex, predicted = predict(boo_1, 
    newdata = wings_training))
boo_training_table
```

```
##       predicted
## actual   F   M
##      F 387   2
##      M   4 266
```

```r

boo_test_table <- table(actual = wings_test$Sex, predicted = predict(boo_1, 
    newdata = wings_test))
boo_test_table
```

```
##       predicted
## actual   F   M
##      F 178  16
##      M  13 122
```


Which results in 99.0895 percent correct classification for the training data and 91.1854 percent correct classification for the test data.

## Neural networks

`nnet()` in `nnet` library.


```r
nnet_1 <- nnet(formula = Sex ~ ., data = wings_training, size = 10, decay = 0)
```

```
## # weights:  281
## initial  value 517.724579 
## iter  10 value 219.371612
## iter  20 value 157.244571
## iter  30 value 109.183875
## iter  40 value 97.961447
## iter  50 value 91.664049
## iter  60 value 88.864830
## iter  70 value 87.449315
## iter  80 value 84.159023
## iter  90 value 78.857635
## iter 100 value 74.015248
## final  value 74.015248 
## stopped after 100 iterations
```

```r

nnet_training_table <- table(actual = wings_training$Sex, predicted = predict(nnet_1, 
    type = "class"))
nnet_training_table
```

```
##       predicted
## actual   F   M
##      F 372  17
##      M  14 256
```

```r

nnet_test_table <- table(actual = wings_test$Sex, predicted = predict(nnet_1, 
    newdata = wings_test, type = "class"))
nnet_test_table
```

```
##       predicted
## actual   F   M
##      F 183  11
##      M  12 123
```

Which results in 95.2959 percent correct classification for the training data and 93.0091 percent correct classification for the test data. It is important to evaluate some of the tuning parameters for this approach.


## Support vector machines (SVMs)

Using the `svm()` in the `e1071` library.


```r
wings_svm <- svm(Sex ~ ., data = wings_training)

svm_training_table <- table(actual = wings_training$Sex, predicted = predict(wings_svm, 
    type = "class"))
svm_training_table
```

```
##       predicted
## actual   F   M
##      F 383   6
##      M   6 264
```

```r

svm_test_table <- table(actual = wings_test$Sex, predicted = predict(wings_svm, 
    newdata = wings_test, type = "class"))
svm_test_table
```

```
##       predicted
## actual   F   M
##      F 183  11
##      M   7 128
```


Which results in 98.1791 percent correct classification for the training data and 94.5289 percent correct classification for the test data. It is important to evaluate some of the tuning parameters for this approach.


## random forests

Using the `randomForest()` in the `randomForest` library.


```r
wings_rf <- randomForest(Sex ~ ., data = wings_training)
wings_rf
```

```
## 
## Call:
##  randomForest(formula = Sex ~ ., data = wings_training) 
##                Type of random forest: classification
##                      Number of trees: 500
## No. of variables tried at each split: 5
## 
##         OOB estimate of  error rate: 11.23%
## Confusion matrix:
##     F   M class.error
## F 350  39      0.1003
## M  35 235      0.1296
```

```r

rf_training_table <- table(actual = wings_training$Sex, predicted = predict(wings_rf, 
    type = "class"))
rf_training_table
```

```
##       predicted
## actual   F   M
##      F 350  39
##      M  35 235
```

```r

rf_test_table <- table(actual = wings_test$Sex, predicted = predict(wings_rf, 
    newdata = wings_test, type = "class"))
rf_test_table
```

```
##       predicted
## actual   F   M
##      F 180  14
##      M  16 119
```


Which results in 88.7709 percent correct classification for the training data and 90.8815 percent correct classification for the test data. It is important to evaluate some of the tuning parameters for this approach.

