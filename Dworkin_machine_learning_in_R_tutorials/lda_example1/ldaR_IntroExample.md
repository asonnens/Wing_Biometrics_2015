# A brief tutorial into 'classification' analyses starting with linear discriminant analysis

Written by Ian Dworkin -April 14 2014.

This is a brief tutorial into classifiying observations into classes for purposes of future predictions. In particular we will utilize so-called **supervised learning** approaches like _linear discriminants analysis_ (lda) get and then quadratic discriminants analysis (qda) and a number of other approaches like svm in a future tutorial.

First we call the libraries we want to use

```r
require(MASS)  # contains lda() and qda() functions
```


Now we read in data. Here we are using 15, 2-dimensional landmarks on *Drosophila* wings

```r
wings <- read.csv("~/Dropbox/WingBiometrics2014/PredationWingsALL.csv", h = T)
```


### Some notes on this data.
1. The landmarks were acquired "manually" from images using the ImageJ plug-in from Chris Klingenberg. This is not the 58 dimensional data (landmarks and semi-landmarks) from the semi-automated [wing machine](http://bio.fsu.edu/~dhoule/Software/) software from the lab of Dr. David Houle that we normally use these days (2008-present).
2. These images are from an analysis of phenotypic selection of juvenile mantids on *Drosophila melanogaster*.
3. For the purposes of this overview, we only care about the shape data and sex of flies. We are going to ignore other variables (like survivorship)
4. Side refers to Left, Right or unknown (from wings that fell off when the flies were being eaten by predator). We should take care of this to avoid pseudo-replication, but for the moment we will ignore it for the moment, but will come back to it. Do not ignore it for a real analysis.
5. Since we are starting by examining sex, we should also include size as variable, or at least use size corrected shape variables. Again, for the purposes of getting this all started, we are ignoring it.


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


In particular we are interested in how many observations we have for each sex.

```r
table(wings$Sex)
```

```
## 
##   F   M 
## 583 405
```

So we have 583 females and 405 males in this data set.

## linear discriminant analysis.

Normally we would start by subsetting the data into a 'training' set and a 'validation set'. We will do that below, but first let's go over some of the basics.

For convenience and clarity I am seperating out the left and right hand sides of the equation. There are other ways of doing this as well. Keep in mind that for linear discrimant analysis the left and right hand side of the equation are 'reversed' relative to a general linear model. THat is the categorical variable is now on the LHS.

Putting together right hand side of the model (to make life easier for inputting)

```r
RHS <- paste("ProcCoord", 1:26, sep = "", collapse = " + ")
```

It is worth taking a look at what has been output

```r
RHS
```

```
## [1] "ProcCoord1 + ProcCoord2 + ProcCoord3 + ProcCoord4 + ProcCoord5 + ProcCoord6 + ProcCoord7 + ProcCoord8 + ProcCoord9 + ProcCoord10 + ProcCoord11 + ProcCoord12 + ProcCoord13 + ProcCoord14 + ProcCoord15 + ProcCoord16 + ProcCoord17 + ProcCoord18 + ProcCoord19 + ProcCoord20 + ProcCoord21 + ProcCoord22 + ProcCoord23 + ProcCoord24 + ProcCoord25 + ProcCoord26"
```


Now we put this all together into a formula

```r
DiscrimEqn <- as.formula(paste("Sex ~ ", RHS, sep = "", collapse = ""))
```

Once again let's take a look at it


```r
DiscrimEqn
```

```
## Sex ~ ProcCoord1 + ProcCoord2 + ProcCoord3 + ProcCoord4 + ProcCoord5 + 
##     ProcCoord6 + ProcCoord7 + ProcCoord8 + ProcCoord9 + ProcCoord10 + 
##     ProcCoord11 + ProcCoord12 + ProcCoord13 + ProcCoord14 + ProcCoord15 + 
##     ProcCoord16 + ProcCoord17 + ProcCoord18 + ProcCoord19 + ProcCoord20 + 
##     ProcCoord21 + ProcCoord22 + ProcCoord23 + ProcCoord24 + ProcCoord25 + 
##     ProcCoord26
```


#### Now we perform the lda
 
For the same model, we will need to run it twice. For some reason when we set the cross-validation flag to CV=T, it does not provide the vector of discriminations 

```r
linDiscrim_0 <- lda(DiscrimEqn, data = wings, CV = F)
```


Since the supervised variable we are interested in (sex) has only two levels, the lda will only produce one vector that discriminantes (we are only discriminating in a single dimension).


```r
linDiscrim_0$scaling
```

```
##                  LD1
## ProcCoord1   -34.489
## ProcCoord2   177.933
## ProcCoord3    16.768
## ProcCoord4  -179.750
## ProcCoord5    -6.883
## ProcCoord6   109.968
## ProcCoord7  -112.176
## ProcCoord8    87.017
## ProcCoord9    45.013
## ProcCoord10 -105.772
## ProcCoord11  -43.168
## ProcCoord12   86.987
## ProcCoord13  -80.005
## ProcCoord14   35.794
## ProcCoord15   12.771
## ProcCoord16   51.418
## ProcCoord17  102.529
## ProcCoord18  -40.635
## ProcCoord19    6.644
## ProcCoord20 -116.080
## ProcCoord21  -81.341
## ProcCoord22  247.453
## ProcCoord23   -8.985
## ProcCoord24  267.234
## ProcCoord25   46.402
## ProcCoord26   12.366
```


We run the same model, only this time with CV=T.


```r
linDiscrim_1 <- lda(DiscrimEqn, data = wings, CV = T)
```


Note that we set flag for cross validation to TRUE. This is simply leave-one-out cross-validation, which evaluates the predictive performance, but it is not as good as a true training and test set like we will perform later.

We can look at the calls for the classifications (not true sex of fly, but what lda is calling it).

```r
linDiscrim_1$class
```

```
##   [1] F F F F F F F F F F F F F F F F F F F F F F F F M F F F F F F F F F F
##  [36] F F F F M F F F F F F F F F F M F F F F F M F F F F F F M F F M F F F
##  [71] F F F M F F F M M M F F F F F F F F F F F F F F F F F F F F F F F F F
## [106] F F F F M F F F F F F F F F F F F F F F F F F F M F F F F F F F F F F
## [141] F F F F F F M F F M F F F F F F F F F F F F F M M F F F F F F F F F F
## [176] F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F F
## [211] F F F F F F F F F F F F F F F F F F F F F F F F F F M F F F F F F F F
## [246] M M F F M F F F F F F F F F F F F F F F F F F F F F F F M M M M M M M
## [281] M M M M M M M M M M M M M M M M M M M M M M M M M M M M M M M M M M M
## [316] M M M M M M M M M M M M M M M M M M M M M M M M M M M M M M M M M M F
## [351] F M M M M M M M M M M M F M M F F M M M M M M M M M F M M F M M M M M
## [386] M M M F M M M M M M M M M F M M M M F M M M M F M M M M F M M M M M M
## [421] M M M M M M M M M M M M M M M M M M M M M M M M M M M M M M M M M M M
## [456] M M M M M M M F F F M M M M F F F F F M F F F F F F F F F F F F F F F
## [491] F F F F F F F F F F F F F F F F F M F F F F F F F F F F F F F F F F F
## [526] F F F F F F F F F F F F F F F F F F F F M M F F F F F F F F F F F F F
## [561] F F F F F M F F F F F F F F F F F F F F F F F F F F F F F F F F F F F
## [596] F F F F F F F F F F F F F F F M F F F F F F F F F F F F F F F F F F F
## [631] F F F F F F F F F F F F F F F F F F M F F F F M M F F F F F F F F F F
## [666] F F F F F F F F M F F F F F F F F F M F F F F F F F F F F F F F F F F
## [701] F F F F F F F F F F F F F F M F F F F F F F F F F F F M M M M M M M M
## [736] M M M M M M M M M M M M M M M M M M M M M M M M M M M M M F M M M M M
## [771] F F M M M M M M M M F M M M F M M M M M M M M M M M M M M M M M M M M
## [806] M M M M M M M M M M M M M F M M M M M M M M M M M M M M M M M M F F F
## [841] M M F M M M M M M F M M M M M F M M M M M M M M F M M M M M M M M M M
## [876] M M M M M M M M M M M M M M M M F F F F F F F F F F F F F F F F F F F
## [911] F F F F F F F F F F F F F F F F F F F F F F F F F F F F M F F F M M F
## [946] M M M M M M M M M M M M M M M M M M M M M M M M M M M M M M M M M M M
## [981] M M M M M M M M
## Levels: F M
```

and the posterior probabilities for the classifications (we have not specified a prior, so it is flat by default I believe). Here we will just look at the first and last 10 observations.


```r
head(linDiscrim_1$posterior, n = 10)
```

```
##         F         M
## 1  0.9892 0.0107886
## 2  0.9692 0.0308273
## 3  0.9883 0.0116947
## 4  0.9804 0.0195656
## 5  0.9978 0.0022387
## 6  0.9976 0.0023530
## 7  0.9998 0.0001596
## 8  0.9996 0.0003768
## 9  0.9924 0.0075874
## 10 0.9981 0.0018559
```

```r
tail(linDiscrim_1$posterior, n = 10)
```

```
##             F      M
## 979 0.0008157 0.9992
## 980 0.0520666 0.9479
## 981 0.0060566 0.9939
## 982 0.2180089 0.7820
## 983 0.0132210 0.9868
## 984 0.0007304 0.9993
## 985 0.0120224 0.9880
## 986 0.0034559 0.9965
## 987 0.0064043 0.9936
## 988 0.0018918 0.9981
```


We can ask which ones did lda classify correctly, as compared to the true sex.


```r
linDiscrim_1$class == wings$Sex
```

```
##   [1]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
##  [12]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
##  [23]  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
##  [34]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE
##  [45]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE
##  [56]  TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE
##  [67] FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE
##  [78] FALSE FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
##  [89]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [100]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE
## [111]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [122]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE
## [133]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [144]  TRUE  TRUE  TRUE FALSE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE
## [155]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE
## [166]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [177]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [188]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [199]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [210]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [221]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [232]  TRUE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE
## [243]  TRUE  TRUE  TRUE FALSE FALSE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE
## [254]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [265]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [276]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [287]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [298]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [309]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [320]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [331]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [342]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE  TRUE
## [353]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE
## [364]  TRUE  TRUE FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [375]  TRUE  TRUE FALSE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE
## [386]  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [397]  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE
## [408]  TRUE FALSE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE
## [419]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [430]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [441]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [452]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [463] FALSE FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [474]  TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [485]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [496]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [507]  TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [518]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [529]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [540]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE  TRUE  TRUE  TRUE
## [551]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [562]  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [573]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [584]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [595]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [606]  TRUE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE
## [617]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [628]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [639]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE
## [650]  TRUE  TRUE  TRUE  TRUE FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE
## [661]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [672]  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [683]  TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [694]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [705]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE
## [716]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [727]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [738]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [749]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [760]  TRUE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE
## [771] FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE
## [782]  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [793]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [804]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [815]  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [826]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [837]  TRUE FALSE FALSE FALSE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE
## [848]  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE
## [859]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE
## [870]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [881]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [892]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [903]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [914]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [925]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [936]  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE FALSE  TRUE FALSE  TRUE
## [947]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [958]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [969]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
## [980]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
```


Importantly we can ask what proportion were correctly identified this way. Again, keep in mind this only used leave-one-out cross-validation, so it is probably not sufficient.


```r
prob_correct <- sum(linDiscrim_1$class == wings$Sex)/length(wings$Sex)
prob_correct
```

```
## [1] 0.9362
```

So the lda with 1-fold cross-validation gives us 93.6235 percent correct identification.

## Below I am just checking some of the math. It can be ignored for now

Computing the discriminant function coefficient vector, a For the two group problem.
First compute means for males and females


```r
males <- wings[wings$Sex == "M", 9:34]
females <- wings[wings$Sex == "F", 9:34]
male.means <- colMeans(males)
female.means <- colMeans(females)
```


Calculate pooled covariance matrix

```r
pooled.cov <- cov(wings[, 9:34])  # This is S
S.inv <- solve(pooled.cov)  # If you set the tolerance too low, all hell breaks loose, so for this example I am using a covariance matrix of full rank,
```


Now compute the discriminant function coefficient vector, a = S[-1]

```r
a = S.inv %*% (male.means - female.means)

# rescaling it to compare to the original vector
a.prime = a * (34.488708/31.918649)
cbind(a.prime, linDiscrim_0$scaling)
```

```
##                           LD1
## ProcCoord1   -34.489  -34.489
## ProcCoord2   177.933  177.933
## ProcCoord3    16.768   16.768
## ProcCoord4  -179.750 -179.750
## ProcCoord5    -6.883   -6.883
## ProcCoord6   109.968  109.968
## ProcCoord7  -112.176 -112.176
## ProcCoord8    87.017   87.017
## ProcCoord9    45.013   45.013
## ProcCoord10 -105.772 -105.772
## ProcCoord11  -43.168  -43.168
## ProcCoord12   86.987   86.987
## ProcCoord13  -80.005  -80.005
## ProcCoord14   35.794   35.794
## ProcCoord15   12.771   12.771
## ProcCoord16   51.418   51.418
## ProcCoord17  102.529  102.529
## ProcCoord18  -40.635  -40.635
## ProcCoord19    6.644    6.644
## ProcCoord20 -116.080 -116.080
## ProcCoord21  -81.341  -81.341
## ProcCoord22  247.453  247.453
## ProcCoord23   -8.985   -8.985
## ProcCoord24  267.234  267.234
## ProcCoord25   46.402   46.402
## ProcCoord26   12.366   12.366
```

