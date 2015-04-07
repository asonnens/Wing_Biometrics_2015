##This script includes several sets of analysis
##    1) analysis of wingmachine coordinates using 
##       principal components and a variety of statistical
##       learning methods, classifying by sex and genotype
##    2) re-analysis of same dataset using methods/settings
##       that were also used for Biocat
##    3) Code for wing blur 
##    4) Code for wing effect plot
##    4) Code for lda plots
##    5) heatmap of confusion matrices from Biocat and wingmachine
##       Biocat data is entered manually (was not recorded in a format 
##       compatible with R)

#necessary packages/libraries
#########################################################################################
require(MASS)     # contains lda() and qda() functions
require(sampling) # easier way to generate subsets of data
require(mda)      # flexible discriminant analysis and mixture discriminant analysis
require(class)    # k nearest neighbours (knn)
require(adabag)   # bagging()
require(ada)      # boosting function, ada()
require(nnet)     # neural net function, nnet()
require(e1071)    # svm()
require(randomForest) # randomForest() 
require(plyr)
#########################################################################################
 
#1.
#########################################################################################
#entering all 4X Olympus wing data, with subsetting
wings <- read.table("../data/Olympus_4X_coords.tsv", h=T) #all olympus 4X wings
sam_wings <- data.frame( wings[ wings$genotype == "samw" ,], row.names=c(1:416) ) #sam wings, both sexes
left_all_wings    <- data.frame( wings[ wings$side == "L" ,], row.names=c(1:1130) ) #left wings only
left_sam_wings    <- data.frame( wings[ wings$genotype == "samw" & wings$side == "L",], row.names=c(1:208) ) #left sam wings only
left_female_wings <- data.frame( wings[ wings$sex == "F" & wings$side == "L",], row.names=c(1:559) ) #left female wings only

#wings from other datasets
all_wings <- read.table("../data/all_wings.tsv", h = T)
wings_olympus_2x <- read.table("../data/Olympus_2X_coords.tsv", h=T)
wings_leica_4x <- read.table("../data/Leica_4X_coords.tsv", h=T)
wings_leica_2x <- read.table("../data/Leica_2X_coords.tsv", h=T)
all_wings_sam <- data.frame(all_wings[all_wings$genotype == "samw" ,], row.names=c(1:1638) )
o2x_sam_wings <- data.frame(wings_olympus_2x[wings_olympus_2x$genotype == "samw" ,], row.names=c(1:413) ) 
l4x_sam_wings <- data.frame(wings_leica_4x[wings_leica_4x$genotype == "samw" ,], row.names=c(1:398) )
l2x_sam_wings <- data.frame(wings_leica_2x[wings_leica_2x$genotype == "samw" ,], row.names=c(1:411) )

all_wings_oly <- data.frame(all_wings[all_wings$microscope == "oly" ,], row.names = c(1:4525))
all_wings_oly4X <- data.frame(all_wings_oly[all_wings_oly$magnification == "4X" ,], row.names = c(1:2264))
all_wings_oly4X_sam <- data.frame(all_wings_oly4X[all_wings_oly4X$genotype == "samw" ,], row.names = c(1:416))
all_wings_oly2X <- data.frame(all_wings_oly[all_wings_oly$magnification == "2X" ,], row.names = c(1:2261))
all_wings_oly2X_sam <- data.frame(all_wings_oly2X[all_wings_oly2X$genotype == "samw" ,], row.names = c(1:413))
all_wings_lei <- data.frame(all_wings[all_wings$microscope == "lei" ,], row.names = c(1:4422))
all_wings_lei4X <- data.frame(all_wings_lei[all_wings_lei$magnification == "4X" ,], row.names = c(1:2165))
all_wings_lei4X_sam <- data.frame(all_wings_lei4X[all_wings_lei4X$genotype == "samw" ,], row.names = c(1:398))
all_wings_lei2X <- data.frame(all_wings_lei[all_wings_lei$magnification == "2X" ,], row.names = c(1:2257))
all_wings_4X <- data.frame(all_wings[all_wings$magnification == "4X" ,], row.names = c(1:4429))
all_wings_4Xsam <- data.frame(all_wings_4X[all_wings_4X$genotype == "samw" ,], row.names = c(1:814))
all_wings_2X <- data.frame(all_wings[all_wings$magnification == "2X" ,], row.names = c(1:4518))

##this creates a second dataframe, with left and right wings from individual flies averaged into one value
wing_average <- aggregate( wings[,9:105], list( ID=wings$ID, genotype=wings$genotype, sex=wings$sex ), mean)
sam_wing_average <- aggregate( sam_wings[,9:105], list( ID=sam_wings$ID, genotype=sam_wings$genotype, sex=sam_wings$sex ), mean)
#averaging for combined datasets
all_wings_average <- aggregate( all_wings[,9:105], list( ID=all_wings$ID, genotype=all_wings$genotype, sex=all_wings$sex ), mean)
all_wings_sam_average <- aggregate( all_wings_sam[,9:105], list( ID=all_wings_sam$ID, genotype=all_wings_sam$genotype, sex=all_wings_sam$sex ), mean)
all_wings_4Xsam_average <- aggregate( all_wings_4Xsam[,9:105], list( ID=all_wings_4Xsam$ID, genotype=all_wings_4Xsam$genotype, sex=all_wings_4Xsam$sex ), mean)
all_wings_oly4Xsam_average <- aggregate( all_wings_oly4X_sam[,9:105], list( ID=all_wings_oly4X_sam$ID, genotype=all_wings_oly4X_sam$genotype, sex=all_wings_oly4X_sam$sex ), mean)
all_wings_oly2Xsam_average <- aggregate( all_wings_oly2X_sam[,9:105], list( ID=all_wings_oly2X_sam$ID, genotype=all_wings_oly2X_sam$genotype, sex=all_wings_oly2X_sam$sex ), mean)
all_wings_lei4Xsam_average <- aggregate( all_wings_lei4X_sam[,9:105], list( ID=all_wings_lei4X_sam$ID, genotype=all_wings_lei4X_sam$genotype, sex=all_wings_lei4X_sam$sex ), mean)
#averaging for other datasets
wing_average_olympus_2x <- aggregate( wings_olympus_2x[,9:105], list( ID=wings_olympus_2x$ID, genotype=wings_olympus_2x$genotype, sex=wings_olympus_2x$sex ), mean)
wing_average_leica_4x <- aggregate( wings_leica_4x[,9:105], list( ID=wings_leica_4x$ID, genotype=wings_leica_4x$genotype, sex=wings_leica_4x$sex ), mean)
wing_average_leica_2x <- aggregate( wings_leica_2x[,9:105], list( ID=wings_leica_2x$ID, genotype=wings_leica_2x$genotype, sex=wings_leica_2x$sex ), mean)
sam_wing_average_olympus_2x <- aggregate( o2x_sam_wings[,9:105], list( ID=o2x_sam_wings$ID, genotype=o2x_sam_wings$genotype, sex=o2x_sam_wings$sex ), mean)
sam_wing_average_leica_4x <- aggregate( l4x_sam_wings[,9:105], list( ID=l4x_sam_wings$ID, genotype=l4x_sam_wings$genotype, sex=l4x_sam_wings$sex ), mean)
sam_wing_average_leica_2x <- aggregate( l2x_sam_wings[,9:105], list( ID=l2x_sam_wings$ID, genotype=l2x_sam_wings$genotype, sex=l2x_sam_wings$sex ), mean)


#Append PCs to new dataframe
#all wings from olympus 4X dataset			   
wings_mean_pc <- prcomp(wing_average[,4:99])
summary(wings_mean_pc)
head(wings_mean_pc$x[,1:58])
average_wings <- data.frame(wing_average, wings_mean_pc$x[,1:58])
#sam wings from olympus 4X dataset
sam_wings_mean_pc <- prcomp(sam_wing_average[,4:99])
sam_average_wings <- data.frame(sam_wing_average, sam_wings_mean_pc$x[,1:58])


#appending PCs to dataframe for wings from other datasets

#combined datasets
#all wings
all_wings_mean_pc <- prcomp(all_wings_average[,4:99])
all_wings_average_wings <- data.frame(all_wings_average, all_wings_mean_pc$x[,1:58])
#all sam wings
all_wings_sam_mean_pc <- prcomp(all_wings_sam_average[,4:99])
all_wings_sam_average_wings <- data.frame(all_wings_sam_average, all_wings_sam_mean_pc$x[,1:58])
#all 4X sam wings
all_wings_4Xsam_average_mean_pc <- prcomp(all_wings_4Xsam_average[,4:99])
all_4Xsam <- data.frame(all_wings_4Xsam_average, all_wings_4Xsam_average_mean_pc$x[,1:58])
#all oly4X sam wings
all_wings_oly4Xsam_average_mean_pc <- prcomp(all_wings_oly4Xsam_average[,4:99])
all_oly4Xsam <- data.frame(all_wings_oly4Xsam_average, all_wings_oly4Xsam_average_mean_pc$x[,1:58])
#all oly2X sam wings
all_wings_oly2Xsam_average_mean_pc <- prcomp(all_wings_oly2Xsam_average[,4:99])
all_oly2Xsam <- data.frame(all_wings_oly2Xsam_average, all_wings_oly2Xsam_average_mean_pc$x[,1:58])
#all lei4X sam wings
all_wings_lei4Xsam_average_mean_pc <- prcomp(all_wings_lei4Xsam_average[,4:99])
all_lei4Xsam <- data.frame(all_wings_lei4Xsam_average, all_wings_lei4Xsam_average_mean_pc$x[,1:58])


#individual datasets
#olympus 2x wings
olympus_2x_wings_mean_pc <- prcomp(wing_average_olympus_2x[,4:99])
olympus_2x_average_wings <- data.frame(wing_average_olympus_2x, olympus_2x_wings_mean_pc$x[,1:58])
#sam olympus 2x wings
sam_olympus_2x_wings_mean_pc <- prcomp(sam_wing_average_olympus_2x[,4:99])
olympus_2x_sam_average_wings <- data.frame(sam_wing_average_olympus_2x, sam_olympus_2x_wings_mean_pc$x[,1:58])
#all wings
leica_4x_wings_mean_pc <- prcomp(wing_average_leica_4x[,4:99])
leica_4x_average_wings <- data.frame(wing_average_leica_4x, leica_4x_wings_mean_pc$x[,1:58])
#sam wings
sam_leica_4x_wings_mean_pc <- prcomp(sam_wing_average_leica_4x[,4:99])
leica_4x_sam_average_wings <- data.frame(sam_wing_average_leica_4x, sam_leica_4x_wings_mean_pc$x[,1:58])
#all wings
leica_2x_wings_mean_pc <- prcomp(wing_average_leica_2x[,4:99])
leica_2x_average_wings <- data.frame(wing_average_leica_2x, leica_2x_wings_mean_pc$x[,1:58])
#sam wings
sam_leica_2x_wings_mean_pc <- prcomp(sam_wing_average_leica_2x[,4:99])
leica_2x_sam_average_wings <- data.frame(sam_wing_average_leica_2x, sam_leica_2x_wings_mean_pc$x[,1:58])

#appending PCs, including centroid- just olympus 4X dataset				
#all wings
wings_mean_pc_cent <- prcomp(wing_average[,4:100])
summary(wings_mean_pc_cent)
head(wings_mean_pc_cent$x[,1:58])
average_wings_cent <- data.frame(wing_average, wings_mean_pc_cent$x[,1:58])
#averaged sam wings
sam_wings_mean_pc_cent <- prcomp(sam_wing_average[,4:100])
sam_average_wings_cent <- data.frame(sam_average_wings, sam_wings_mean_pc_cent$x[,1:58])

#sam wings from olympus 4X dataset 
#here wings are not averaged, because these subsets are meant for comparison to BioCAT
sam_pc <- prcomp(sam_wings[,9:104]) #excluding centroid
sam_data <- data.frame(sam_wings, sam_pc$x[,1:58])
#left wings, both sexes
left_wings_pc <- prcomp(left_all_wings[,9:104])
left_all_wings_data <- data.frame(left_all_wings, left_wings_pc$x[,1:58])
left_wings_pc_cent <- prcomp(left_all_wings[,9:105])
left_all_wings_cent <- data.frame(left_all_wings, left_wings_pc_cent$x[,1:58])
#left wings, females
left_female_wings_pc <- prcomp(left_female_wings[,9:104])
left_female_wings_data <- data.frame(left_female_wings, left_female_wings_pc$x[,1:58])
left_female_wings_pc_cent <- prcomp(left_female_wings[,9:105])
left_female_wings_cent <- data.frame(left_female_wings, left_female_wings_pc_cent$x[,1:58])
#left wings, sam
left_sam_wings_pc <- prcomp(left_sam_wings[,9:104])
left_sam_wings_data <- data.frame(left_sam_wings, left_sam_wings_pc$x[,1:58])
left_sam_wings_pc_cent <- prcomp(left_sam_wings[,9:105])
left_sam_wings_cent <- data.frame(left_sam_wings, left_sam_wings_pc_cent$x[,1:58])


#########################################################################################
#Using morphometric PCs to predict sex and genotype
#########################################################################################


#setting tables
#########################################################################################
#tables for sex prediction - olympus 4X dataset
tabw <- table(average_wings$sex)
tabwc <- table(average_wings_cent$sex)
tablw <- table(left_all_wings_data$sex)
tablwc <- table(left_all_wings_cent$sex)
tabls <- table(left_sam_wings_data$sex)
tablsc <- table(left_sam_wings_cent$sex)
tab_sam <- table(sam_wings$sex)
tab_a_sam <- table(sam_average_wings$sex)
tabw 
tabwc
tablw
tablwc
tabls
tablsc
tab_sam
tab_a_sam

#tables for sex prediction - other datasets
#combined datasets
tab_all <- table(all_wings_average$sex)
tab_all_4Xsam <- table(all_wings_4Xsam_average$sex)
tab_alloly4Xsam <- table(all_wings_oly4Xsam_average$sex)
tab_alloly2Xsam <- table(all_wings_oly2Xsam_average$sex)
tab_alllei4Xsam <- table(all_wings_lei4Xsam_average$sex)
#individual datasets
tabwo2x <- table(olympus_2x_average_wings$sex)
tabwl4x <- table(leica_4x_average_wings$sex)
tabwl2x <- table(leica_2x_average_wings$sex)
tabwo2x
tabwl4x
tabwl2x
#sam averaged wings
tab_a_wo2x <- table(olympus_2x_sam_average_wings$sex)
tab_a_wl4x <- table(leica_4x_sam_average_wings$sex)
tab_a_wl2x <- table(leica_2x_sam_average_wings$sex)
tab_a_wo2x
tab_a_wl4x
tab_a_wl2x


#tables for genotype prediction, all genotypes, - olympus 4X dataset
tabw_g <- table(average_wings$genotype)
tabwc_g <- table(average_wings_cent$genotype)
tablw_g <- table(left_female_wings_data$genotype)
tablwc_g <- table(left_female_wings_cent$genotype)
tabw_g
tabwc_g
tablw_g
tablwc_g
#########################################################################################

#testing and training sets for sex-- strata is wrapped in separate function to save space
#########################################################################################
strata_2var <- function(dataset, variable_name, variable_table, percentage, variable_pos, dataframe_length = 162){
#this function stratifies a dataset by two variables, 
#and makes testing and training sets based on the strata
#testing and training sets returned as a list
  start_position <- dataframe_length - 57
  my_strata <- strata(dataset, stratanames = c(variable_name), method = "srswor", size = 
    c(round(variable_table[1] *(percentage)), round(variable_table[2] * (percentage))))
  my_data <- dataset[, c(variable_pos, start_position:dataframe_length)]
  my_training <- my_data[rownames(my_data) %in% my_strata$ID_unit, ]
  my_testing <- my_data[!rownames(my_data) %in% my_strata$ID_unit, ] 
  my_train_test_info <- list(my_training, my_testing)
  return(my_train_test_info)
}

#creating testing and training sets, genotype
#########################################################################################
strata_5var <- function(dataset, variable_name, variable_table, percentage, variable_pos, dataframe_length = 162){
#this function stratifies a dataset by five variables, 
#and makes testing and training sets based on the strata
#testing and training sets returned as a list
  start_position <- dataframe_length - 57
  my_strata <- strata(dataset, stratanames = c(variable_name), method = "srswor", 
    size = c(round(variable_table[1] *(percentage)), 
	  round(variable_table[2] * (percentage)),
	  round(variable_table[3] * (percentage)),
	  round(variable_table[4] * (percentage)),
	  round(variable_table[5] * (percentage))))
  my_data <- dataset[, c(variable_pos, start_position:dataframe_length)]
  my_training <- my_data[rownames(my_data) %in% my_strata$ID_unit, ]
  my_testing <- my_data[!rownames(my_data) %in% my_strata$ID_unit, ] 
  my_train_test_info <- list(my_training, my_testing)
  return(my_train_test_info)
}

#########################################################################################

#Function for re-sampling the training and testing set, for obtaining error 
#estimates. This is intended for functions run on default settings (below)
#functions with variable parameters have unique re-sampling functions
#all of the resampling functions output:
# 1) average classification accuracy
# 2) stdev of classification accuracy 
# 3) list of classification accuracy from each file 

resample_default <- function(function_name, dataset, datatable, reps){
	out_list <- list()
    for(i in 1:reps){
        myresult <- function_name(dataset, datatable)
		out_list <- unlist(c(out_list, myresult))
     }
	mean_score <- mean(out_list)
	error_score <- sd(out_list)
	return_list <- list(mean_score,error_score,out_list) 
	return(return_list)
}

#this resample function is specifically for training on one dataset, testing on another
#the only way it really differs from the default is that it requires more arguments
resample_default_compare <- function(function_name, dataset1, dataset2, datatable1, datatable2, reps){
	out_list <- list()
    for(i in 1:reps){
        myresult <- function_name(dataset1, dataset2, datatable1, datatable2)
		out_list <- unlist(c(out_list, myresult))
     }
	mean_score <- mean(out_list)
	error_score <- sd(out_list)
	return_list <- list(mean_score,error_score,out_list) 
	return(return_list)
}


#Tests of how well various statistical learning methods can classify wings by sex
#These are done on default parameters (lda,qda,fda,mda,bagging, neural network)
#these are done within the sam (wildtype) genotype, left and right wings
#these use shape information only (no centroid)

#linear discriminant analysis
lda_repeat_function_sex <- function(dataset, datatable){
    lda_wings_sex <- strata_2var(dataset, "sex", datatable, 2/3, 3, 158)
	train_set <- lda_wings_sex[[1]]
	test_set <- lda_wings_sex[[2]]
        linDiscrim <- lda(formula=sex ~.,data = train_set, tol = 1.0e-8, CV=FALSE) #modified tolerance
    lda_training_table <- table(actual = train_set$sex,
                          predicted = predict(linDiscrim, newdata=train_set)$class)                         
    lda_test_table <- table(actual = test_set$sex,
                          predicted = predict(linDiscrim, newdata=test_set)$class)  
    prediction_success <- 100*sum(diag(lda_test_table)/sum(lda_test_table)) 
	return(prediction_success)
	}

#linear discriminant analysis comparison function
#trains LDA on one dataset, uses model to predict classfications for another	
compare_lda_repeat_function_sex <- function(dataset1, dataset2, datatable1, datatable2){
    lda_wings_sex_trainset <- strata_2var(dataset1, "sex", datatable1, 2/3, 3, 158)
	lda_wings_sex_testset <- strata_2var(dataset2, "sex", datatable2, 2/3, 3, 158)
	train_set <- lda_wings_sex_trainset[[1]]
	test_set <- lda_wings_sex_testset[[2]]
        linDiscrim <- lda(formula=sex ~.,data = train_set, tol = 1.0e-8, CV=FALSE) #modified tolerance
    lda_training_table <- table(actual = train_set$sex,
                          predicted = predict(linDiscrim, newdata=train_set)$class)                         
    lda_test_table <- table(actual = test_set$sex,
                          predicted = predict(linDiscrim, newdata=test_set)$class)  
    prediction_success <- 100*sum(diag(lda_test_table)/sum(lda_test_table)) 
	return(prediction_success)
	}

#quadratic discriminant analysis
qda_repeat_function_sex <- function(dataset, datatable){
    qda_wings_sex <- strata_2var(dataset, "sex", datatable, 2/3, 3, 158)
	train_set <- qda_wings_sex[[1]]
	test_set <- qda_wings_sex[[2]]
        QDiscrim <- qda(formula=sex ~.,data = train_set, tol = 1.0e-8, CV=FALSE) #modified tolerance
    qda_training_table <- table(actual = train_set$sex,
                          predicted = predict(QDiscrim, newdata=train_set)$class)                         
    qda_test_table <- table(actual = test_set$sex,
                          predicted = predict(QDiscrim, newdata=test_set)$class)  
    prediction_success <- 100*sum(diag(qda_test_table)/sum(qda_test_table)) 
	return(prediction_success)
	}

#flexible discriminant analysis
fda_repeat_function_sex <- function(dataset, datatable){
    fda_wings_sex <- strata_2var(dataset, "sex", datatable, 2/3, 3, 158)
	train_set <- fda_wings_sex[[1]]
	test_set <- fda_wings_sex[[2]]
    fda_model <- fda(formula=sex ~.,data = train_set)
    fda_training_table <- table(actual = train_set$sex,
                          predicted = predict(fda_model, newdata=train_set, type="class"))                     
    fda_test_table <- table(actual = test_set$sex,
                          predicted = predict(fda_model, newdata=test_set, type="class"))
    prediction_success <- 100*sum(diag(fda_test_table)/sum(fda_test_table))
	return(prediction_success)
	}

#multiple discriminant analysis
mda_repeat_function_sex <- function(dataset, datatable){
    mda_wings_sex <- strata_2var(dataset, "sex", datatable, 2/3, 3, 158)
	train_set <- mda_wings_sex[[1]]
	test_set <- mda_wings_sex[[2]]
    mda_model <- mda(formula=sex ~.,data = train_set)
    mda_training_table <- table(actual = train_set$sex,
                          predicted = predict(mda_model, newdata=train_set, type="class"))                     
    mda_test_table <- table(actual = test_set$sex,
                          predicted = predict(mda_model, newdata=test_set, type="class"))
    prediction_success <- 100*sum(diag(mda_test_table)/sum(mda_test_table))
	return(prediction_success)
	}

#bootstrap aggregation
bagging_repeat_function_sex <- function(dataset, datatable){
    bagging_wings_sex <- strata_2var(dataset, "sex", datatable, 2/3, 3, 158)
	train_set <- bagging_wings_sex[[1]]
	test_set <- bagging_wings_sex[[2]]
	bag_model <- bagging(formula=sex ~.,data = train_set)
	bag_training_table <- table(actual = train_set$sex,
                          predicted = predict(bag_model, newdata=train_set)$class)
	bag_test_table <- table(actual = test_set$sex,
                          predicted = predict(bag_model, newdata=test_set)$class)
    prediction_success <- 100*sum(diag(bag_test_table)/sum(bag_test_table))
    return(prediction_success)
}

#neural network
nnet_repeat_function_sex <- function(dataset, datatable){
    nnet_wings_sex <- strata_2var(dataset, "sex", datatable, 2/3, 3, 158)
	train_set <- nnet_wings_sex[[1]]
	test_set <- nnet_wings_sex[[2]]
	nnet_model <- nnet(formula=sex ~.,data = train_set, size=10, decay=0)
    nnet_training_table <- table(actual = train_set$sex,
                          predicted = predict(nnet_model, type="class"))                         
    nnet_test_table <- table(actual = test_set$sex,
                          predicted = predict(nnet_model, newdata=test_set, type="class"))
    prediction_success <- 100*sum(diag(nnet_test_table)/sum(nnet_test_table))
	return(prediction_success)
	}
	
#running above functions on default parameters
#currently set to a very LOW number of repetitions (5). 
#1000 reps of the bagging function takes several hours to run

#using olympus 4X dataset, averaged sams wings
lda_score_CI_sex <- resample_default(lda_repeat_function_sex, sam_average_wings, tab_a_sam, 5) 
lda_score_CI_sex
qda_score_CI_sex <- resample_default(qda_repeat_function_sex, sam_average_wings, tab_a_sam, 5)
qda_score_CI_sex
fda_score_CI_sex <- resample_default(fda_repeat_function_sex, sam_average_wings, tab_a_sam, 5)
fda_score_CI_sex
mda_score_CI_sex <- resample_default(mda_repeat_function_sex, sam_average_wings, tab_a_sam, 5)
mda_score_CI_sex
bag_score_CI_sex <- resample_default(bagging_repeat_function_sex, sam_average_wings, tab_a_sam, 5)
bag_score_CI_sex
nnet_score_CI_sex <- resample_default(nnet_repeat_function_sex, sam_average_wings, tab_a_sam, 5)
nnet_score_CI_sex

#using LDA on other datasets
#Olympus 2X
lda_score_CI_sex_olympus2x <- resample_default(lda_repeat_function_sex, olympus_2x_average_wings, tab_a_wo2x, 5) 
lda_score_CI_sex_olympus2x
#Leica 4X
lda_score_CI_sex_leica4x <- resample_default(lda_repeat_function_sex, leica_4x_average_wings, tab_a_wl4x, 5) 
lda_score_CI_sex_leica4x
#Leica 2X
lda_score_CI_sex_leica2x <- resample_default(lda_repeat_function_sex, leica_2x_average_wings, tab_a_wl2x, 5) 
lda_score_CI_sex_leica2

#using LDA on combined datasets
#Olympus 4X
lda_score_CI_sex_oly4x_combined <- resample_default(lda_repeat_function_sex, all_oly4Xsam , tab_alloly4Xsam, 5) 
lda_score_CI_sex_oly4x_combined

#using an LDA trained on one dataset to classify other datasets

#using combined datasets
#train with oly4X, test with lei2X
#THIS values output here seem very unlikely
lda_alloly4X_alloly2X <- resample_default_compare(compare_lda_repeat_function_sex, all_oly4Xsam, all_oly2Xsam, tab_alloly4Xsam, tab_alloly2Xsam, 5)
lda_alloly4X_alloly2X
lda_alloly4X_alllei4X <- resample_default_compare(compare_lda_repeat_function_sex, all_oly4Xsam, all_lei4Xsam, tab_alloly4Xsam, tab_alllei4Xsam, 5)
lda_alloly4X_alllei4X
lda_all4X_alloly4X <- resample_default_compare(compare_lda_repeat_function_sex, all_4Xsam, all_oly4Xsam, tab_all_4Xsam, tab_alloly4Xsam, 5)
lda_all4X_alloly4X
lda_all4X_alllei4X <- resample_default_compare(compare_lda_repeat_function_sex, all_4Xsam, all_lei4Xsam, tab_all_4Xsam, tab_alllei4Xsam, 5)
lda_all4X_alllei4X


#using individual datasets
#train with oly4X, test with oly2X
lda_oly4X_oly2X <- resample_default_compare(compare_lda_repeat_function_sex, sam_average_wings, olympus_2x_average_wings, tab_a_sam, tab_a_wo2x, 5)
lda_oly4X_oly2X
#train with oly4X, test with lei4X
lda_oly4X_lei4X <- resample_default_compare(compare_lda_repeat_function_sex, sam_average_wings, leica_4x_average_wings, tab_a_sam, tab_a_wl4x, 5)
lda_oly4X_lei4X
#train with oly4X, test with lei2X
lda_oly4X_lei2X <- resample_default_compare(compare_lda_repeat_function_sex, sam_average_wings, leica_2x_average_wings, tab_a_sam, tab_a_wl2x, 5)
lda_oly4X_lei2X
#train with lei4X, test with oly4X
lda_lei4X_oly4X <- resample_default_compare(compare_lda_repeat_function_sex, leica_4x_average_wings, sam_average_wings, tab_a_wl4x, tab_a_sam, 5)
lda_lei4X_oly4X


#Tests of how well various statistical learning methods can classify wings by genotype
#These are done on default parameters (lda,qda,fda,mda,bagging, neural network)
#these are done with all wings, left and right wings averaged
#these use shape information only (no centroid)

#linear discriminant analysis
lda_repeat_function_genotype <- function(dataset, datatable){
    lda_wings_genotype <- strata_5var(dataset, "genotype", datatable, 2/3, 2, 158)
	train_set <- lda_wings_genotype[[1]]
	test_set <- lda_wings_genotype[[2]]
        linDiscrim <- lda(formula=genotype ~.,data = train_set, tol = 1.0e-8, CV=FALSE) #modified tolerance
    lda_training_table <- table(actual = train_set$genotype,
                          predicted = predict(linDiscrim, newdata=train_set)$class)                         
    lda_test_table <- table(actual = test_set$genotype,
                          predicted = predict(linDiscrim, newdata=test_set)$class)  
    prediction_success <- 100*sum(diag(lda_test_table)/sum(lda_test_table)) 
	return(prediction_success)
	}

#quadratic discriminant analysis	
qda_repeat_function_genotype <- function(dataset, datatable){
    qda_wings_genotype <- strata_5var(dataset, "genotype", datatable, 2/3, 2, 158)
	train_set <- qda_wings_genotype[[1]]
	test_set <- qda_wings_genotype[[2]]
        QDiscrim <- qda(formula=genotype ~.,data = train_set, tol = 1.0e-8, CV=FALSE) #modified tolerance
    qda_training_table <- table(actual = train_set$genotype,
                          predicted = predict(QDiscrim, newdata=train_set)$class)                         
    qda_test_table <- table(actual = test_set$genotype,
                          predicted = predict(QDiscrim, newdata=test_set)$class)  
    prediction_success <- 100*sum(diag(qda_test_table)/sum(qda_test_table)) 
	return(prediction_success)
	}

#flexible discriminant analysis
fda_repeat_function_genotype <- function(dataset, datatable){
    fda_wings_genotype <- strata_5var(dataset, "genotype", datatable, 2/3, 2, 158)
	train_set <- fda_wings_genotype[[1]]
	test_set <- fda_wings_genotype[[2]]
    fda_model <- fda(formula=genotype ~.,data = train_set)
    fda_training_table <- table(actual = train_set$genotype,
                          predicted = predict(fda_model, newdata=train_set, type="class"))                     
    fda_test_table <- table(actual = test_set$genotype,
                          predicted = predict(fda_model, newdata=test_set, type="class"))
    prediction_success <- 100*sum(diag(fda_test_table)/sum(fda_test_table))
	return(prediction_success)
	}

#multiple discriminant analysis
mda_repeat_function_genotype <- function(dataset, datatable){
    mda_wings_genotype <- strata_5var(dataset, "genotype", datatable, 2/3, 2, 158)
	train_set <- mda_wings_genotype[[1]]
	test_set <- mda_wings_genotype[[2]]
    mda_model <- mda(formula=genotype ~.,data = train_set)
    mda_training_table <- table(actual = train_set$genotype,
                          predicted = predict(mda_model, newdata=train_set, type="class"))                     
    mda_test_table <- table(actual = test_set$genotype,
                          predicted = predict(mda_model, newdata=test_set, type="class"))
    prediction_success <- 100*sum(diag(mda_test_table)/sum(mda_test_table))
	return(prediction_success)
	}

#bootstrap aggregation
bagging_repeat_function_genotype <- function(dataset, datatable){
    bagging_wings_genotype <- strata_5var(dataset, "genotype", datatable, 2/3, 2, 158)
	train_set <- bagging_wings_genotype[[1]]
	test_set <- bagging_wings_genotype[[2]]
	bag_model <- bagging(formula=genotype ~.,data = train_set)
	bag_training_table <- table(actual = train_set$genotype,
                          predicted = predict(bag_model, newdata=train_set)$class)
	bag_test_table <- table(actual = test_set$genotype,
                          predicted = predict(bag_model, newdata=test_set)$class)
    prediction_success <- 100*sum(diag(bag_test_table)/sum(bag_test_table))
    return(prediction_success)
	}

#neural network
nnet_repeat_function_genotype <- function(dataset, datatable){
    nnet_wings_genotype <- strata_5var(dataset, "genotype", datatable, 2/3, 2, 158)
	train_set <- nnet_wings_genotype[[1]]
	test_set <- nnet_wings_genotype[[2]]
	nnet_model <- nnet(formula=genotype ~.,data = train_set, size=10, decay=0)
    nnet_training_table <- table(actual = train_set$genotype,
                          predicted = predict(nnet_model, type="class"))                         
    nnet_test_table <- table(actual = test_set$genotype,
                          predicted = predict(nnet_model, newdata=test_set, type="class"))
    prediction_success <- 100*sum(diag(nnet_test_table)/sum(nnet_test_table))
	return(prediction_success)
	}
	

	
#running above functions on default parameters
#currently set to a very LOW number of repetitions (5). 
#1000 reps of the bagging function takes several hours to run

lda_score_CI_genotype <- resample_default(lda_repeat_function_genotype, average_wings, tabw_g, 5) 
lda_score_CI_genotype
qda_score_CI_genotype <- resample_default(qda_repeat_function_genotype, average_wings, tabw_g, 5) 
qda_score_CI_genotype
fda_score_CI_genotype <- resample_default(fda_repeat_function_genotype, average_wings, tabw_g, 5) 
fda_score_CI_genotype
mda_score_CI_genotype <- resample_default(mda_repeat_function_genotype, average_wings, tabw_g, 5) 
mda_score_CI_genotype
bagging_score_CI_genotype <- resample_default(bagging_repeat_function_genotype, average_wings, tabw_g, 1000) 
bagging_score_CI_genotype
nnet_score_CI_genotype <- resample_default(nnet_repeat_function_genotype, average_wings, tabw_g, 5) 
nnet_score_CI_genotype

#Tests of how well various statistical learning methods can classify wings by sex
#see script Wing_Biometrics_settings.r
#Functions with tweaked parameters
#knn: for specific values of k
#random forest: variable numbers of trees
#SVM: kernel shapes
#############################remind Anne to explain small error
#function for looping K nearest neighbours
knn_accuracy_loop_sex <- function(dataset, datatable,k_val = 4){ 
  knn_wings_sex <- strata_2var(dataset, "sex", datatable, 2/3, 3, 158)
  train_set <- knn_wings_sex[[1]]
  test_set <- knn_wings_sex[[2]]
  knn_fun <- knn(train = train_set[, -1], test = test_set[, -1], cl = train_set$sex, k = k_val, prob=TRUE);
  knn_table <- table(predicted=knn_fun, actual = test_set$sex);
  accuracy <- 100 * sum(diag(knn_table)/sum(knn_table));
  knn_info <- list(accuracy)
  return(knn_info)
}

knn_accuracy_loop_genotype <- function(dataset, datatable,k_val =32){ 
  knn_wings_genotype <- strata_5var(dataset, "genotype", datatable, 2/3, 2, 158)
  train_set <- knn_wings_genotype[[1]]
  test_set <- knn_wings_genotype[[2]]
  knn_fun <- knn(train = train_set[, -1], test = test_set[, -1], cl = train_set$genotype, k = k_val, prob=TRUE);
  knn_table <- table(predicted=knn_fun, actual = test_set$genotype);
  accuracy <- 100 * sum(diag(knn_table)/sum(knn_table));
  knn_info <- list(accuracy)
  return(knn_info)
}


#random forest loop function for sex
#100 trees generally has the best results
#makes a training and testing table
#outputs overall prediction success percentage
rf_sex_fun <- function(dataset, datatable, trees = 100){
    rf_wings_sex <- strata_2var(dataset, "sex", datatable, 2/3, 3, 158)
	rf_training <- rf_wings_sex[[1]]
	rf_test <- rf_wings_sex[[2]]
	random_forest_model <- randomForest(sex ~., data = rf_training, ntree = trees)
	rf_training_table <- table(actual = rf_training$sex,
                       predicted = predict(random_forest_model, type="class"))
	rf_test_table <- table(actual = rf_test$sex,
                       predicted = predict(random_forest_model, newdata=rf_test, type="class"))
	prediction_success <- 100 * sum(diag(rf_test_table)/sum(rf_test_table))
	return(prediction_success)
}

rf_genotype_fun <- function(dataset, datatable, trees = 1000){
    rf_wings_genotype <- strata_5var(dataset, "genotype", datatable, 2/3, 2, 158)
	rf_training <- rf_wings_genotype[[1]]
	rf_test <- rf_wings_genotype[[2]]
	random_forest_model <- randomForest(genotype ~., data = rf_training, ntree = trees)
	rf_training_table <- table(actual = rf_training$genotype,
                       predicted = predict(random_forest_model, type="class"))
	rf_test_table <- table(actual = rf_test$genotype,
                       predicted = predict(random_forest_model, newdata=rf_test, type="class"))
	prediction_success <- 100 * sum(diag(rf_test_table)/sum(rf_test_table))
	return(prediction_success)
}

#support vector machines
svm_sex_fun <- function(dataset, datatable, kernel_type = "sigmoid"){
    svm_wings_sex <- strata_2var(dataset, "sex", datatable, 2/3, 3, 158)
	svm_training <- svm_wings_sex[[1]]
	svm_test <- svm_wings_sex[[2]]
	svm_model <- svm(sex ~., data = svm_training, kernel = kernel_type)
	svm_training_table <- table(actual = svm_training$sex,
                       predicted = predict(svm_model, type="class"))
	svm_test_table <- table(actual = svm_test$sex,
                       predicted = predict(svm_model, newdata=svm_test, type="class"))
	prediction_success <- 100 * sum(diag(svm_test_table)/sum(svm_test_table))
	return(prediction_success)
}

svm_genotype_fun <- function(dataset, datatable, kernel_type = "radial"){
    svm_wings_genotype <- strata_5var(dataset, "genotype", datatable, 2/3, 2, 158)
	svm_training <- svm_wings_genotype[[1]]
	svm_test <- svm_wings_genotype[[2]]
	svm_model <- svm(genotype ~., data = svm_training, kernel = kernel_type)
	svm_training_table <- table(actual = svm_training$genotype,
                       predicted = predict(svm_model, type="class"))
	svm_test_table <- table(actual = svm_test$genotype,
                       predicted = predict(svm_model, newdata=svm_test, type="class"))
	prediction_success <- 100 * sum(diag(svm_test_table)/sum(svm_test_table))
	return(prediction_success)
}

#running classification methods with tweaked settings
#k nearest neighbours
knn_sex_info <- resample_default(knn_accuracy_loop_sex, sam_average_wings, tab_a_sam, 5) #best k value is 3-5
knn_genotype_info <- resample_default(knn_accuracy_loop_genotype, average_wings, tabw_g, 5) #best k value 32
#random forest
rf_sex_info <- resample_default(rf_sex_fun, sam_average_wings, tab_a_sam, 5) #best accuracy with 100 trees 
rf_genotype_info <- resample_default(rf_genotype_fun, average_wings, tabw_g, 5) #best accuracy with 1000 trees
#support vector machines
svm_sex_info <- resample_default(svm_sex_fun, sam_average_wings, tab_a_sam, 5) #best k value is 3-5
svm_genotype_info <- resample_default(svm_genotype_fun, average_wings, tabw_g, 5)


#Biocat comparisons
#this differs from the above classifications in that it's done on a small subset of 
#wing images--it's also performed on wing data with and without centroid, to test how 
#including size information changes classification accuracy
#only algorithms used are random forest and SVM
#########################################################################################

#support vector machines-- linear kernel
svm_sex_comparison <- function(dataset, datatable, kernel_type = "linear"){
    svm_wings_sex <- strata_2var(dataset, "sex", datatable, 2/3, 2, 162)
	svm_training <- svm_wings_sex[[1]]
	svm_test <- svm_wings_sex[[2]]
	svm_model <- svm(sex ~., data = svm_training, kernel = kernel_type)
	svm_training_table <- table(actual = svm_training$sex,
                       predicted = predict(svm_model, type="class"))
	svm_test_table <- table(actual = svm_test$sex,
                       predicted = predict(svm_model, newdata=svm_test, type="class"))
	prediction_success <- 100 * sum(diag(svm_test_table)/sum(svm_test_table))
	return(prediction_success)
}

svm_genotype_comparison <- function(dataset, datatable, kernel_type = "linear"){
    svm_wings_genotype <- strata_5var(dataset, "genotype", datatable, 2/3, 1, 162)
	svm_training <- svm_wings_genotype[[1]]
	svm_test <- svm_wings_genotype[[2]]
	svm_model <- svm(genotype ~., data = svm_training, kernel = kernel_type)
	svm_training_table <- table(actual = svm_training$genotype,
                       predicted = predict(svm_model, type="class"))
	svm_test_table <- table(actual = svm_test$genotype,
                       predicted = predict(svm_model, newdata=svm_test, type="class"))
	prediction_success <- 100 * sum(diag(svm_test_table)/sum(svm_test_table))
	return(prediction_success)
}

#left SAM wings, for sex comparisons
#no centroid
svm_sex_comparison_info <- resample_default(svm_sex_comparison, left_sam_wings_data, tabls, 5) 
#with centroid
svm_sex_comparison_info_centroid <- resample_default(svm_sex_comparison, left_sam_wings_cent, tablsc, 5)
#left female wings, for sex comparisons
#no centroid
svm_genotype_comparison_info <- resample_default(svm_genotype_comparison, left_female_wings_data, tablw_g, 5)
#with centroid
svm_genotype_comparison_info_centroid <- resample_default(svm_genotype_comparison, left_female_wings_cent, tablwc_g, 5)


#random forest comparison, for sex
#returns results for 10 trees, 1000 trees on each sample of training and test set
rf_sex_comparison <- function(dataset, datatable){
    trees = list(10, 1000)
    rf_wings_sex <- strata_2var(dataset, "sex", datatable, 2/3, 2, 158)
	rf_training <- rf_wings_sex[[1]]
	rf_test <- rf_wings_sex[[2]]
	results_list <- list()
	for (i in trees){
	    random_forest_model <- randomForest(sex ~., data = rf_training, ntree = i)
	    rf_training_table <- table(actual = rf_training$sex,
                       predicted = predict(random_forest_model, type="class"))
	    rf_test_table <- table(actual = rf_test$sex,
                       predicted = predict(random_forest_model, newdata=rf_test, type="class"))
	    prediction_success <- 100 * sum(diag(rf_test_table)/sum(rf_test_table))
		results_list <- list(unlist(results_list), prediction_success)
		}
	return(results_list)
}

#random forest comparison, for genotype
#returns results for 10 trees, 1000 trees on each sample of training and test set
rf_genotype_comparison <- function(dataset, datatable){
    trees = list(10, 1000)
    rf_wings_genotype <- strata_5var(dataset, "genotype", datatable, 2/3, 1, 162) #different dimensions, not average_wings dataframe
	rf_training <- rf_wings_genotype[[1]]
	rf_test <- rf_wings_genotype[[2]]
	results_list <- list()
	for (i in trees){
	    random_forest_model <- randomForest(genotype ~., data = rf_training, ntree = i)
	    rf_training_table <- table(actual = rf_training$genotype,
                       predicted = predict(random_forest_model, type="class"))
	    rf_test_table <- table(actual = rf_test$genotype,
                       predicted = predict(random_forest_model, newdata=rf_test, type="class"))
	    prediction_success <- 100 * sum(diag(rf_test_table)/sum(rf_test_table))
		results_list <- list(unlist(results_list), prediction_success)
		}
	return(results_list)
}

#re-sample the 10 trees, 1000 trees random forest comparison functions
resample_comparison <- function(function_name, dataset, datatable, reps){
	out_list1 <- list()
	out_list2 <- list()
    for(i in 1:reps)
	{
        myresult <- function_name(dataset, datatable)
		out_list1 <- unlist(c(out_list1, myresult[[1]]))
		out_list2 <- unlist(c(out_list2, myresult[[2]]))
     }
	mean_score1 <- mean(out_list1)
	mean_score2 <- mean(out_list2)
	error_score1 <- sd(out_list1)
	error_score2 <- sd(out_list2)
	return_list <- list(c(mean_score1, error_score1), out_list1, c(mean_score2,error_score2), out_list2) 
	return(return_list)
}


#left SAM wings, for sex comparisons
#no centroid
rf_sex_info_left_sam <- resample_comparison(rf_sex_comparison, left_sam_wings_data, tabls, 5) 
rf_sex_info_left_sam #first result is 10 trees, percentage and error, second result is 1000 trees, percentage and error
#with centroid
rf_sex_info_left_sam_centroid <- resample_comparison(rf_sex_comparison, left_sam_wings_cent, tablsc, 1000)
rf_sex_info_left_sam_centroid

#left female wings, for genotype comparisons
#no centroid
rf_genotype_info_left <- resample_comparison(rf_genotype_comparison, left_female_wings_data, tablw_g, 5) 
rf_genotype_info_left #first result is 10 trees, percentage and error, second result is 1000 trees, percentage and error
#with centroid
rf_genotype_info_left_centroid <- resample_comparison(rf_genotype_comparison, left_female_wings_cent, tablwc_g, 5) 
rf_genotype_info_left_centroid  




#Code for figures
#########################################################################################
#3. WingBlur function - WP 2014
#########################################################################################
WingBlur <- function ( PrcCrds, wingcol="#add8e632", winglwd=3, winglty=1, wingpch=1, wingcex=0.9, meanline=TRUE, meancol="blue", meanlty=1, meanlwd=2 ) {
  # set up required coord vectors
  from <- c( 6, 47, 5, 29, 5, 28, 27, 26, 25, 4, 24, 23, 3, 22, 2, 21, 20, 19, 1, 18, 17, 16, 15, 14, 12, 33, 32, 7, 31, 8, 48, 37, 9, 36, 8, 35, 34, 30, 11, 42, 10, 41, 40, 39, 38, 10, 11, 46, 45, 44, 43 )
  to <- c( 47, 5, 29, 30, 28, 27, 26, 25, 4, 24, 23, 3, 22, 2, 21, 20, 19, 1, 18, 17, 16, 15, 14, 13, 33, 32, 7, 31, 1, 48, 7, 9, 36, 8, 35, 34, 2, 11, 42, 10, 41, 40, 39, 38, 3, 9, 46, 45, 44, 43, 4 )
  x.vals <- seq(1, 95, 2)		#these pull out the odd...
  y.vals <- seq(2, 96, 2)		#...and even columns...
  shape.X <- PrcCrds[,x.vals]	# make x coord matrix
  shape.Y <- PrcCrds[,y.vals]	# and one for y
  LMx.mean <- colMeans(PrcCrds)[x.vals]	# make mean x vector
  LMy.mean <- colMeans(PrcCrds)[y.vals]	# and one for y

  # set up the plot
  plot( LMx.mean, LMy.mean, xlim=c(-0.22,0.27), ylim=c(-0.22,0.27), type="n", xaxt="n", yaxt="n", ann=FALSE)

  # nested loops to draw each line segment for each row in shape.X matrices
  for (j in 1:nrow(shape.X)){
    for (i in 1:length(from)) {
      lines( c(shape.X[j,from[i]], shape.X[j,to[i]]), c(shape.Y[j,from[i]], shape.Y[j,to[i]]), col=wingcol, lty=winglty, lwd=3 )
    }	}

  # loop to draw the mean lines if requested
  if (meanline ==TRUE) {
    for (i in 1:length(from)) {
      lines( c(LMx.mean[from[i]], LMx.mean[to[i]]), c(LMy.mean[from[i]], LMy.mean[to[i]]), col=meancol, lty=meanlty, lwd=meanlwd )
    }	}

  # finally, this drop in the mean LM/semi-LM points
  points(LMx.mean, LMy.mean, cex=wingcex, pch=wingpch)
}

# wingblur variation plot
par( mfrow= c(1, 1))
WingBlur( wings[, 9:104] )
text( 0, 0.25, "Variation among 4x" )
text( 0, 0.23, "Olympus wing images" )
text( 0, -0.2, "scale factor = 1")


#4. Wing Effect plots
#########################################################################################
WingEffect <- function ( meanshape, effectplus, effectminus, wingcol=c("black", "black", "red"), winglwd=c(2, 2, 2), winglty=c(1, 1, 1), scale.factor=2, meanline=FALSE, wingcex=c(0,0,0), winglabel="", add=FALSE, wingframe=T, scale.display=T, wingpoints=F, wingpch=c(16,1,1) ) {

  # set up required coord vectors
  from <- c( 6, 47, 5, 29, 5, 28, 27, 26, 25, 4, 24, 23, 3, 22, 2, 21, 20, 19, 1, 18, 17, 16, 15, 14, 12, 33, 32, 7, 31, 8, 48, 37, 9, 36, 8, 35, 34, 30, 11, 42, 10, 41, 40, 39, 38, 10, 11, 46, 45, 44, 43 )
  to <- c( 47, 5, 29, 30, 28, 27, 26, 25, 4, 24, 23, 3, 22, 2, 21, 20, 19, 1, 18, 17, 16, 15, 14, 13, 33, 32, 7, 31, 1, 48, 7, 9, 36, 8, 35, 34, 2, 11, 42, 10, 41, 40, 39, 38, 3, 9, 46, 45, 44, 43, 4 )
  x.vals <- seq(1, 95, 2)	# indices to specify x and y coords
  y.vals <- seq(2, 96, 2)
  LMx.mean <- meanshape[x.vals]	# apply x & y indices to mean data
  LMy.mean <- meanshape[y.vals]
  # x & y indices to effectplus vector
  x.plus <- LMx.mean + (effectplus[x.vals] * scale.factor)
  y.plus <- LMy.mean + (effectplus[y.vals] * scale.factor)
  # x & y indices to effectminus vector
  x.minus <- LMx.mean - (effectminus[x.vals] * scale.factor)
  y.minus <- LMy.mean - (effectminus[y.vals] * scale.factor)

  # plot mean shape
  if ( add == FALSE ) { # opens blank plot
    plot( LMx.mean, LMy.mean, xlim=c(-0.22,0.27), ylim=c(-0.22,0.27), type="n", xaxt="n", yaxt="n", ann=FALSE, frame=wingframe)
  }
  if (meanline == TRUE){
    # this for loop draws line segments for the mean shape
    for (i in 1:length(from)) {
      lines( c(LMx.mean[from[i]], LMx.mean[to[i]]), c(LMy.mean[from[i]], LMy.mean[to[i]]), col=wingcol[1], lty=winglty[1], lwd=winglwd[2])
    }	}
  if (wingpoints == TRUE) {
    points( x.plus, y.plus, cex=wingcex[2], pch=wingpch[2], col=wingcol[2])
    points( x.minus, y.minus, cex=wingcex[3], pch=wingpch[3], col=wingcol[3])
  }

  # this drop in the mean LM/semi-LM points
  points(LMx.mean, LMy.mean, cex=wingcex[1], pch=wingpch[1], col=wingcol[1])

  # this for loop draws the mean+effect line segments
  for (i in 1:length(from)) {
    lines( c(x.plus[from[i]], x.plus[to[i]]), c(y.plus[from[i]], y.plus[to[i]]), col=wingcol[2], lty=winglty[2], lwd=winglwd[2] )
    lines( c(x.minus[from[i]], x.minus[to[i]]), c(y.minus[from[i]], y.minus[to[i]]), col=wingcol[3], lty=winglty[3], lwd=winglwd[3] )
  }
  # adding an annotation of the effect scaling factor
  if (scale.display == T) {
    text( 0.025, -0.215, paste("scale factor =", as.character(scale.factor)) )
  }
  text(0, 0.25, winglabel) # optional centre-top label
}


# calculate tangent approximation for tangent approximates Procrustes Distance (Euclidean Distance)
PD <- function(x) {
  sqrt(t(x)%*%x)}


wings_id <- aggregate( wings[, 7:104], by=list( genotype=wings$genotype, sex=wings$sex, microscope=wings$microscope, magnification=wings$magnification, wings$ID ), mean )[,-5]

procoords <- wings_id[, 7:102] * matrix( rep( rep( c(-1, 1), 48), 1140), ncol=96, byrow = T)

meanshape <- colMeans( procoords )

levels( wings_id$genotype )

egfr_diff <- colMeans( procoords[ wings_id$genotype == "samw", ] ) - colMeans( procoords[ wings_id$genotype == "egfr", ] )
mam_diff <- colMeans( procoords[ wings_id$genotype == "samw", ] ) - colMeans( procoords[ wings_id$genotype == "mam", ] )
star_diff <- colMeans( procoords[ wings_id$genotype == "samw", ] ) - colMeans( procoords[ wings_id$genotype == "star", ] )
tkv_diff <- colMeans( procoords[ wings_id$genotype == "samw", ] ) - colMeans( procoords[ wings_id$genotype == "tkv", ] )

##### Figure 2
# mean effects of mutation plots
par( mfrow= c(2, 2), mar=c(1,1,1,1))
WingEffect( meanshape, egfr_diff, egfr_diff, scale.factor=3, winglabel="Egfr"  )
text( 0, 0.2, paste( "PD=", round(PD(egfr_diff), 4)) )
WingEffect( meanshape, mam_diff, mam_diff, scale.factor=3, winglabel="mam" )
text( 0, 0.2, paste( "PD=", round(PD(mam_diff), 4)) )
WingEffect( meanshape, star_diff, star_diff, scale.factor=3, winglabel="S"  )
text( 0, 0.2, paste( "PD=", round(PD(star_diff), 4)) )
WingEffect( meanshape, tkv_diff, tkv_diff, scale.factor=3, winglabel="tkv" )
text( 0, 0.2, paste( "PD=", round(PD(tkv_diff), 4)) )



#5. LDA plots
#########################################################################################
#plotting LDA results , from training set only

#single LDA on all wings, no repetitions
#training and test set
wings_genotype_train_test <- strata_5var(average_wings, "genotype", tabw_g, 2/3, 2, 158) #variable position and table length are different for the averaged dataset
wings_wg_training <- wings_genotype_train_test[[1]]
wings_wg_test <- wings_genotype_train_test[[2]]
#model
linDiscrim_all <- lda(formula=genotype ~.,data = wings_wg_training, tol = 1.0e-8, CV=FALSE) #modified tolerance
all_lda_training_table <- table(actual = wings_wg_training$genotype,
                                predicted = predict(linDiscrim_all, newdata=wings_wg_training)$class)
all_lda_training_table
all_lda_test_table <- table(actual = wings_wg_test$genotype,
                            predicted = predict(linDiscrim_all, newdata=wings_wg_test)$class)
all_lda_test_table
#accuracy rate
100*sum(diag(all_lda_test_table)/sum(all_lda_test_table))


aLDA <-  data.frame(wings_wg_training, as.matrix(wings_wg_training[,2:59]) %*% linDiscrim_all$scaling ) #matrix-multiplying the scaling matrix by the LM's gives a vector of individual LD 'scores' for each LDA
genocols <- c( "#00008B99", "#B8860B99", "#00000099", "#00640099", "#8B000099")
genopch <- c( 21, 25, 23, 24, 22 )
par( mfrow= c(1, 1), mar=c(5,5,1,1))
plot( aLDA$LD1, aLDA$LD2, pch=genopch[aLDA$genotype], bg=genocols[aLDA$genotype], col=genocols[aLDA$genotype], xlab = "1st Discriminant Function", ylab = "2nd Discriminant Function", cex=1.5 )
# x = lda 1, y = lda2
legend( -7, 4.7, bty = "n", pch = genopch, pt.bg=c("#00008B", "#B8860B", "#000000", "#006400", "#8B0000"), col=c("#00008B", "#B8860B", "#000000", "#006400", "#8B0000"), c("Egfr", "mam", "sam", "star", "tkv"), xjust=0 )


#plotting LDA results, test set only
bLDA <-  data.frame(wings_wg_test, as.matrix(wings_wg_test[,2:59]) %*% linDiscrim_all$scaling ) #matrix-multiplying the scaling matrix by the LM's gives a vector of individual LD 'scores' for each LDA
genocols <- c( "#00008B99", "#B8860B99", "#00000099", "#00640099", "#8B000099")
genopch <- c( 21, 25, 23, 24, 22 )
par( mfrow= c(1, 1), mar=c(5,5,1,1))
plot( bLDA$LD1, bLDA$LD2, pch=genopch[bLDA$genotype], col=genocols[bLDA$genotype], bg=genocols[bLDA$genotype], xlab = "1st Discriminant Function", ylab = "2nd Discriminant Function", cex=1.5 ) # x = lda 1, y = lda2

legend(3.6, 4.9, bty = "n", pch = genopch, pt.bg=c("#00008B", "#B8860B", "#000000", "#006400", "#8B0000"), col=c("#00008B", "#B8860B", "#000000", "#006400", "#8B0000"), c("Egfr", "mam", "sam", "star", "tkv"), xjust=0)

 
#6. Heat map plots of confusion matrices
#########################################################################################

myHeatMap <- function( dim, names, data, xLab="Actual", yLab="Predicted", Main="" ){
	myCols <- colorRampPalette(c("#ffffff7f", "#ffff007f", "#ff00007f"), space="rgb")
	image( 1:dim, 1:dim, data[,dim:1], col=myCols(100), xaxt='n', yaxt='n', xlab=xLab, ylab=yLab, main=Main )
	abline(h=0:dim+0.5, col="grey")
	abline(v=0:dim+0.5, col="grey")
	box( lwd=1 )
	text( 1:dim, rep(dim:1, each=dim), sub('^0$', '',data) )
	axis(1, at=1:dim, labels=names, cex.axis=0.8, font=3 )
	axis(2, at=1:dim, labels=rev(names), cex.axis=0.8, las=1, font=3 )
}



#biocat information -- see supplementary methods for raw data
Biocat_RF_10_confusion <- matrix(c(8,0,5,0,17,0,8,0,22,0,1,4,23,0,2,0,4,0,26,0,7,5,5,0,13),ncol=5, byrow = TRUE)
#this data came from my Biocat Notes, see same folder
# 
Biocat_RF_10_conf_perc <- (Biocat_RF_10_confusion / rowSums( Biocat_RF_10_confusion ) )*100



### landmark data
wings_info <- rf_genotype_comparison( left_female_wings_data, tablwc_g)[[3]]

wings_info_perc <- 100*( wings_info / rowSums( wings_info ) )


### combining the two
wings_info_combined <- wings_info_perc
rownames(wings_info_combined) <- c("wm_egfr", "wm_mam", "wm_sam", "wm_star", "wm_tkv")
Biocat_info_combined <- Biocat_RF_10_conf_perc
rownames(Biocat_info_combined ) <- c("biocat_egfr", "biocat_mam", "biocat_sam", "biocat_star", "biocat_tkv")
RF_10_combined <- rbind(wings_info_combined, Biocat_info_combined)
RF_10_combined


### code-block for Figure6

layout( matrix( c(1,3,2,3), 2, byrow=T), T)

myHeatMap( dim=5, names=c("Egfr","mam","sam","Star","tkv"), data=round(Biocat_RF_10_conf_perc, 0), Main="Biocat" )
	mtext( "(a)", 1, 1, at=-.5 )

myHeatMap( dim=5, names=c("Egfr","mam","sam","Star","tkv"), data=round(wings_info_perc, 0), Main="Landmarks" )
	mtext( "(b)", 1, 1, at=-.5 )

myCols <- colorRampPalette(c("#ffffff7f", "#ffff007f", "#ff00007f"), space="rgb")

rotate = function(mat){ t(mat[nrow(mat):1,,drop=FALSE]) }
RF_10_combined_rot <- round( rotate( RF_10_combined ), 0)

image( 1:5, 1:10, RF_10_combined_rot[,10:1], col=myCols(100), xaxt='n', yaxt='n', xlab='Actual', ylab='' )
		abline(h=0:10+0.5, col="grey")
		abline(v=0:5+0.5, col="grey")
		box( lwd=1 )
		text( 1:5, rep(10:1, each=5), sub('^0$', '', RF_10_combined_rot) )
	axis(1, at=1:5, labels=c("Egfr","mam","sam","Star","tkv"), cex.axis=0.8, font=3 )
	axis(2, at=1:10, labels=c( c("lm_Egfr", "lm_mam", "lm_sam", "lm_Star", "lm_tkv"), c("biocat_Egfr", "biocat_mam", "biocat_sam", "biocat_Star", "biocat_tkv")), cex.axis=0.8, las=1, font=3 )
	mtext( "(c)", 1, 1, at=-.5 )
	



