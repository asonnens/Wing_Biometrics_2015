Scripts:

Wing_Biometrics_2014.r
	Script for analyzing landmark and semi-landmark coordinate data from the Wing Biometrics Database.

Wing_Biometrics_settings.r
	Script used for optimizing settings on 'svm', 'knn', and 'randomForest'
	These parameters were applied in Wing_Biometrics_2014.r

These scripts were both written for R version 3.1.0

These scripts are intended to run on the file wingBiometry_oly_4X_coords.tsv, in the data subfolder.

Although the information in this file is identical to the information in the file Olympus_4X_coords.tsv, 
all files in the Landmark_Semilandmark_Coordinates folder have an additional column indicating the dimensions 
of the splined images.

Before applying these scripts to other datafiles, check column paramters-- in many cases, they will need to be adjusted by one.