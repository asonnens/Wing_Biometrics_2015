Data:
The five data files are all in tab seperated format. They can be opened in any text editor, and should easily transfer to spreadsheets. 
All five files are parsed output from CPR version OI.
The files Leica_2X_coords.tsv, Leica_4X_coords.tsv, Olympus_2X_coords.tsv, and Olympus_4X_coords.tsv were all superimposed separately in CPR.
all_wings.tsv represents output from all four datasets, jointly superimposed in CPR. 

Scripts:

Wing_Biometrics_2014.r
	Script for analyzing landmark and semi-landmark coordinate data from the Wing Biometrics Database.

Wing_Biometrics_settings.r
	Script used for optimizing settings on 'svm', 'knn', and 'randomForest'
	These parameters were applied in Wing_Biometrics_2014.r

These scripts were both written for R version 3.1.0
R packages used include:
plyr package (V.1.8.1)
MASS package (V. 7.3-33)
adabag package (V. 3.2)
randomForest package (V. 4.6-7)
e1071 package (V. 1.6-3),
nnet package (V. 7.3-8)
class package (V. 7.3-10)

These scripts are intended to run on tsv files, in the data subfolder.

Many aspects of these scripts were based off of the training scripts in  Dworkin_machine_learing_in_R_tutorials.
