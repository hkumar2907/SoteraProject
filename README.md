# SoteraProject
Sotera project at Cleveland Clinic

File descriptions:
1. visiData.R: Loading data, calculating the primary outcomes from the reference paper, and getting Gaussian distributions for patients
2. SoteraAnalysis.qmd: figures code from Rowan, some figures comparing Hidden Markov model transition values for different clinical variables
3. ML.R: fitting the hidden Markov model to Gaussian means, plots comparing hidden Markov model transition mean differences
4. missingData.R: plotting the missing data in the original dataset
5. soteraPlots.R: some random plot; not important lol
6. windowsAnalysis.R: generating the windows for the SpO2 (different codes for it, with different overlaps. calculating the summary features for each window, plotting sample SpO2 traces for each window, feature correlation, and plotting the association between features, subsetting the different windows into stable/deteriorating, etc,. and associated numbers/plots
7. randomForest.R: combining the different window types into the overall dataset, train/test partition for the dataset, tuning the random forest parameters, running the random forest model, testing the random forest model, calculating random forest model evaluations, calculating alaram burden, false positive windows calculations
8. aggregateHAR.R: figures for HAR not sotera paper
