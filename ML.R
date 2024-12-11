library(dplyr)
library(plyr)
library(forecast)
library(randomForest)
library(xts)
library(ggplot2)
library(depmixS4)
#selecting training vs testing data
trainData <- pVisi_data[2:407]
testData <- pVisi_data[407:544] 
#ARIMA forecasting 
data1 <- pVisi_data[[1]]
trainData <- data1$SPO2[1:(length(data1$SPO2)*0.75)]
testData <- data1$SPO2[(length(data1$SPO2)*0.75):(length(data1$SPO2))]
trainData <- ts(trainData[(length(trainData)*0.5):(length(trainData)*0.75)], frequency = 4)
model2 <- auto.arima(trainData, d = 2, trace = T)
plot.ts(model2$residuals)
myforecast <- forecast(model2, level=c(95), h=50)
mypredict <- predict(model2, 50)
prediction <- c(trainData, mypredict$pred)
plot(prediction)
plot(myforecast)
#HMM Model for G = 3 ----
g3Patients <- rownames(gaussianG)[gaussianG$G == 3] #filtering all the G3 patients
summaryg3 <- as.data.frame(matrix(data = NA, nrow = length(g3Patients), ncol = 6))
transitionMatrix <- as.data.frame(matrix(data = NA, nrow = length(g3Patients), ncol = 9))
summaryg3_2 <- as.data.frame(matrix(data = NA, nrow = length(g3Patients), ncol = 6))
transitionMatrix_3 <- as.data.frame(matrix(data = NA, nrow = length(g3Patients), ncol = 9))
colnames(summaryg3_2) <- c("St1Intercept", "St2Intercept", "St3Intercept", "St1SD", "St2SD", "St3SD")
for (x in 253:length(g3Patients)) {
  j <- as.numeric(g3Patients[x])
  print(j)
  dataH <- as.data.frame(gaussianData[[j]])
  dataH$SYS <- as.numeric(pVisi_data[[j]][,4])
  dataH$DIA <- as.numeric(pVisi_data[[j]][,5])
  covariates <- baseline[baseline$randomization_number == names(pVisi_data[j]),c(1,61)]
  dataH$age <- rep(covariates$age, length.out = nrow(dataH))
  dataH$bmi <- covariates$bmi
  #HModel <- depmix(dataH$V1~1, nstates = gaussianG$G[j], ntimes = length(dataH$V1), transition = ~ pVisi_data[[j]][,2] , family = gaussian(), prior = ~ age, initdata = dataH[3,])
  dataH <- na.omit(dataH)
  #HModel <- depmix(dataH$V1~1, nstates = gaussianG$G[j], ntimes = length(dataH$V1), transition = ~ SYS + DIA, data = dataH)
  HModel <- depmix(dataH$V1~1, nstates = gaussianG$G[j], ntimes = length(dataH$V1))
  HFit <- fit(HModel, verbose = FALSE,  )
  #summaryg3[x,] <- c(summary(HFit)[,1], summary(HFit)[,2])
  #transitionMatrix[x,] <- HFit@trDens
  summaryg3_2[x,] <- c(summary(HFit)[,1], summary(HFit)[,2])
  transitionMatrix_3[x,] <- HFit@trDens
  #ggplot(aes(x=1:nrow(dataH), y=V1, colour=as.factor(state)), data=dataH)+ geom_point() + labs(title = "Patient 1", x = "time", y = "SPO2")  
}
for (x in 2:length(g3Patients)){
  o <- order(summaryg3_2[x,1:3])
  summaryg3[x,] <- summaryg3_2[x,c(o,o+3)]
  transitionMatrix[x,] <- transitionMatrix_3[x,unlist(lapply(o, function(x) (x - 1) * 3 + 1:3))]
}
summaryg3$ID <- as.numeric(gaussianG$ID[as.numeric(g3Patients)])
colnames(transitionMatrix) <- c("onetoone","1->2","1->3", "2->1","twototwo","2->3","3->1","3->2","threetothree")
summaryg3 <- na.omit(summaryg3)
transitionMatrix <- na.omit(transitionMatrix)
summaryg3 <- left_join(summaryg3, baseline, by = c("ID" = "randomization_number"))
summaryg3$gender <- factor(summaryg3$gender, levels = c(0, 1), labels = c("male", "female"))
summaryg3 <- summaryg3 %>% filter(!is.na(gender)) 
summaryg3$onetoone <- transitionMatrix$onetoone
summaryg3$twototwo <- transitionMatrix$twototwo
summaryg3$threetothree <- transitionMatrix$threetothree
summaryg3$smoking_status <- ifelse(summaryg3$smoking_status == 0, "non-smoking", "smoking")
summaryg3$bmi_category <- ifelse(summaryg3$bmi >= 25, "BMI >= 25", "BMI < 25")
summaryg3$age_category <- ifelse(summaryg3$age >= 50, "age >= 50", "age < 50")
colnames(summaryg3)[1:6] <- c("St1Intercept", "St2Intercept", "St3Intercept", "St1SD", "St2SD", "St3SD")
#HMM Model for G = 2 ----
g2Patients <- rownames(gaussianG)[gaussianG$G == 2] #filtering all the G3 patients
summaryg2 <- as.data.frame(matrix(data = NA, nrow = length(g2Patients), ncol = 4))
transitionMatrixg2 <- as.data.frame(matrix(data = NA, nrow = length(g2Patients), ncol = 4))
colnames(summaryg2) <- c("St1Intercept", "St2Intercept", "St1SD", "St2SD")
forecasts <- as.data.frame(matrix(data = NA, nrow = length(g2Patients), ncol = 10))
for (x in 2:length(g2Patients)) {
  j <- as.numeric(g2Patients[x])
  print(x)
  dataH <- as.data.frame(gaussianData[[j]])
  #SIMPLE MODEL
    HModel <- depmix(dataH$V1~1, nstates = gaussianG$G[j], ntimes = length(dataH$V1))
  #ADDING SYSTOLIC and DIASTOLIC BP
    # dataH$SYS <- as.numeric(pVisi_data[[j]][,4])
    # dataH$DIA <- as.numeric(pVisi_data[[j]][,5])
    # dataH <- na.omit(dataH)
    # HModel <- depmix(dataH$V1~1, nstates = gaussianG$G[j], ntimes = length(dataH$V1), transition = ~ SYS + DIA, data = dataH)
  #ADDING PRIOR
    # dataH <- na.omit(dataH)
    # dataH$age <- baseline[baseline$randomization_number == names(pVisi_data[j]),c(1,61)]$age
    # dataH$bmi <- baseline[baseline$randomization_number == names(pVisi_data[j]),c(1,61)]$bmi
    # initialData <- as.data.frame(dataH[,2:3])
    # HModel <- mix(V1~1, nstates = gaussianG$G[j], data = dataH, prior=~age+bmi,initdata = dataH, family = multinomial("identity"))
  HFit <- fit(HModel, sverbose = FALSE,  )
  summaryg2[x,] <- c(summary(HFit)[,1], summary(HFit)[,2])
  transitionMatrixg2[x,] <- HFit@trDens
  #making future predictions
  # transition_mat <- rbind(getpars(getmodel(HFit, "transition", 1)),getpars(getmodel(HFit, "transition", 2)))
  # priorVec <- as.numeric(posterior(HFit)[length(dataH$V1), -1])
  # predRState <- c(getpars(getmodel(HFit, "response", 1))[1], getpars(getmodel(HFit, "response", 2))[1])
  # for (i in 1:ncol(forecasts)) {
  #   priorVec <- priorVec %*% transition_mat
  #   forecasts[x,i] <- sum(predRState * priorVec)
  # }
  #ggplot(aes(x=1:nrow(dataH), y=V1, colour=as.factor(state)), data=dataH)+ geom_point() + labs(title = "Patient 1", x = "time", y = "SPO2")  
}
#ordering the summary matrix by Mean 1 being smallest, and mean 2 being largest
#transitionMatrixg2 <- transitionMatrixg2[2:length(transitionMatrixg2),]
swappingRows <- which(summaryg2$St1Intercept > summaryg2$St2Intercept) #now S2 is larger
summaryg2[swappingRows,] <- summaryg2[swappingRows, c(2,1,4,3)]
transitionMatrixg2[swappingRows,] <- transitionMatrixg2[swappingRows, c(4,3,2,1)]
summaryg2$ID <- as.numeric(gaussianG$ID[as.numeric(g2Patients)])
summaryg2 <- left_join(summaryg2, baseline, by = c("ID" = "randomization_number"))
summaryg2$gender <- factor(summaryg2$gender, levels = c(0, 1), labels = c("male", "female"))
summaryg2 <- summaryg2 %>% filter(!is.na(gender)) 
summaryg2$onetoone <- transitionMatrixg2$V1[-1]
summaryg2$twototwo <- transitionMatrixg2$V4[-1]
summaryg2$diff2minus1 <- summaryg2$St2Intercept - summaryg2$St1Intercept 
summaryg2$diff22minus11 <- summaryg2$twototwo - summaryg2$onetoone
summaryg2$smoking_status <- ifelse(summaryg2$smoking_status == 0, "non-smoking", "smoking")
summaryg2$bmi_category <- ifelse(summaryg2$bmi >= 30, "BMI >= 30", "BMI < 30")
summaryg2$age_category <- ifelse(summaryg2$age >= median(summaryg2$age), "age >= 50", "age < 50")
# forecasts[1,] <- 1:10
# forecasts <- as.data.frame(t(forecasts))
# #plotting forecasts
# time <- (length(dataH$V1))+(1:10)
# plot(time, forecasts[,249], type = "l", col = "blue", xlab = "Time", ylab = "Forecasted Values", main = "HMM Forecasts")
# lines(1:length(dataH$V1), dataH$V1, col = "red", lty = 2)

#Random Forest forecasting
data1 <- pVisi_data[[1]]
lagged_data <- lag(data1$SPO2)  # Create lagged data

data("AirPassengers")
ts_data <- AirPassengers
ts_df <- data.frame(Date = index(ts_data), Passengers = coredata(ts_data))
ts_df$Date <- as.Date(ts_df$Date)
ts_xts <- xts(ts_df$Passengers, order.by = ts_df$Date)
