#Libraries ----
library(dplyr)
library(zoo)
library(tidyr)
library(rstatix)
library(ggplot2)
library(tidyverse)
#library(GGally)
library(mclust)
library(viridis) 
library(datawizard)
library(patchwork)
#Loading data ----
pOJama <- read.csv("primary_outcomes_jama_revision.csv") #loading their primary outcome data 
finalStudyPop <- read.csv("Final_study_pop_20190107.csv") #loading final study population data, start date and end date on this file doesn't make sense
visi_data <- read.csv("all_factor_visi_data.csv")
colnames(visi_data)[13] <- "StudyID" #changing the last column name to just StudyID
visi_data$StudyID <- gsub("-", "", visi_data$StudyID) #removing dashes as that it how it is in the final excel
visi_data[visi_data==''] <- NA
pVisi_data <- split(visi_data, visi_data$StudyID) #splitting into a list sorted by ID
baseline <- read.csv("baseline_20190121.csv")
intraop <- read.csv("intraop_20190121.csv")
intraop$pacu_arrival <- strptime(intraop$pacu_arrival, "%d%b%Y:%H:%M:%S")
intraop$pacu_discharge_eligible <- strptime(intraop$pacu_discharge_eligible, "%d%b%Y:%H:%M:%S")
intraop$pacu_discharge <- strptime(intraop$pacu_discharge, "%d%b%Y:%H:%M:%S")
intraop$pacuEligibleHours <- as.numeric(difftime(intraop$pacu_discharge_eligible, intraop$pacu_arrival, units = "hours"))
intraop$Duration.of.surgery..hours. <- as.numeric(intraop$Duration.of.surgery..hours.)
intraop$pacu_los_hr <- as.numeric(intraop$pacu_los_hr)
#Calculating the primary outcomes ----
primaryOutcomes <- as.data.frame(matrix(ncol = 6, nrow = length(pVisi_data)))
colnames(primaryOutcomes) <- c("StudyID", "Total hypoxia (hour)", "Hypoxia (min per hour)", "Total reading time" ,"Total Gap", "AUC below 90")
for (j in 1:length(pVisi_data)) {
  data1 <- pVisi_data[[j]]
  #removing all XXs
  data1 <- data1[data1$SPO2 != "XX",] #removing all XXs in the SPO2
  data1$SPO2 <- as.numeric(data1$SPO2) #turning spo2 is numeric
  data1 <- data1 %>% drop_na(SPO2)
  data1 <- data1[order(data1$dev_reading_tm, decreasing = FALSE),]
  data1$dev_reading_tm <- strptime(data1$dev_reading_tm, "%d%b%Y:%H:%M:%S")
  hypox_dura_min <- 0 #"As a sensitivity analysis, the total duration of hypoxemia over a 48-hour period (excluding gaps) was also calculated."
  totalgap <- 0 #number of missing datapoints
  totalReading <- 0
  print(j)
  for (i in 2:nrow(data1)) {
    difft <- as.numeric(difftime(data1$dev_reading_tm[i], data1$dev_reading_tm[i-1]))
    totalReading <- totalReading + difft
    if (difft > 60) { 
      totalgap = totalgap + difft
    }
    if (data1$SPO2[i] < 90) { #if spo2 below 90, then a time period added to it
      hypox_dura_min <- hypox_dura_min + difft
    } 
    # else if (dewqata1$SPO2[i-1] < 90) { #if current spo2 not below 90, but previos below then a time period added to it
    #   hypox_dura_min <- hypox_dura_min + difft
    # }
  }
  pVisi_data[[j]] <- data1
  totalReading <- totalReading/60 #converting seconds to minutes
  hypox_dura_min <- hypox_dura_min/60 #converting seconds to minutes #equal to hypox_dura_min
  totalgap <- totalgap/60 #converting seconds to minutes
  primaryOutcomes$StudyID[j] <- names(pVisi_data)[[j]]
  primaryOutcomes$`Hypoxia (min per hour)`[j] <- hypox_dura_min/(totalReading - totalgap) *60 #was calculated as (totalminutes of SpO2<90%)/(total SpO2 readingminutes − total gap)
  primaryOutcomes$`Total reading time`[j] <- totalReading/60 #converting minutes to hours #equal to total_duration_hr
  primaryOutcomes$`Total Gap`[j] <- totalgap
  primaryOutcomes$`Total hypoxia (hour)`[j] <- hypox_dura_min #converting hour to minute
  #data1time <- as.numeric(difftime(data1$dev_reading_tm, min(data1$dev_reading_tm), units = "mins"))
  #caTools::trapz(data1time[which(data1$SPO2 < 90)], (90 - data1$SPO2[which(data1$SPO2 < 90)]))
  vals <- (90 - data1$SPO2)
  vals[vals < 0] <- 0
  primaryOutcomes$`AUC below 90`[j] <- caTools::trapz(as.numeric(difftime(data1$dev_reading_tm, min(data1$dev_reading_tm), units = "mins")), vals)
}
distribution <- as.data.frame(matrix(ncol = 21, nrow = length(pVisi_data)))
colnames(distribution) <- c("StudyID", "Min", "1st Q", "Median", "Mean", "3rd Q", "Max", "B_Min", "B_1st Q", "B_Median", "B_Mean", "B_3rd Q", "B_Max","A_Min", "A_1st Q", "A_Median", "A_Mean", "A_3rd Q", "A_Max", "Lower whisker", "Upper whisker")
distribution$StudyID <- names(pVisi_data)
#getting distribution below 90 and above 90 ----
for (j in 1:length(pVisi_data)) {
  print(j)
  data1 <- pVisi_data[[j]]
  distribution[j, 2:7] <- summary(data1$SPO2)
  distribution[j, 8:13] <- summary(data1[data1$SPO2 < 90, ]$SPO2)
  distribution[j, 14:19] <- summary(data1[data1$SPO2 >= 90, ]$SPO2)
  distribution[j, 20] <- boxplot.stats(data1$SPO2)$stats[1]
  distribution[j, 21] <- boxplot.stats(data1$SPO2)$stats[5]
}
distribution$AMedianMinusBMedian <- distribution$A_Median - distribution$B_Median
for (j in 1:length(pVisi_data)) {
  print(j)
  print(which(!is.na(pVisi_data[[j]]$BP_DIA)))
}
#PLOTTING
id <- "376010035"
data1 <- pVisi_data[[id]]
p1 <- boxplot(data1$SPO2, data1[data1$SPO2 < 90, ]$SPO2, data1[data1$SPO2 >= 90, ]$SPO2, main = "G = 1", names = c("Overall", "Below 90", "Above 90"), horizontal = TRUE, xlab = "SPO2")
id <- "376010009"
data1 <- pVisi_data[[id]]
p2 <- boxplot(data1$SPO2, data1[data1$SPO2 < 90, ]$SPO2, data1[data1$SPO2 >= 90, ]$SPO2, main = "G = 2", names = c("Overall", "Below 90", "Above 90"), horizontal = TRUE, xlab = "SPO2")
id <- "376010001"
data1 <- pVisi_data[[id]]
p3 <- boxplot(data1$SPO2, data1[data1$SPO2 < 90, ]$SPO2, data1[data1$SPO2 >= 90, ]$SPO2, main = "G = 3", names = c("Overall", "Below 90", "Above 90"), horizontal = TRUE, xlab = "SPO2")
p1+p2+p3
#identify_outliers(data1$SPO2)

#setting hypoxia threshold to their lower quartile bar ----
primaryOutcomesLowWhisker <- as.data.frame(matrix(ncol = 6, nrow = length(pVisi_data)))
colnames(primaryOutcomesLowWhisker) <- c("StudyID", "Total hypoxia (hour)", "Hypoxia (min per hour)", "Total reading time" ,"Total Gap", "AUC below 90")
for (j in 1:length(pVisi_data)) {
  data1 <- pVisi_data[[j]]
  #removing all XXs
  hypox_dura_min <- 0 #"As a sensitivity analysis, the total duration of hypoxemia over a 48-hour period (excluding gaps) was also calculated."
  totalgap <- 0 #number of missing datapoints
  totalReading <- 0
  print(j)
  for (i in 2:nrow(data1)) {
    difft <- as.numeric(difftime(data1$dev_reading_tm[i], data1$dev_reading_tm[i-1]))
    totalReading <- totalReading + difft
    if (difft > 60) { 
      totalgap = totalgap + difft
    }
    if (data1$SPO2[i] < distribution$`Lower whisker`[j]) { #if spo2 below 90, then a time period added to it
      hypox_dura_min <- hypox_dura_min + difft
    } 
  }
  pVisi_data[[j]] <- data1
  totalReading <- totalReading/60 #converting seconds to minutes
  hypox_dura_min <- hypox_dura_min/60 #converting seconds to minutes #equal to hypox_dura_min
  totalgap <- totalgap/60 #converting seconds to minutes
  primaryOutcomesLowWhisker$StudyID[j] <- names(pVisi_data)[[j]]
  primaryOutcomesLowWhisker$`Hypoxia (min per hour)`[j] <- hypox_dura_min/(totalReading - totalgap) *60 #was calculated as (totalminutes of SpO2<90%)/(total SpO2 readingminutes − total gap)
  primaryOutcomesLowWhisker$`Total reading time`[j] <- totalReading/60 #converting minutes to hours #equal to total_duration_hr
  primaryOutcomesLowWhisker$`Total Gap`[j] <- totalgap
  primaryOutcomesLowWhisker$`Total hypoxia (hour)`[j] <- hypox_dura_min #converting hour to minute
  vals <- (90 - data1$SPO2)
  vals[vals < 0] <- 0
  primaryOutcomesLowWhisker$`AUC below 90`[j] <- caTools::trapz(as.numeric(difftime(data1$dev_reading_tm, min(data1$dev_reading_tm), units = "mins")), vals)
}
boxplot(primaryOutcomes$`Hypoxia (min per hour)`, primaryOutcomesLowWhisker$`Hypoxia (min per hour)`, main = "Hypoxia", names = c("below 90", "Below lower whisker"), horizontal = TRUE, xlab = "min/hour")
#covariate analysis ----
ggpairs(intraop[-c(354,98,449),c(9,14,15)])
ggpairs(as.data.frame(cbind(intraop[-c(354,98,449),c(9,14,15)], primaryOutcomesLowWhisker$`Hypoxia (min per hour)`)))

#Gaussian model ----
data3 <- pVisi_data[[7]]$SPO2
fit = Mclust(data3,G = NULL, modelNames="V")
plot(fit, what="density", main="", xlab="SPO2") 
rug(data3)
plot(fit, what = c("classification"))
grid()
plot(fit, what = "BIC")

gaussianG <- as.data.frame(matrix(data = NA, nrow = length(pVisi_data), ncol = 8))
colnames(gaussianG) <- c("ID", "G", "Mean 1", "Mean 2", "Mean 3", "Variance 1", "Variance 2", "Variance 3")
gaussianG$ID <- names(pVisi_data)
gaussianClassification <-  vector("list", length = length(pVisi_data))
gaussianData <-  vector("list", length = length(pVisi_data))
names(gaussianClassification) <- names(pVisi_data)
names(gaussianData) <- names(pVisi_data)
for (j in 2:length(pVisi_data)) {
  print(j)
  data1 <- pVisi_data[[j]]$SPO2
  fit <- Mclust(data1,G = 1:3, modelNames="E")
  gaussianG$G[j] <- fit$G
  gaussianG[j, 3:(fit$G + 2)] <- fit$parameters$mean
  gaussianG[j, 6:(fit$G + 5)] <- fit$parameters$variance$sigmasq
  gaussianClassification[[j]] <- fit$classification
  gaussianData[[j]] <- fit$data
}
hist(gaussianG$V1, main = "Gaussian G value", xlab = "G", xlim = c(1,4))
#patient 7 (G = 1), patient 264 (G = 2), patient 120 (G = 3), patient 460 (G = $)
#PLOTTING
j = 120
data1 <- pVisi_data[[j]]$SPO2

hist(data1, main = "Patient 7 (G = 1)", xlab = "SPO2")

plot(gaussian[[j]], what="density",  xlab="SPO2") 
rug(data1)

id <- "376010035"
data1 <- pVisi_data[[id]]
p1 <- boxplot(data1$SPO2, data1[data1$SPO2 < 90, ]$SPO2, data1[data1$SPO2 >= 90, ]$SPO2, main = "G = 1", names = c("Overall", "Below 90", "Above 90"), horizontal = TRUE, xlab = "SPO2")
id <- "376010009"
data1 <- pVisi_data[[id]]
p2 <- boxplot(data1$SPO2, data1[data1$SPO2 < 90, ]$SPO2, data1[data1$SPO2 >= 90, ]$SPO2, main = "G = 2", names = c("Overall", "Below 90", "Above 90"), horizontal = TRUE, xlab = "SPO2")
id <- "376010001"
data1 <- pVisi_data[[id]]
p3 <- boxplot(data1$SPO2, data1[data1$SPO2 < 90, ]$SPO2, data1[data1$SPO2 >= 90, ]$SPO2, main = "G = 3", names = c("Overall", "Below 90", "Above 90"), horizontal = TRUE, xlab = "SPO2")
p1+p2+p3

j = "376010104"
dataClust <- data.frame(Index = 1:length(gaussianClassification[[j]]),Data = gaussianData[[j]],Classification = gaussianClassification[[j]])
infoText<- paste(names(gaussianG[gaussianG$ID == j,3:8]), round(gaussianG[gaussianG$ID == j,3:8], 2), sep = ": ", collapse = "\n")
ggplot(data = dataClust, aes(x = Index, y = Data)) + geom_line(aes(group = 1)) + geom_point(aes(color = factor(Classification))) +
  labs(title = "G = 1", x = "Time", y = "SPO2", color = "Classification") + theme_minimal()  + scale_color_viridis(discrete = TRUE) +
  geom_hline(yintercept = 90, linetype = "dashed", color = "blue", size = 1) +
  annotate("text", x = 1500, y = 85, label =  infoText, color = "black", size = 2)
j <- "376010009"
dataClust <- data.frame(Index = 1:length(gaussianClassification[[j]]),Data = gaussianData[[j]],Classification = gaussianClassification[[j]])
infoText<- paste(names(gaussianG[gaussianG$ID == j,3:8]), round(gaussianG[gaussianG$ID == j,3:8], 2), sep = ": ", collapse = "\n")
ggplot(data = dataClust, aes(x = Index, y = Data)) + geom_line(aes(group = 1)) + geom_point(aes(color = factor(Classification))) +
  labs(title = "G = 2", x = "Time", y = "SPO2", color = "Classification") + theme_minimal()  + scale_color_viridis(discrete = TRUE) +
  geom_hline(yintercept = 90, linetype = "dashed", color = "blue", size = 1) +
  annotate("text", x = 450, y = 83, label =  infoText, color = "black", size = 2)
j <- "376010002"
dataClust <- data.frame(Index = 1:length(gaussianClassification[[j]]),Data = gaussianData[[j]],Classification = gaussianClassification[[j]])
infoText<- paste(names(gaussianG[gaussianG$ID == j,3:8]), round(gaussianG[gaussianG$ID == j,3:8], 2), sep = ": ", collapse = "\n")
ggplot(data = dataClust, aes(x = Index, y = Data)) + geom_line(aes(group = 1)) + geom_point(aes(color = factor(Classification))) +
  labs(title = "G = 3", x = "Time", y = "SPO2", color = "Classification") + theme_minimal()  + scale_color_viridis(discrete = TRUE) +
  geom_hline(yintercept = 90, linetype = "dashed", color = "blue", size = 1) +
  annotate("text", x = 400, y = 86, label =  infoText, color = "black", size = 2)

gaussianG$ID  <- as.numeric(gaussianG$ID)
gaussian_clinical <- left_join(gaussianG, baseline, by = c("ID" = "randomization_number"))
ggplot(data = gaussian_clinical, aes(x = age, y = `Mean 1`, color = as.factor(gender))) + geom_point()
ggplot(data = gaussian_clinical, aes(x = bmi, y = `Mean 1`, color = as.factor(gender))) + geom_point() + xlim(14, 40)
ggplot(data = gaussian_clinical %>% dplyr::filter(!is.na(`Mean 3`)), aes(x = as.factor(smoking_status), y = `Mean 1`)) + geom_boxplot()
kruskal.test(`Mean 1`~gender, data = gaussian_clinical)
#G = 3 analysis
g3Patients <- rownames(gaussianG)[gaussianG$G == 3] #filtering all the G3 patients
percentage <- as.data.frame(matrix(data = NA, nrow = length(g3Patients), ncol = 4))
colnames(percentage) <- c("ID", "1", "2", "3")
percentage$ID <- g3Patients
for (j in 206:length(g3Patients)) {
  # Calculate percentage of values below 90% in each group
  x <- as.numeric(g3Patients[j])
  percentage[j, 2:4] <- tapply(gaussianData[[x]], gaussianClassification[[x]], function(x) {mean(x < 90) * 100})
}



#Gaussian split across 12 hours ----
gaussianG_2 <- as.data.frame(matrix(data = NA, nrow = length(pVisi_data), ncol = 8))
colnames(gaussianG_2) <- c("ID", "G", "Mean 1", "Mean 2", "Mean 3", "Variance 1", "Variance 2", "Variance 3")
gaussianG_2$ID <- names(pVisi_data)
gaussianClassification_2 <-  vector("list", length = length(pVisi_data))
gaussianData_2 <-  vector("list", length = length(pVisi_data))
names(gaussianClassification_2) <- names(pVisi_data)
names(gaussianData_2) <- names(pVisi_data)
gaussianG_3 <- as.data.frame(matrix(data = NA, nrow = length(pVisi_data), ncol = 8))
colnames(gaussianG_3) <- c("ID", "G", "Mean 1", "Mean 2", "Mean 3", "Variance 1", "Variance 2", "Variance 3")
gaussianG_3$ID <- names(pVisi_data)
gaussianClassification_3 <-  vector("list", length = length(pVisi_data))
gaussianData_3 <-  vector("list", length = length(pVisi_data))
names(gaussianClassification_3) <- names(pVisi_data)
names(gaussianData_3) <- names(pVisi_data)
for (j in 497:length(pVisi_data)) {
  #pVisiData split so that df 1 has the 1st 12 hours
  print(j)
  start_time <- pVisi_data[[j]]$dev_reading_tm[1]
  end_time <- start_time + 12 * 60 * 60
  first12 <- pVisi_data[[j]]$SPO2[pVisi_data[[j]]$dev_reading_tm <= end_time]
  #restTime <- pVisi_data[[j]]$SPO2[pVisi_data[[j]]$dev_reading_tm > end_time]
  restTime <- pVisi_data[[j]]$SPO2
  #setting to equal variance
    fit <- Mclust(first12,G = 1:3, modelNames="E")
    gaussianG_2$G[j] <- fit$G
    gaussianG_2[j, 3:(fit$G + 2)] <- fit$parameters$mean
    gaussianG_2[j, 6:(fit$G + 5)] <- fit$parameters$variance$sigmasq
    gaussianClassification_2[[j]] <- fit$classification
    gaussianData_2[[j]] <- fit$data
  #restTime
    fit <- Mclust(pVisi_data[[j]]$SPO2,G = 1:3, modelNames="E")
    gaussianG_3$G[j] <- fit$G
    gaussianG_3[j, 3:(fit$G + 2)] <- fit$parameters$mean
    gaussianG_3[j, 6:(fit$G + 5)] <- fit$parameters$variance$sigmasq
    gaussianClassification_3[[j]] <- fit$classification
    gaussianData_3[[j]] <- fit$data
    
}

gaussianG_2$which <- "first12"
gaussianG_3$which <- "allTime"
#swapping rows for G = 2
gaussian2_3_2 <- gaussian2_3[gaussian2_3$G == 2,]
swappingRows <- which(gaussian2_3_2$`Mean 1` > gaussian2_3_2$`Mean 2`) 
gaussian2_3_2[swappingRows,] <- gaussian2_3_2[swappingRows,c(1,2,4,3,5,7,6,8,9)]
#visualization
gaussian2_3 <- merge(gaussianG_2, gaussianG_3, by = "ID")

ggplot(gaussian2_3) + geom_boxplot(aes(x = which, y = `Mean 1`, fill = which), na.rm = TRUE)

ggplot(gaussian2_3_2) + geom_boxplot(aes(x = which, y = `Mean 1`, fill = which), na.rm = TRUE) 
ggplot(gaussian2_3_2) + geom_boxplot(aes(x = which, y = `Mean 2`, fill = which), na.rm = TRUE) 

ggplot(gaussian2_3, aes(x = G.x, y = G.y)) + geom_count() + theme_bw() + 
  labs(title = "Equal variance (G in 1st 12 hours vs all time)", x = "1st twelve hours", y = "All time points") +
  geom_text(aes(label = ..n..), stat = "sum", hjust = 1.5, size = 4) 
onlyG1 <- gaussian2_3[((gaussian2_3$G.x == '1') & (gaussian2_3$G.y == '1')),]
boxplot(onlyG1$`Mean 1.x`, onlyG1$`Mean 1.y`, main  = "Comparison of G = 1 SPO2 values", xlab = "1st 12 hours vs rest", ylab = "SPO2 Mean")
stripchart(onlyG1$`Mean 1.x`, method = "jitter", pch = 19, col = 4, add = TRUE) 
#2nd model
j = "376010001"
dataClust <- data.frame(Index = 1:length(gaussianClassification_2[[j]]),Data = gaussianData_2[[j]],Classification = gaussianClassification_2[[j]])
infoText<- paste(names(gaussianG_2[gaussianG_2$ID == j,3:8]), round(gaussianG_2[gaussianG_2$ID == j,3:8], 2), sep = ": ", collapse = "\n")
ggplot(data = dataClust, aes(x = Index, y = Data)) + geom_line(aes(group = 1)) + geom_point(aes(color = factor(Classification))) +
  labs(title = j, x = "Time", y = "SPO2", color = "Classification") + theme_minimal()  + scale_color_viridis(discrete = TRUE) +
  geom_hline(yintercept = 90, linetype = "dashed", color = "blue", size = 1) +
  annotate("text", x = 100, y = 88, label =  infoText, color = "black", size = 2)
#3rd model
j = "376010001"
dataClust <- data.frame(Index = 1:length(gaussianClassification_3[[j]]),Data = gaussianData_3[[j]],Classification = gaussianClassification_3[[j]])
infoText<- paste(names(gaussianG_3[gaussianG_3$ID == j,3:8]), round(gaussianG_3[gaussianG_3$ID == j,3:8], 2), sep = ": ", collapse = "\n")
ggplot(data = dataClust, aes(x = Index, y = Data)) + geom_line(aes(group = 1)) + geom_point(aes(color = factor(Classification))) +
  labs(title = j, x = "Time", y = "SPO2", color = "Classification") + theme_minimal()  + scale_color_viridis(discrete = TRUE) +
  geom_hline(yintercept = 90, linetype = "dashed", color = "blue", size = 1) +
  annotate("text", x = 100, y = 88, label =  infoText, color = "black", size = 2)
