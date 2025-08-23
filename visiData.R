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
#baseline filtered ----
baseline_filtered <- baseline
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
    # else if (dewqata1$SPO2[i-1] < 90) { #if current spo2 not below 90, but previous below then a time period added to it
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
#SPO2 distribution
plot(pVisi_data[[93]]$dev_reading_tm, pVisi_data[[93]]$SPO2, type = "l",
     main = "Sample Patient SPO2", xlab = "Time", ylab = "SPO2 Level")
# plotting primary outcomes ----
primaryOutcomes$totalHypoxiaActualHour <- primaryOutcomes$`Total hypoxia (hour)`/60
p1 <- ggplot(primaryOutcomes, aes(x = factor(1), y = `Total hypoxia (hour)`)) +
  geom_boxplot() +
  labs(y = "Total Hypoxemia (hours)") +
  theme_light() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())
p2 <- ggplot(primaryOutcomes, aes(x = factor(1), y = `Hypoxia (min per hour)`)) +
  geom_boxplot() +
  labs(y = "Hypoxemia burden (min per hour)") +
  theme_light() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())
p3 <- ggplot(primaryOutcomes, aes(x = factor(1), y = `Total reading time`)) +
  geom_boxplot() +
  labs(y = "Total Reading Time (hours)") +
  theme_light() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())
p3 + p1 + p2 + plot_annotation(title = "Boxplots of Primary Outcomes")
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
# data3 <- pVisi_data[[7]]$SPO2
# fit = Mclust(data3,G = NULL, modelNames="V")
# plot(fit, what="density", main="", xlab="SPO2") 
# rug(data3)
# plot(fit, what = c("classification"))
# grid()
# plot(fit, what = "BIC")

gaussianG <- as.data.frame(matrix(data = NA, nrow = length(pVisi_data), ncol = 8))
colnames(gaussianG) <- c("ID", "G", "Mean 1", "Mean 2", "Mean 3", "Variance 1", "Variance 2", "Variance 3")
gaussianG$ID <- names(pVisi_data)
gaussianClassification <-  vector("list", length = length(pVisi_data))
gaussianData <-  vector("list", length = length(pVisi_data))
names(gaussianClassification) <- names(pVisi_data)
names(gaussianData) <- names(pVisi_data)
for (j in 1:length(pVisi_data)) {
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
gaussianG3Patients <- gaussianG[gaussianG$G == 3,]
gaussianG2Patients <- gaussianG[gaussianG$G == 2,]
#patient 7 (G = 1), patient 264 (G = 2), patient 120 (G = 3), patient 460 (G = $)
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
#Gaussian Plotting ----
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
#showing individual traces and histograms
j = "376010104"
dataClust <- data.frame(Index = 1:length(gaussianClassification[[j]]),Data = gaussianData[[j]],Classification = gaussianClassification[[j]])
infoText<- paste(names(gaussianG[gaussianG$ID == j,3:8]), round(gaussianG[gaussianG$ID == j,3:8], 2), sep = ": ", collapse = "\n")
ggplot(data = dataClust, aes(x = Index, y = Data)) + geom_line(aes(group = 1)) + geom_point(aes(color = factor(Classification))) +
  labs(title = "G = 1", x = "Time", y = "SPO2", color = "Classification") + theme_minimal()  + scale_color_viridis(discrete = TRUE) +
  geom_hline(yintercept = 90, linetype = "dashed", color = "blue", size = 1) +
  annotate("text", x = 1500, y = 85, label =  infoText, color = "black", size = 2)
ggplot(data = dataClust, aes(x = Data, fill = factor(Classification))) + 
  geom_histogram(binwidth = 1, alpha = 0.7, fill = "grey60", color = "grey40") + 
  labs(title = "G = 1", x = "SPO2", y = "Count", fill = "Classification") + 
  theme_minimal() + 
  scale_fill_viridis(discrete = TRUE)
j <- "376010026"
dataClust <- data.frame(Index = 1:length(gaussianClassification[[j]]),Data = gaussianData[[j]],Classification = gaussianClassification[[j]])
infoText<- paste(names(gaussianG[gaussianG$ID == j,3:8]), round(gaussianG[gaussianG$ID == j,3:8], 2), sep = ": ", collapse = "\n")
ggplot(data = dataClust, aes(x = Index, y = Data)) + geom_line(aes(group = 1)) + geom_point(aes(color = factor(Classification)), size = 1) +
  labs(title = "G = 2", x = "Time", y = "SPO2", color = "Classification") + theme_minimal()  + scale_color_viridis(discrete = TRUE) 
  geom_hline(yintercept = 90, linetype = "dashed", color = "blue", size = 1) +
  annotate("text", x = 450, y = 83, label =  infoText, color = "black", size = 2)
ggplot(data = dataClust, aes(x = Data, fill = factor(Classification))) + 
  geom_histogram(binwidth = 1, alpha = 0.7, fill = "grey60", color = "grey40") + 
  labs(title = "G = 2", x = "SPO2", y = "Count", fill = "Classification") + 
  theme_minimal() + 
  scale_fill_viridis(discrete = TRUE)
j <- "376010019"
dataClust <- data.frame(Index = 1:length(gaussianClassification[[j]]),Data = gaussianData[[j]],Classification = gaussianClassification[[j]])
infoText<- paste(names(gaussianG[gaussianG$ID == j,3:8]), round(gaussianG[gaussianG$ID == j,3:8], 2), sep = ": ", collapse = "\n")
ggplot(data = dataClust, aes(x = Index, y = Data)) + geom_line(aes(group = 1)) + geom_point(aes(color = factor(Classification)), size = 1) +
  labs(title = "G = 3", x = "Time", y = "SPO2", color = "Classification") + theme_minimal()  + scale_color_viridis(discrete = TRUE) 
  geom_hline(yintercept = 90, linetype = "dashed", color = "blue", size = 1) +
  annotate("text", x = 400, y = 86, label =  infoText, color = "black", size = 2)
ggplot(data = dataClust, aes(x = Data, fill = factor(Classification))) + 
  geom_histogram(binwidth = 1, alpha = 0.7, fill = "grey60", color = "grey40") + 
  labs(title = "G = 3", x = "SPO2", y = "Count", fill = "Classification") + 
  theme_minimal() + 
  scale_fill_viridis(discrete = TRUE)
# weird patient
j <- "376010065"
dataClust <- data.frame(Index = 1:length(gaussianClassification[[j]]),Data = gaussianData[[j]],Classification = gaussianClassification[[j]])
ggplot(data = dataClust, aes(x = Data, fill = factor(Classification))) + 
  geom_histogram(binwidth = 1, alpha = 0.7, fill = "grey60", color = "grey40") + 
  labs(title = "False positive patient", x = "SPO2", y = "Count", fill = "Classification") + 
  theme_minimal() + 
  scale_fill_viridis(discrete = TRUE)
ggplot(data = dataClust, aes(x = Index, y = Data)) + geom_line(aes(group = 1)) + geom_point(aes(color = factor(Classification)), size = 1) +
  labs(title = "False positive patient", x = "Time", y = "SPO2", color = "Classification") + theme_minimal()  + scale_color_viridis(discrete = TRUE) 
  geom_hline(yintercept = 90, linetype = "dashed", color = "blue", size = 1) 
  
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

#intra-op 
intraop <- intraop[which(intraop$randomization_number %in% names(pVisi_data)),]
summary(factor(intraop$Type.of.Surgery))

# histogram of ggplot distributions
# Create histogram
gaussian_long <- gaussianG3Patients %>%
  pivot_longer(cols = starts_with("Mean"), names_to = "Mean_Type", values_to = "SpO2")
# Create histogram with 75% transparency
ggplot(gaussian_long, aes(x = SpO2, fill = Mean_Type)) +
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.65) +  # Set transparency to 75%
  geom_vline(xintercept = 90, linetype = "dashed", color = "black", linewidth = 1) +  # Vertical line at 90
  scale_fill_manual(values = c("#8E7CC3", "#76AFAE", "#F7E67C")) +  # Colors for means
  labs(x = "Mean SpO2", y = "Frequency", fill = "Mean Type") +  
  theme_minimal() + 
  theme(legend.position = "top")
gaussian_long <- gaussianG2Patients %>%
  pivot_longer(cols = starts_with("Mean"), names_to = "Mean_Type", values_to = "SpO2")
# Create histogram with 75% transparency
ggplot(gaussian_long, aes(x = SpO2, fill = Mean_Type)) +
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.65) +  # Set transparency to 75%
  geom_vline(xintercept = 90, linetype = "dashed", color = "black", linewidth = 1) +  # Vertical line at 90
  scale_fill_manual(values = c("#8E7CC3", "#F7E67C")) +  # Colors for means
  labs(x = "Mean SpO2", y = "Frequency", fill = "Mean Type") +  
  theme_minimal() + 
  theme(legend.position = "top")
#boxplots
gaussian_long <- gaussianG3Patients %>%
  pivot_longer(cols = starts_with("Mean"), names_to = "Mean_Type", values_to = "SpO2")
ggplot(gaussian_long, aes(x = Mean_Type, y = SpO2, fill = Mean_Type)) +
  geom_boxplot(alpha = 0.65) +
  geom_hline(yintercept = 90, linetype = "dashed", color = "black", linewidth = 1) +
  scale_fill_manual(values = c("#8E7CC3", "#76AFAE", "#F7E67C")) +
  labs(x = "Mean Type", y = "SpO2", fill = "Mean Type") +
  theme_minimal() + ylim(70, 100) +
  theme(legend.position = "none") + coord_flip()
gaussian_long <- gaussianG2Patients %>%
  pivot_longer(cols = starts_with("Mean"), names_to = "Mean_Type", values_to = "SpO2")
ggplot(gaussian_long, aes(x = Mean_Type, y = SpO2, fill = Mean_Type)) +
  geom_boxplot(alpha = 0.65) +
  geom_hline(yintercept = 90, linetype = "dashed", color = "black", linewidth = 1) +
  scale_fill_manual(values = c("#8E7CC3", "#F7E67C")) +
  labs(x = "Mean Type", y = "SpO2", fill = "Mean Type") +
  theme_minimal() + ylim(70, 100) +
  theme(legend.position = "none") + coord_flip()


#Gaussian thresholds equal probabilities ----
gaussianG2Patients$midpoint <- rowMeans(gaussianG2Patients[, c("Mean 1", "Mean 2")])
sum(gaussianG2Patients$midpoint >= 87 & gaussianG2Patients$midpoint <= 93)
sum(gaussianG2Patients$midpoint >= 87 & gaussianG2Patients$midpoint <= 93)/nrow(gaussianG2Patients)

sum(gaussianG2Patients$`Mean 1` >= 87 & gaussianG2Patients$`Mean 1` <= 93)
sum(gaussianG2Patients$`Mean 1` >= 87 & gaussianG2Patients$`Mean 1` <= 93)/nrow(gaussianG2Patients)

gaussianG3Patients$mid1_2 <- rowMeans(gaussianG3Patients[, c("Mean 1", "Mean 2")])
gaussianG3Patients$mid2_3 <- rowMeans(gaussianG3Patients[, c("Mean 2", "Mean 3")])

sum(gaussianG3Patients$mid1_2 >= 87 & gaussianG3Patients$mid1_2 <= 93)
sum(gaussianG3Patients$mid1_2 >= 87 & gaussianG3Patients$mid1_2 <= 93)/nrow(gaussianG3Patients)
# Create summary data
summary_data <- data.frame(
  Statistic = c("Min", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max"),
  `G2 Patients Midpoint` = c(77.23, 90.39, 92.79, 92.09, 94.63, 98.83),
  `G3 Patients Mid1_2` = c(77.23, 90.27, 92.33, 91.73, 94.18, 98.50),
  `G3 Patients Mid2_3` = c(88.23, 94.27, 95.84, 95.48, 97.05, 99.50)
)

# Display the table
kable(summary_data, caption = "Summary Statistics for Gaussian Groups")
#plotting thresholds
#G = 2
ggplot(gaussianG2Patients, aes(x = midpoint)) +
  geom_histogram(binwidth = 1, fill = "pink3", position = "identity", alpha = 0.65) +  # Set transparency to 75%
  geom_vline(xintercept = 90, linetype = "dashed", color = "black", linewidth = 1) + 
  labs(x = "Thresholds", y = "Frequency") +  
  theme_minimal() + 
  theme(legend.position = "top") + coord_flip()
ggplot(gaussianG2Patients, aes(y = midpoint)) +
  geom_boxplot(binwidth = 1, fill = "pink3", position = "identity", alpha = 0.65) +  # Set transparency to 75%
  geom_hline(yintercept = 90, linetype = "dashed", color = "black", linewidth = 1) + 
  labs(y = "Thresholds") +  
  theme_minimal() + 
  theme(legend.position = "top")
#G = 3
gaussian_long <- gaussianG3Patients %>%
  pivot_longer(cols = starts_with("mid"), names_to = "type", values_to = "midpoint")
# Create histogram with 75% transparency
ggplot(gaussian_long, aes(x = midpoint, fill = type)) +
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.65) +  # Set transparency to 75%
  geom_vline(xintercept = 90, linetype = "dashed", color = "black", linewidth = 1) +  # Vertical line at 90
  scale_fill_manual(values = c("purple2", "orange3")) +  # Colors for means
  labs(x = "Thresholds", y = "Frequency", fill = "Mean Type") +  
  theme_minimal() + 
  theme(legend.position = "top") + coord_flip()
sggplot(gaussian_long, aes(x = type, y = midpoint, fill = type)) +
  geom_boxplot(alpha = 0.65) +  # Set transparency to 75%
  geom_hline(yintercept = 90, linetype = "dashed", color = "black", linewidth = 1) +  # Horizontal line at 90
  scale_fill_manual(values = c("purple2", "orange3")) +  # Colors for means
  labs(x = "Mean Type", y = "Thresholds", fill = "Mean Type") +  
  theme_minimal() + 
  theme(legend.position = "none")  # Legend not needed as x-axis already shows type

#plotting thresholds with clinical probabilities
gaussianG2Patients <- cbind(
  gaussianG2Patients,
  baseline[match(gaussianG2Patients$ID, baseline$randomization_number), ]
)
gaussianG3Patients <- cbind(
  gaussianG3Patients,
  baseline[match(gaussianG3Patients$ID, baseline$randomization_number), ]
)
gaussianG2Patients$stopbang_risk <- as.factor(gaussianG2Patients$stopbang_risk)
gaussianG2Patients$asthma <- as.factor(gaussianG2Patients$asthma)
gaussianG2Patients$smoking_status <- as.factor(gaussianG2Patients$smoking_status)
gaussianG2Patients$chronic_pain <- as.factor(gaussianG2Patients$chronic_pain)
gaussianG2Patients$gender <- as.factor(gaussianG2Patients$gender)
gaussianG2Patients$asa <- as.factor(gaussianG2Patients$asa)

gaussianG3Patients$stopbang_risk <- as.factor(gaussianG3Patients$stopbang_risk)
gaussianG3Patients$asthma <- as.factor(gaussianG3Patients$asthma)
gaussianG3Patients$smoking_status <- as.factor(gaussianG3Patients$smoking_status)
gaussianG3Patients$chronic_pain <- as.factor(gaussianG3Patients$chronic_pain)

ggplot(gaussianG2Patients, aes(x = asthma, y = midpoint, fill = asthma)) +
  stat_compare_means(method = "t.test", label = "p.format") +  
  geom_boxplot() +
  labs(title = "Threshold distribution by Asthma", x = "Asthma: Yes or No",
       y = "Midpoint threshold") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(gaussianG2Patients, aes(x = asa, y = midpoint, fill = asa)) +
  stat_compare_means(method = "anova", label = "p.format") +  
  geom_boxplot() +
  labs(title = "Threshold distribution by asa G = 2", x = "ASA",
       y = "Midpoint threshold") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(gaussianG2Patients, aes(x = smoking_status, y = midpoint, fill = smoking_status)) +
  stat_compare_means(method = "t.test", label = "p.format") +  
  geom_boxplot() +
  labs(title = "Threshold distribution by smoking status", x = "Smoking: Yes or No",
       y = "Midpoint threshold") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(gaussianG2Patients, aes(x = chronic_pain, y = midpoint, fill = chronic_pain)) +
  stat_compare_means(method = "t.test", label = "p.format") +  
  geom_boxplot() +
  labs(title = "Threshold distribution by chronic_pain", x = "chronic_pain: Yes or No",
       y = "Midpoint threshold") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(gaussianG2Patients, aes(x = stopbang_risk, y = midpoint, fill = stopbang_risk)) +
  stat_compare_means(method = "t.test", label = "p.format") +  
  geom_boxplot() +
  labs(title = "Threshold distribution by STOPBANG Risk", x = "Risk: Yes or No",
    y = "Midpoint threshold") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(gaussianG2Patients, aes(x = gender, y = midpoint, fill = gender)) +
  stat_compare_means(method = "t.test", label = "p.format") +  
  geom_boxplot() +
  labs(title = "Threshold distribution by Gender", x = "Gender",
       y = "Midpoint threshold") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(gaussianG2Patients, aes(x = age, y = midpoint)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  labs(title = "Threshold distribution vs Age",
    x = "Age", y = "Midpoint") + theme_minimal()
ggplot(gaussianG2Patients, aes(x = bmi, y = midpoint)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  labs(title = "Threshold distribution vs BMI",
       x = "BMI", y = "Midpoint") + theme_minimal()

ggplot(gaussianG3Patients |> dplyr::filter(!is.na(smoking_status)), aes(x = smoking_status, y = mid2_3, fill = smoking_status)) +
  stat_compare_means(method = "t.test", label = "p.format") +  
  geom_boxplot() +
  labs(title = "Threshold distribution by smoking status G = 3", x = "Smoking: Yes or No",
       y = "Midpoint threshold") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# gaussian thresholds P<0.05 ----
gaussianG3Patients$threshold_0.05 <- qnorm(0.05, mean = gaussianG3Patients$`Mean 3`, sd = sqrt(gaussianG3Patients$`Variance 3`))
gaussianG2Patients$threshold_0.05 <- qnorm(0.05, mean = gaussianG2Patients$`Mean 2`, sd = sqrt(gaussianG2Patients$`Variance 2`))
ggplot(gaussianG2Patients, aes(x = threshold_0.05)) +
  geom_histogram(binwidth = 1, fill = "pink3", position = "identity", alpha = 0.65) +  # Set transparency to 75%
  geom_vline(xintercept = 90, linetype = "dashed", color = "black", linewidth = 1) + 
  labs(x = "Thresholds", y = "Frequency", title = "G = 2 P<0.05 threshold") +  
  theme_minimal() + 
  theme(legend.position = "top")
ggplot(gaussianG3Patients, aes(x = threshold_0.05)) +
  geom_histogram(binwidth = 1, fill = "pink3", position = "identity", alpha = 0.65) +  # Set transparency to 75%
  geom_vline(xintercept = 90, linetype = "dashed", color = "black", linewidth = 1) + 
  labs(x = "Thresholds", y = "Frequency", title = "G = 3 P<0.05 threshold") +  
  theme_minimal() + 
  theme(legend.position = "top")

# gaussian difference between thresholds ----
gaussianG2Patients$diffThresholds <- gaussianG2Patients$midpoint - gaussianG2Patients$threshold_0.05
ggplot(gaussianG2Patients, aes(x = midpoint)) +
  geom_histogram(binwidth = 1, fill = "grey60", position = "identity", alpha = 0.65) +  # Set transparency to 75%
  labs(x = "Thresholds", y = "equal - p<0.05", title = "G = 2 difference in thresholds") +  
  theme_minimal() + 
  theme(legend.position = "top")
bw <- 1 # specify width of bins in histogram
p1 <- ggplot(gaussianG2Patients, aes(midpoint)) +
  geom_histogram(binwidth = bw, fill = "black", col = "yellow" ) +
  stat_function(fun = function(x) {
    dnorm(x, mean = mean(gaussianG2Patients$midpoint, na.rm = TRUE),
          sd = sd(gaussianG2Patients$midpoint, na.rm = TRUE) ) *
      length(gaussianG2Patients$midpoint) * bw},
    geom = "area", alpha = 0.5, fill = "lightblue", col = "blue") +
  labs(x = "threshold",
       title = "")
p2 <- ggplot(gaussianG2Patients, aes(sample = midpoint)) +
  geom_qq() +
  geom_qq_line(col = "red") +
  labs(y = "",
       title = "")
p3 <- ggplot(gaussianG2Patients, aes(x = midpoint, y = "")) +
  geom_violin(fill = "cornsilk") +
  geom_boxplot(width = 0.2) +
  stat_summary(fun = mean, geom = "point", shape = 16, col = "red") +
  labs(y = "", x = "", title = "Boxplot with Violin")
p1 + (p2 / p3 + plot_layout(heights = c(2, 1))) +
  plot_annotation(title = "equal probability threshold")
bw <- 1 # specify width of bins in histogram
p1 <- ggplot(gaussianG2Patients, aes(threshold_0.05)) +
  geom_histogram(binwidth = bw, fill = "black", col = "yellow" ) +
  stat_function(fun = function(x) {
    dnorm(x, mean = mean(gaussianG2Patients$threshold_0.05, na.rm = TRUE),
          sd = sd(gaussianG2Patients$threshold_0.05, na.rm = TRUE) ) *
      length(gaussianG2Patients$threshold_0.05) * bw},
    geom = "area", alpha = 0.5, fill = "lightblue", col = "blue") +
  labs(x = "threshold",
       title = "")
p2 <- ggplot(gaussianG2Patients, aes(sample = threshold_0.05)) +
  geom_qq() +
  geom_qq_line(col = "red") +
  labs(y = "",
       title = "")
p3 <- ggplot(gaussianG2Patients, aes(x = threshold_0.05, y = "")) +
  geom_violin(fill = "cornsilk") +
  geom_boxplot(width = 0.2) +
  stat_summary(fun = mean, geom = "point", shape = 16, col = "red") +
  labs(y = "", x = "", title = "Boxplot with Violin")
p1 + (p2 / p3 + plot_layout(heights = c(2, 1))) +
  plot_annotation(title = "p<0.05 threshold")
# sample gussian distribution ----
patient_index <- "376010001"
# Pull stored data
data1 <- gaussianData[[patient_index]]
means <- as.numeric(gaussianG[1, 3:(2 + G)])
variance <- as.numeric(gaussianG[1, 6])
props <- table(gaussianClassification[[patient_index]]) / length(gaussianClassification[[patient_index]])
# Create x-axis range
x_vals <- seq(min(data1, na.rm = TRUE) - 2, max(data1, na.rm = TRUE) + 2, length.out = 1000)
# Build Gaussian components
G <- 3
gauss_df <- data.frame(
  x = rep(x_vals, G),
  density = unlist(lapply(1:G, function(k) {
    props[k] * dnorm(x_vals, mean = means[k], sd = sqrt(variance))
  })),
  Cluster = factor(rep(paste0("Cluster ", 1:G), each = length(x_vals)))
)
thresholds <- rowMeans(embed(sort(means), 2))
vline_df <- data.frame(x = thresholds)
cluster_colors <- c("Cluster 1" = "#8E7CC3",   # soft purple
                    "Cluster 2" = "#76AFAE",   # soft teal
                    "Cluster 3" = "#F7E67C") 
# Create histogram with Gaussian overlays
ggplot(data.frame(SpO2 = data1), aes(x = SpO2)) +
  geom_histogram(aes(y = ..density..), binwidth = 1, fill = "gray70", color = "white") +
  geom_line(data = gauss_df, aes(x = x, y = density, color = Cluster), size = 1) +
  labs(
    title = paste("Patient", patient_index, "| G =", G),
    x = "SpO2", y = "Density"
  ) +
  geom_vline(data = vline_df, aes(xintercept = x), linetype = "dashed", size = 0.8, color = "black") +
  scale_color_manual(values = cluster_colors) +
  theme_minimal(base_size = 14) 
