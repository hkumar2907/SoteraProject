#libraries ----
library(randomForest)
library(Epi)
library(caret)
library(dplyr)
library(ggplot2)
library(pROC)
library(PRROC)
# previous RF ----
#testSet <- sample(1:543, round(543 * 0.3))
#trainSet <- setdiff(1:543, testSet)
#summary_windowsTest <- na.omit(summary_windows_cont[summary_windows_cont$patientID %in% testSet, ])
#summary_windowsTrain <- na.omit(summary_windows_cont[summary_windows_cont$patientID %in% trainSet, ])
ggplot(summary_windowsTrain, aes(x = "", y = pre_hyp)) +
  geom_boxplot(fill = "skyblue", color = "black", outlier.color = "red", outlier.shape = 16) +
  labs(title = "Pre-hypoxia window distribution in train set", x = "", y = "Pre-hypoxia window") + 
  theme_minimal()
# feature model 5 hypoxia groups
set.seed(1111)
summary_windows_cont <- summary_windows_cont %>%
  mutate(post_hyp_grouped = case_when(post_hyp >= 0 & post_hyp <= 0.2 ~ "1", 
                           post_hyp > 0.2 & post_hyp <= 0.4 ~ "2",
                           post_hyp > 0.4 & post_hyp <= 0.6 ~ "3", 
                           post_hyp > 0.6 & post_hyp <= 0.8 ~ "4",
                           post_hyp > 0.8 & post_hyp <= 1   ~ "5"))
summary_windows_cont$post_hyp_grouped <- as.factor(summary_windows_cont$post_hyp_grouped)
c <- randomForest(x = summary_windowsTrain[,8:25], y = summary_windowsTrain$post_hyp_grouped, ntree = 1000) 
pred <- predict(c, newdata = summary_windowsTest[,8:25]) 
confusion_mtx <- table(summary_windowsTest$post_hyp_grouped, pred) 
confusion_mtx
# feature model true/false hypoxia
summary_windows_cont <- summary_windows_cont %>% 
  mutate(post_hyp_binary = case_when(post_hyp >= 0.9 ~ "TRUE", TRUE ~ "FALSE"))
summary_windows_cont$post_hyp_binary <- as.factor(summary_windows_cont$post_hyp_binary)
summary_windowsTest <- na.omit(summary_windows_cont[summary_windows_cont$patientID %in% testSet, ])
summary_windowsTrain <- na.omit(summary_windows_cont[summary_windows_cont$patientID %in% trainSet, ])
set.seed(1111)
c2 <- randomForest(x = summary_windowsTrain[,8:25], y = summary_windowsTrain$post_hyp_binary, ntree = 1000) 
pred2 <- predict(c2, newdata = summary_windowsTest[,8:25]) 
confusion_mtx <- table(summary_windowsTest$post_hyp_binary, pred1) 
confusion_mtx
# feature model true/false hypoxia + longest_seq
c3 <- randomForest(x = summary_windowsTrain[,c(8:25, 48)], y = summary_windowsTrain$post_hyp_binary, ntree = 1000) 
pred3 <- predict(c3, newdata = summary_windowsTest[,c(8:25, 48)]) 
confusion_mtx <- table(summary_windowsTest$post_hyp_binary, pred3) 
confusion_mtx
# feature model true/false hypoxia + longest_seq + maximum gap
c4 <- randomForest(x = summary_windowsTrain[,c(8:25, 48, 49)], y = summary_windowsTrain$post_hyp_binary, ntree = 1000) 
pred4 <- predict(c4, newdata = summary_windowsTest[,c(8:25, 48, 49)]) 
confusion_mtx <- table(summary_windowsTest$post_hyp_binary, pred4) 
confusion_mtx
# feature model true/false hypoxia  + maximum gap
c5 <- randomForest(x = summary_windowsTrain[,c(8:25, 49)], y = summary_windowsTrain$post_hyp_binary, ntree = 1000) 
pred5 <- predict(c5, newdata = summary_windowsTest[,c(8:25, 49)]) 
confusion_mtx <- table(summary_windowsTest$post_hyp_binary, pred5) 
confusion_mtx


# caret package r ----
#removing deteroriating windows
#summary_windows_cont <- summary_windows_cont[!(summary_windows_cont$pre_hyp < 0.9 & summary_windows_cont$post_hyp > 0.9), ]
#summary_windows_rfmodel <- rbind(stable_summary_windows, detoriating_summary_windows)
#summary_windows_rfmodel <- rbind(stable_summary_windows_10, detoriating_summary_windows_10)
summary_windows_rfmodel <- rbind(stable_summary_windows_50, detoriating_summary_windows_50)
summary_windows_rfmodel <- rbind(stable_summary_windows50_90minus, detoriating_summary_windows50_90minus)
summary_windows_rfmodel <- rbind(summaries_stable_10, summaries_deter_10)
summary_windows_rfmodel <- rbind(summaries_stable_30, summaries_deter_30)
summary_windows_rfmodel <- rbind(summaries_stable_50, summaries_deter_50)
summary_windows_rfmodel <- rbind(summaries_stable_70, summaries_deter_70)
summary_windows_rfmodel <- rbind(summaries_stable_90, summaries_deter_90)


  #model5_10perc <- model_5
  #model5_30perc <- model_5
  model5_50perc <- model_5
# train test partition
  patientIDs <- unique(summary_windows_rfmodel$patientID)
  set.seed(11)
  trainSet <- createDataPartition(patientIDs, p=0.7, list=FALSE)
  testSet <- setdiff(patientIDs, trainSet)
# building dataset
summary_windows_rfmodel <- summary_windows_rfmodel %>%
  mutate(post_hyp_grouped = case_when(post_hyp >= 0 & post_hyp <= 0.2 ~ "1", 
                                      post_hyp > 0.2 & post_hyp <= 0.4 ~ "2",
                                      post_hyp > 0.4 & post_hyp <= 0.6 ~ "3", 
                                      post_hyp > 0.6 & post_hyp <= 0.8 ~ "4",
                                      post_hyp > 0.8 & post_hyp <= 1   ~ "5"))
summary_windows_rfmodel$post_hyp_grouped <- as.factor(summary_windows_rfmodel$post_hyp_grouped)
summary_windows_rfmodel <- summary_windows_rfmodel %>% 
  mutate(post_hyp_binary = case_when(post_hyp >= 0.9 ~ "TRUE", TRUE ~ "FALSE"))
summary_windows_rfmodel$post_hyp_binary <- factor(summary_windows_rfmodel$post_hyp_binary, levels = c(TRUE, FALSE))
summary_windowsTest <- na.omit(summary_windows_rfmodel[summary_windows_rfmodel$patientID %in% testSet, ])
summary_windowsTrain <- na.omit(summary_windows_rfmodel[summary_windows_rfmodel$patientID %in% trainSet, ])

# RF, caret w/ 5 binary first features only ----
columns <- c(colnames(summary_windowsTrain)[8:25], "post_hyp_binary")
summary_windowsTrain2 <- summary_windowsTrain[,columns]
summary_windowsTest2 <- summary_windowsTest[,columns]
model_1 = train(post_hyp_binary ~ ., data=summary_windowsTrain2, 
                method='rf', trControl = trainControl(method = "cv", number = 5))
pred <- predict(model_1, newdata = summary_windowsTest2)
confusionMatrix(pred, summary_windowsTest2$post_hyp_binary)
confusionMatrix(data = pred, reference = summary_windowsTest2$post_hyp_binary, mode = "prec_recall")
dim(summary_windowsTrain2)
dim(summary_windowsTest2)
importance <- varImp(model_1, scale = FALSE)
plot(importance, main = "RF, caret w/ 5 cv, first features only")
# RF, caret w/ 5 binary first features + adding longest seq + maximum gap ----
columns <- c(colnames(summary_windowsTrain)[8:25], "post_hyp_binary", "longest_hypoxic_seq", "maximumGap")
summary_windowsTrain2 <- summary_windowsTrain[,columns]
summary_windowsTest2 <- summary_windowsTest[,columns]
model_2 = train(post_hyp_binary ~ ., data=summary_windowsTrain2, 
                method='rf', trControl = trainControl(method = "cv", number = 5))
pred <- predict(model_2, newdata = summary_windowsTest2)
confusionMatrix(pred, summary_windowsTest2$post_hyp_binary)
importance <- varImp(model_2, scale = FALSE)
plot(importance, main = "RF, caret w/ 5 cv, first features + longest seq + maximum gap")
# RF, caret w/ 5 binary first features + adding longest seq ----
columns <- c(colnames(summary_windowsTrain)[8:25], "post_hyp_binary", "longest_hypoxic_seq")
summary_windowsTrain2 <- summary_windowsTrain[,columns]
summary_windowsTest2 <- summary_windowsTest[,columns]
model_3 = train(post_hyp_binary ~ ., data=summary_windowsTrain2, 
                method='rf', trControl = trainControl(method = "cv", number = 5))
pred <- predict(model_3, newdata = summary_windowsTest2)
confusionMatrix(pred, summary_windowsTest2$post_hyp_binary)
confusionMatrix(data = pred, reference = summary_windowsTest2$post_hyp_binary, mode = "prec_recall")
dim(summary_windowsTrain2)
dim(summary_windowsTest2)
importance <- varImp(model_3, scale = FALSE)
plot(importance, main = "RF, caret w/ 5 cv, first features + longest seq")
# RF, caret w/ 5 grouped first features + adding longest seq----
columns <- c(colnames(summary_windowsTrain)[8:25], "post_hyp_grouped", "longest_hypoxic_seq")
summary_windowsTrain2 <- summary_windowsTrain[,columns]
summary_windowsTest2 <- summary_windowsTest[,columns]
model_4 = train(post_hyp_grouped ~ ., data=summary_windowsTrain2, 
                method='rf', trControl = trainControl(method = "cv", number = 5))
pred <- predict(model_4, newdata = summary_windowsTest2)
confusionMatrix(pred, summary_windowsTest2$post_hyp_grouped)
importance <- varImp(model_4, scale = FALSE)
plot(importance, main = "RF, caret w/ 5 cv, grouped, first features + longest seq")
# RF, caret w/ 5 binary first features + adding longest seq, balance windows with 2:1 ----
columns <- c(colnames(summary_windowsTrain)[8:25], "post_hyp_binary", "longest_hypoxic_seq")
summary_windowsTrain2 <- summary_windowsTrain %>% group_by(post_hyp_binary) %>%
  slice_sample(n = 2*min(table(summary_windowsTrain$post_hyp_binary))) 
summary_windowsTrain2 <- summary_windowsTrain2[,columns]
summary_windowsTest2 <- summary_windowsTest %>% group_by(post_hyp_binary) %>%
  slice_sample(n = 2*min(table(summary_windowsTest$post_hyp_binary)))
summary_windowsTest2 <- summary_windowsTest2[,columns]
model_4 = train(post_hyp_binary ~ ., data=summary_windowsTrain2, 
                method='rf', trControl = trainControl(method = "cv", number = 5))
pred <- predict(model_4, newdata = summary_windowsTest2)
confusionMatrix(pred, summary_windowsTest2$post_hyp_binary)
confusionMatrix(data = pred, reference = summary_windowsTest2$post_hyp_binary, mode = "prec_recall")
dim(summary_windowsTrain2)
dim(summary_windowsTest2)
importance <- varImp(model_4, scale = FALSE)
plot(importance, main = "RF, caret w/ 5 cv, binary, 2:1, first features + longest seq")
# testing biased model on unbiased model
columns <- c(colnames(summary_windowsTrain)[8:25], "post_hyp_binary", "longest_hypoxic_seq")
summary_windowsTest2 <- summary_windowsTest[,columns]
pred <- predict(model_4, newdata = summary_windowsTest2)
confusionMatrix(pred, summary_windowsTest2$post_hyp_binary)
confusionMatrix(data = pred, reference = summary_windowsTest2$post_hyp_binary, mode = "prec_recall")
dim(summary_windowsTrain2)
dim(summary_windowsTest2)
importance <- varImp(model_4, scale = FALSE)
plot(importance, main = "RF, caret w/ 5 cv, binary, 2:1, first features + longest seq, biased tested on unbiased")
# RF, caret w/ 5 binary first features only, balance windows with 2 non-hypoxic : 1 hypoxic windows ----
columns <- c(colnames(summary_windowsTrain)[8:25], "post_hyp_binary")
summary_windowsTrain2 <- summary_windowsTrain %>% group_by(post_hyp_binary) %>%
  slice_sample(n = 2*min(table(summary_windowsTrain$post_hyp_binary))) 
summary_windowsTrain2 <- summary_windowsTrain2[,columns]
summary_windowsTest2 <- summary_windowsTest %>% group_by(post_hyp_binary) %>%
  slice_sample(n = 2*min(table(summary_windowsTest$post_hyp_binary)))
summary_windowsTest2 <- summary_windowsTest2[,columns]
model_4 = train(post_hyp_binary ~ ., data=summary_windowsTrain2, 
                method='rf', trControl = trainControl(method = "cv", number = 5))
pred <- predict(model_4, newdata = summary_windowsTest2)
confusionMatrix(pred, summary_windowsTest2$post_hyp_binary)
importance <- varImp(model_4, scale = FALSE)
plot(importance, main = "RF, caret w/ 5 cv, binary, 2:1, first features only")
# RF, caret w/ 5 binary first features + adding longest seq, balance windows with 1:1 windows ----
columns <- c(colnames(summary_windowsTrain)[8:25], "post_hyp_binary", "longest_hypoxic_seq")
summary_windowsTrain2 <- summary_windowsTrain %>% group_by(post_hyp_binary) %>%
  slice_sample(n = min(table(summary_windowsTrain$post_hyp_binary))) 
summary_windowsTrain2 <- summary_windowsTrain2[,columns]
summary_windowsTest2 <- summary_windowsTest %>% group_by(post_hyp_binary) %>%
  slice_sample(n = min(table(summary_windowsTest$post_hyp_binary)))
summary_windowsTest2 <- summary_windowsTest2[,columns]
model_5 = train(post_hyp_binary ~ ., data=summary_windowsTrain2, 
                method='rf', trControl = trainControl(method = "cv", number = 5))
pred <- predict(model_5, newdata = summary_windowsTest2)
confusionMatrix(pred, summary_windowsTest2$post_hyp_binary)
confusionMatrix(data = pred, reference = summary_windowsTest2$post_hyp_binary, mode = "prec_recall")
dim(summary_windowsTrain2)
dim(summary_windowsTest2)
importance <- varImp(model_4, scale = FALSE)
plot(importance, main = "RF, caret w/ 5 cv, binary, 1:1, first features + longest seq")
# testing 1:1 model on 2:1 testing data
summary_windowsTest2 <- summary_windowsTest %>% group_by(post_hyp_binary) %>%
  slice_sample(n = 2*min(table(summary_windowsTest$post_hyp_binary)))
summary_windowsTest2 <- summary_windowsTest2[,columns]
pred <- predict(model_4, newdata = summary_windowsTest2)
confusionMatrix(pred, summary_windowsTest2$post_hyp_binary)

# testing biased model on unbiased model
columns <- c(colnames(summary_windowsTrain)[8:25], "post_hyp_binary", "longest_hypoxic_seq")
summary_windowsTest2 <- summary_windowsTest %>% group_by(post_hyp_binary) %>%
  slice_sample(n = 2*min(table(summary_windowsTest$post_hyp_binary)))
summary_windowsTest2 <- summary_windowsTest2[,columns]
pred <- predict(model_4, newdata = summary_windowsTest2)
confusionMatrix(pred, summary_windowsTest2$post_hyp_binary)
importance <- varImp(model_4, scale = FALSE)
# RF, caret w/ 5 binary first features only, balance windows with ratio of 1:1 windows ----
columns <- c(colnames(summary_windowsTrain)[8:25], "post_hyp_binary")
summary_windowsTrain2 <- summary_windowsTrain %>% group_by(post_hyp_binary) %>%
  slice_sample(n = min(table(summary_windowsTrain$post_hyp_binary))) 
summary_windowsTrain2 <- summary_windowsTrain2[,columns]
summary_windowsTest2 <- summary_windowsTest %>% group_by(post_hyp_binary) %>%
  slice_sample(n = min(table(summary_windowsTest$post_hyp_binary)))
summary_windowsTest2 <- summary_windowsTest2[,columns]
model_4 = train(post_hyp_binary ~ ., data=summary_windowsTrain2, 
                method='rf', trControl = trainControl(method = "cv", number = 5))
pred <- predict(model_4, newdata = summary_windowsTest2)
confusionMatrix(pred, summary_windowsTest2$post_hyp_binary)
importance <- varImp(model_4, scale = FALSE)
plot(importance, main = "RF, caret w/ 5 cv, binary, 1:1, first features only")
# TUNING  ----
columns <- c(colnames(summary_windowsTrain)[8:25], "post_hyp_binary", "longest_hypoxic_seq")
summary_windowsTrain2 <- summary_windowsTrain %>% group_by(post_hyp_binary) %>%
  slice_sample(n = 2*min(table(summary_windowsTrain$post_hyp_binary))) 
summary_windowsTrain2 <- summary_windowsTrain2[,columns]
# tuning n_trees
store_maxtrees <- list()
for (nt in c(250, 300, 350, 400, 450, 500, 550, 600, 800, 1000, 2000)) {
  set.seed(5678)
  print(nt)
  rf_maxtrees <- train(post_hyp_binary~., data = summary_windowsTrain2,
                       method = "rf", 
                       trControl = trainControl(method = "cv", number = 5),
                       ntree = nt)
  key <- toString(nt)
  store_maxtrees[[key]] <- rf_maxtrees
}
results_tree <- resamples(store_maxtrees)
summary(results_tree)
# best n_t = 600
# tuning maxnodes
store_maxnode <- list()
for (mn in c(3, 6, 9, 12, 15, 18, 21, 24, 27, 30)) {
  set.seed(1234)
  print(mn)
  rf_maxnode <- train(post_hyp_binary~., data = summary_windowsTrain2,
                       method = "rf", 
                       trControl = trainControl(method = "cv", number = 5),
                       maxnodes = mn)
  current_iteration <- toString(mn)
  store_maxnode[[current_iteration]] <- rf_maxnode
}
results_node <- resamples(store_maxnode)
summary(results_node)
# best number of nodes = 24
# Tuning plots ----
results_summary <- summary(results_tree)$statistics$Accuracy[, "Mean"]
results_df <- data.frame(val = as.numeric(names(results_summary)),
                              Accuracy = results_summary)
ggplot(results_df, aes(x = val, y = Accuracy)) +
  geom_line(color = "black", linewidth = 0.75) +
  geom_point(color = "black", size = 1) +
  labs(title = "Effect of number of trees on model accuracy",
       x = "Number of trees", y = "Mean Accuracy") + theme_minimal()
results_summary <- summary(results_node)$statistics$Accuracy[, "Mean"]
results_df <- data.frame(val = as.numeric(names(results_summary)),
                              Accuracy = results_summary)
ggplot(results_df, aes(x = val, y = Accuracy)) +
  geom_line(color = "black", linewidth = 0.75) +
  geom_point(color = "black", size = 1) +
  labs(title = "Effect of number of maximum nodes on model accuracy",
       x = "Number of maximum nodes", y = "Mean Accuracy") + theme_minimal()
# Final Random forest model ----
longest_cols <- grep("^longest_hypoxic_seq", colnames(summary_windowsTrain), value = TRUE)
columns <- c(colnames(summary_windowsTrain)[8:25], "post_hyp_binary", longest_cols)
summary_windowsTrain2 <- summary_windowsTrain %>% group_by(post_hyp_binary) %>%
  slice_sample(n = 2*min(table(summary_windowsTrain$post_hyp_binary))) 
summary_windowsTrain2 <- summary_windowsTrain2[,columns]
summary_windowsTest2 <- summary_windowsTest %>% group_by(post_hyp_binary) %>%
  slice_sample(n = 2*min(table(summary_windowsTest$post_hyp_binary)))
summary_windowsTest2 <- summary_windowsTest2[,columns]
model_5 = train(post_hyp_binary~., data = summary_windowsTrain2,
                method = "rf", 
                trControl = trainControl(method = "cv", number = 5),
                ntree = 600, maxnodes = 22)
# model_6 = train(post_hyp_binary~., data = summary_windowsTrain2,
#                 method = "rf", 
#                 trControl = trainControl(method = "cv", number = 5),
#                 ntree = 600, maxnodes = 22)
pred <- predict(model_5, newdata = summary_windowsTest2)
confusionMatrix(pred, summary_windowsTest2$post_hyp_binary)
confusionMatrix(data = pred, reference = summary_windowsTest2$post_hyp_binary, mode = "prec_recall")
dim(summary_windowsTrain2)
dim(summary_windowsTest2)
# calculating AUPRC
pred_probs <- predict(model_5, newdata = summary_windowsTest2, type = "prob")
positive_probs <- pred_probs[, "TRUE"]  
true_labels <- ifelse(summary_windowsTest2$post_hyp_binary == "TRUE", 1, 0)
pr <- pr.curve(scores.class0 = positive_probs[true_labels == "1"],
                 scores.class1 = positive_probs[true_labels == "0"],
                 curve = TRUE)
importance <- varImp(model_5, scale = FALSE)
rownames(importance[["importance"]]) <- gsub("First_", "", rownames(importance[["importance"]]))
plot(importance, col = "deepskyblue4")
# testing biased model on unbiased model
longest_cols <- grep("^longest_hypoxic_seq", colnames(summary_windowsTrain), value = TRUE)
columns <- c(colnames(summary_windowsTrain)[8:25], "post_hyp_binary", longest_cols)
summary_windowsTest2 <- summary_windowsTest[,columns]
pred <- predict(model_5, newdata = summary_windowsTest2)
confusionMatrix(pred, summary_windowsTest2$post_hyp_binary)
confusionMatrix(data = pred, reference = summary_windowsTest2$post_hyp_binary, mode = "prec_recall")
dim(summary_windowsTrain2)
dim(summary_windowsTest2)
# calculating AUPRC
pred_probs <- predict(model_5, newdata = summary_windowsTest2, type = "prob")
positive_probs <- pred_probs[, "TRUE"]  
true_labels <- ifelse(summary_windowsTest2$post_hyp_binary == "TRUE", 1, 0)
pr <- pr.curve(scores.class0 = positive_probs[true_labels == "1"],
               scores.class1 = positive_probs[true_labels == "0"],
               curve = TRUE)
importance <- varImp(model_5, scale = FALSE)
plot(importance)
rownames(importance[["importance"]]) <- gsub("First_", "", rownames(importance[["importance"]]))
plot(importance, main = "RF, caret w/ 5 cv, binary, 2:1, first features + longest seq, biased tested on unbiased")
# calculating ROC values
pred_probs <- predict(model_5, newdata = summary_windowsTest2, type = "prob")
roc_curve <- roc(summary_windowsTest2$post_hyp_binary, pred_probs[,2])
auc_val <- auc(roc_curve)
plot(roc_curve, col = "deepskyblue4")
text(x = 0.25, y = 0.25, 
     labels = paste0("AUC = ", round(auc_val, 3)), 
     col = "black", cex = 1)
# Patient level accuracy calculations
summary_windowsTest2 <- summary_windowsTest
summary_windowsTest2$pred <- predict(model_5, newdata = summary_windowsTest2)
summary_windowsTest2 <- summary_windowsTest2[, c("patientID", "pred", "post_hyp_binary")]
patient_accuracy <- summary_windowsTest2 %>%
  group_by(patientID) %>% summarise(accuracy = mean(pred == post_hyp_binary))
patient_accuracy <- patient_accuracy %>% left_join(summary_windowsTest2 %>%
      group_by(patientID) %>% summarise(n_windows = n()), by = "patientID")
#plotting
ggplot(patient_accuracy, aes(x = reorder(factor(patientID), accuracy), y = accuracy, fill = n_windows)) +
  geom_col() +
  scale_fill_viridis_c(option = "plasma") +  # nice continuous color scale
  labs(title = "Patient-level Accuracy with Number of Windows",
       x = "Patient ID (sorted by accuracy)", y = "Accuracy", fill = "Number of Windows") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))
#calculating alarm burden (number of false positives over total positives)
report_alarm_burden <- function(cm) {
  TP <- cm$table["TRUE", "TRUE"]
  FP <- cm$table["TRUE", "FALSE"]
  total_alarms <- TP + FP
  alarm_burden <- ifelse(total_alarms > 0, FP / total_alarms, NA)
  return(alarm_burden)
}
report_alarm_burden(cm)
#calculating alarm burden (number of positives/total number)
perc_positive <- function(cm) {
  TP <- cm$table["TRUE", "TRUE"]
  FP <- cm$table["TRUE", "FALSE"]
  total_predictions <- sum(cm$table)
  num_predicted_positives <- TP + FP
  alarm_burden <- num_predicted_positives / total_predictions
  return(alarm_burden)
}
# False positive windows ----
# Identify false positive indices
false_positive_indices <- which(pred == TRUE & summary_windowsTest2$post_hyp_binary == FALSE)
false_positive_windows <- summary_windowsTest[false_positive_indices, ]
# plotting patient counts
patient_counts <- as.data.frame(table(false_positive_windows$patientID))
colnames(patient_counts) <- c("PatientID", "numFalsePositives")
patient_counts <- patient_counts %>%
  arrange(desc(numFalsePositives)) %>%
  mutate(CumulativePercent = cumsum(numFalsePositives) / sum(numFalsePositives) * 100)
# Pareto Chart
ggplot(patient_counts, aes(x = reorder(PatientID, -numFalsePositives), y = numFalsePositives)) +
  geom_bar(stat = "identity", fill = "blue", alpha = 0.7) +
  geom_line(aes(y = CumulativePercent * max(numFalsePositives) / 100, group = 1), color = "red", size = 1) +
  geom_point(aes(y = CumulativePercent * max(numFalsePositives) / 100), color = "red") +
  scale_y_continuous(sec.axis = sec_axis(~ . * 100 / max(patient_counts$numFalsePositives), name = "Cumulative Percentage")) +
  labs(title = "Pareto Chart of False Positives per Patient", x = "Patient ID", y = "False Positives") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_flip()
# calculating proportion of false positive windows per patient
false_positive_proportion <- aggregate(time ~ patientID, data = false_positive_windows, FUN = length)
colnames(false_positive_proportion) <- c("patientID", "false_windows")
total_count <- aggregate(time ~ patientID, data = summary_windowsTest, FUN = length)
colnames(total_count) <- c("patientID", "total_windows")
false_positive_proportion <- merge(total_count, false_positive_proportion, by = "patientID", all.x = TRUE)
false_positive_proportion$false_windows[is.na(false_positive_proportion$false_windows)] <- 0
false_positive_proportion$proportion <- false_positive_proportion$false_windows / false_positive_proportion$total_windows
rm(total_count)
false_positive_proportion$randomization_number <- names(pVisi_data)[false_positive_proportion$patientID]
false_positive_proportion <- merge(false_positive_proportion, baseline[,c(1:8, 60:69)])
# Create Pareto Chart with proportion of false windows
false_positive_proportion_small <- false_positive_proportion[false_positive_proportion$false_windows > 0, ]
ggplot(false_positive_proportion_small, aes(x = reorder(patientID, -false_windows), y = false_windows)) +
  geom_bar(stat = "identity", aes(fill = age), alpha = 0.7) +
  geom_line(aes(y = proportion * max(false_positive_proportion_small$false_windows), group = 1), color = "red", size = 1) +
  geom_point(aes(y = proportion * max(false_positive_proportion_small$false_windows)), color = "red") +
  scale_y_continuous(
    sec.axis = sec_axis(~ . * 100 / max(false_positive_proportion_small$false_windows), 
                        name = "% of false windows over total windows")) +
  geom_point(aes(y = false_windows + 5, color = factor(gender)), 
    size = 2, show.legend = TRUE) + 
  scale_color_manual(values = c("pink2", "blue"), labels = c("Female", "Male")) +
  labs( title = "Pareto Chart of False Positives per Patient",
    x = "Patient ID", y = "Number of False Positives") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_flip()
ggplot(false_positive_proportion_small, aes(x = reorder(patientID, -false_windows), y = false_windows)) +
  geom_bar(stat = "identity", aes(fill = bmi), alpha = 0.7) +
  geom_line(aes(y = proportion * max(false_positive_proportion_small$false_windows), group = 1), color = "red", size = 1) +
  geom_point(aes(y = proportion * max(false_positive_proportion_small$false_windows)), color = "red") +
  scale_y_continuous(
    sec.axis = sec_axis(~ . * 100 / max(false_positive_proportion_small$false_windows), 
                        name = "% of false windows over total windows")) +
  geom_point(aes(y = false_windows + 5, color = factor(stopbang_risk)), 
             size = 2, show.legend = TRUE) + 
  scale_color_manual(values = c("orange3", "purple2"), labels = c("High Risk", "Low Risk")) +
  labs( title = "Pareto Chart of False Positives per Patient",
        x = "Patient ID", y = "Number of False Positives") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_flip()
ggplot(false_positive_proportion_small, aes(x = reorder(patientID, -false_windows), y = false_windows)) +
  geom_bar(stat = "identity", aes(fill = bmi), alpha = 0.7) +
  geom_line(aes(y = proportion * max(false_positive_proportion_small$false_windows), group = 1), color = "red", size = 1) +
  geom_point(aes(y = proportion * max(false_positive_proportion_small$false_windows)), color = "red") +
  scale_y_continuous(
    sec.axis = sec_axis(~ . * 100 / max(false_positive_proportion_small$false_windows), 
                        name = "% of false windows over total windows")) +
  geom_point(aes(y = false_windows + 5, color = factor(stopbang_risk)), 
             size = 2, show.legend = TRUE) + 
  scale_color_manual(values = c("orange3", "purple2"), labels = c("High Risk", "Low Risk")) +
  labs( title = "Pareto Chart of False Positives per Patient",
        x = "Patient ID", y = "Number of False Positives") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_flip()
ggplot(false_positive_proportion_small, aes(x = reorder(patientID, -false_windows), y = false_windows)) +
  geom_bar(stat = "identity", aes(fill = factor(stopbang_risk)), alpha = 0.7) +
  scale_fill_manual(values = c("orange3", "purple2"), name = "STOP-BANG Risk", labels = c("High Risk", "Low Risk")) +
  labs( title = "Pareto Chart of False Positives per Patient",
        x = "Patient ID", y = "Number of False Positives") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggplot(false_positive_proportion, aes(x = has_false_positive, y = bmi, fill = has_false_positive)) +
  geom_boxplot(alpha = 0.7) +
  labs(
    title = "BMI Comparison by False Positive Status",
    x = "False Positive Status",
    y = "BMI"
  ) +
  theme_minimal() +
  theme(legend.position = "none")
# running chi-squared test
false_positive_proportion$has_false_positive <- ifelse(false_positive_proportion$false_windows > 0, "Has false", "No False")
chi1 <- table(false_positive_proportion$has_false_positive, false_positive_proportion$gender)
chi1_test <- chisq.test(chi1)
medianAge <- median(as.numeric(baseline$age), na.rm = TRUE)
medianBMI <- median(as.numeric(baseline$bmi), na.rm = TRUE)
false_positive_proportion$age_cat <- ifelse(false_positive_proportion$age >= medianAge, "Old", "Young")
false_positive_proportion$bmi_cat <- ifelse(false_positive_proportion$age >= medianBMI, "highBMI", "lowBMI")
chi2 <- table(false_positive_proportion$has_false_positive, false_positive_proportion$bmi_cat)
chi2_test <- chisq.test(chi2)
chi3 <- table(false_positive_proportion$has_false_positive, false_positive_proportion$stopbang_risk)
chi3_test <- chisq.test(chi3)

mosaicplot(table(false_positive_proportion$has_false_positive, false_positive_proportion$stopbang_risk), 
           main = "Mosaic Plot of Variable1 vs Variable2", 
           color = TRUE)
mosaicplot(table(false_positive_proportion$has_false_positive, false_positive_proportion$gender), 
           main = "Mosaic Plot of Variable1 vs Variable2", 
           color = TRUE)

# gaussian false positive ----
gaussianFalsePositive <- left_join(false_positive_proportion_small, gaussianG,
                         by = c("randomization_number" = "ID"))
# testing on disjoint windows dataset ----
columns <- c(colnames(summary_windowsTrain)[8:25], "post_hyp_binary", "longest_hypoxic_seq")
summary_windows_disjoint <- summary_windows_disjoint %>% 
  mutate(post_hyp_binary = case_when(post_hyp >= 0.9 ~ "TRUE", TRUE ~ "FALSE"))
summary_windows_disjoint$post_hyp_binary <- factor(summary_windows_disjoint$post_hyp_binary, levels = c(TRUE, FALSE))
summary_windowsTestDisjoint <- summary_windows_disjoint[,columns]
summary_windowsTestDisjoint <- summary_windows_disjoint[summary_windows_disjoint$patientID %in% testSet,]

dim(summary_windowsTestDisjoint)
summary_windowsTestDisjoint <- na.omit(summary_windowsTestDisjoint) #removing NAs
dim(summary_windowsTestDisjoint)

pred <- predict(model_5, newdata = summary_windowsTestDisjoint)
confusionMatrix(pred, summary_windowsTestDisjoint$post_hyp_binary)
confusionMatrix(data = pred, reference = summary_windowsTestDisjoint$post_hyp_binary, mode = "prec_recall")
dim(summary_windowsTrain2)
dim(summary_windowsTest2)
importance <- varImp(model_5, scale = FALSE)
plot(importance)
rownames(importance[["importance"]]) <- gsub("First_", "", rownames(importance[["importance"]]))
plot(importance, main = "RF, caret w/ 5 cv, binary, 2:1, first features + longest seq, biased tested on unbiased")
# testing 