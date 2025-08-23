# libraries ----
library(ggplot2)
library(dplyr)
#exploration of missing data ----
missing <- data.frame() 
for (j in 1:length(pVisi_data)) {
  print(j)
  dataJ <- pVisi_data[[j]]
  missing <- rbind(missing, dataJ[dataJ$timeDiff >= 60, c("StudyID", "timeDiff")])
}
#plotting ----
#histogram plot of missing data
ggplot(missing, aes(x = timeDiff)) +
  geom_histogram(binwidth = 15, fill = "deepskyblue4", color = "white") +
  labs(title = "Histogram of Time Difference", x = "Length of gap (seconds)", 
       y = "Count") +
  theme_minimal() + xlim(0, 500)
#proportion of missing data per person
missing_count_per_person <- missing %>%
  group_by(StudyID) %>%
  dplyr::summarise(total_missing = sum(timeDiff), na.rm = TRUE)
missing_count_per_person$total_missing_min <- missing_count_per_person$total_missing/60
ggplot(missing_count_per_person, aes(x = total_missing_min)) +
  geom_histogram(fill = "deepskyblue4", color = "white") +
  labs(title = "Distribution of Missing Data per Person", 
       x = "Total missing data in minutes", 
       y = "Frequency of Patients") +
  theme_minimal()

#exploring interventions by nurses
data_postOp <- data[, sort(unique(c(1:4, grep("^post_op", colnames(data)))))]
data_postOp <- data_postOp[data_postOp$redcap_event_name %in% c("pod_0_after_pacu_arm_1", "device_summary_arm_1"), ]
data_postOp <- data_postOp %>%
  group_by(uid) %>%
  mutate(study_id = ifelse(study_id == "", first(study_id[study_id != ""]), study_id)) %>%
  ungroup()
data_postOp <- data_postOp[data_postOp$redcap_event_name %in% c("pod_0_after_pacu_arm_1"), ]
