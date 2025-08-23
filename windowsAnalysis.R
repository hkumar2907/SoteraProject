# libraries ----
library(dplyr)
library(plyr)
library(forecast)
library(randomForest)
library(xts)
library(depmixS4)
library(mclust)
library(zoo)
library(psych)
library(moments)
library(tidyr)
library(ggpubr)
library(GGally)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(rstatix)
# generating windows 2 ----
windows2 <- list()
for (j in 1:length(pVisi_data)) {
  print(j)
  dataJ <- pVisi_data[[j]]
  dataJ$timeDiff <- c(0, as.numeric(difftime(dataJ$dev_reading_tm[-1], dataJ$dev_reading_tm[-nrow(dataJ)], units = "secs")))
  dataJ$hypoxicPeriods <- dataJ$SPO2 < 90
  dataJ$cumulativeTime <- cumsum(dataJ$timeDiff)
  # Initialize rolling hypoxia indicator
  dataJ$rollingHypoxia <- rep(0, nrow(dataJ))
  # Rolling window logic based on time intervals
  i <- 121
  counter <- 1
  while (i <= nrow(dataJ)) {
    # Calculate the start time of the 30-minute window ending at row i
    windowStartTime <- dataJ$cumulativeTime[i] - 30 * 60  # 30 minutes in seconds
    # Identify rows within the 30-minute window
    windowIndices <- which(dataJ$cumulativeTime >= windowStartTime & dataJ$cumulativeTime <= dataJ$cumulativeTime[i])
    # Calculate total hypoxic time in the window
    hypoxicTime <- sum(dataJ$timeDiff[windowIndices] * dataJ$hypoxicPeriods[windowIndices])
    # Mark as hypoxic if the threshold is met
    if ((hypoxicTime >= 1620) && (dataJ$timeDiff[i] < 60)) { # 27 minutes (90%) * 60 = 1620 seconds
      dataJ$rollingHypoxia[i] <- counter
      # Skip 241 rows (adjust for time intervals)
      i <- i + 241
      counter <- counter + 1
    } 
    else {
      # Otherwise, just move to the next row
      i <- i + 1
    }
  }
  # Filter and group valid hypoxia windows
  if (nrow(dataJ %>% filter(rollingHypoxia != 0)) != 0) {
    windows2[[j]] <- dataJ %>% 
      filter(rollingHypoxia != 0) %>% 
      group_by(rollingHypoxia) %>% 
      summarize(end = dev_reading_tm,
                pre_time = dev_reading_tm - as.difftime(60, units = "mins"), 
                start = dev_reading_tm - as.difftime(30, units = "mins")) %>%
      mutate(person = j)
  }
  else {
    windows2[[j]] <- data.frame()
  }
}

# generating windows  ----
windows <- list()
for (j in 1:length(pVisi_data)) {
  print(j)
  dataJ <- pVisi_data[[j]]
  dataJ$timeDiff <- c(0, as.numeric(difftime(dataJ$dev_reading_tm[-1], dataJ$dev_reading_tm[-nrow(dataJ)], units = "secs")))
  dataJ$hypoxicPeriods <- dataJ$SPO2 < 90
  rleHypoxia <- rle(dataJ$hypoxicPeriods)
  # Create a vector to store the cumulative time for each run
  cumulativeTimeVector <- rep(NA, length(dataJ$hypoxicPeriods))
  start <- 1
  dataJ$validHypoxia <- rep(0, length(dataJ$hypoxicPeriods))
  counter <- 1
  for (i in seq_along(rleHypoxia$lengths)) {
    if (rleHypoxia$values[i]) { 
      end <- start + rleHypoxia$lengths[i] - 1
      cumulativeTimeVector[start:end] <- cumsum(dataJ$timeDiff[start:end])
      if (cumulativeTimeVector[end] >= 1800) {
        dataJ$validHypoxia[end] <- counter
        counter <- counter + 1
      }
    }
    start <- start + rleHypoxia$lengths[i]
  }
  dataJ$validHypoxia <- as.factor( dataJ$validHypoxia)
  if (nrow(dataJ %>% filter(validHypoxia != 0)) != 0) {
    windows[[j]] <- dataJ %>% filter(validHypoxia != 0) %>% group_by(validHypoxia) %>% summarize(end = dev_reading_tm,duration = cumulativeTime,pre_time = dev_reading_tm - as.difftime(30, units = "mins") - cumulativeTime, start = dev_reading_tm - cumulativeTime) %>% mutate(person = j)
    # windows[[j]] <- dataJ %>% filter(validHypoxia != 0) %>% group_by(validHypoxia) %>%
    #   summarize(end = dev_reading_tm,duration = cumulativeTime, 
    #             pre_time = dev_reading_tm - as.difftime(60, units = "mins"),
    #             start = dev_reading_tm - as.difftime(30, units = "mins"), .groups = "drop") %>% mutate(person = j)
  }
  else {
    windows[[j]] <- data.frame()
  }
  pVisi_data[[j]] <- dataJ
}
# generating windows 3 ----
windows3 <- list()
windows3nonHypoxia <- list()
for (j in 1:length(pVisi_data)) {
  print(j)
  dataJ <- pVisi_data[[j]]
  dataJ$timeDiff <- c(0, as.numeric(difftime(dataJ$dev_reading_tm[-1], dataJ$dev_reading_tm[-nrow(dataJ)], units = "secs")))
  dataJ$hypoxicPeriods <- dataJ$SPO2 < 90
  dataJ$cumulativeTime <- cumsum(dataJ$timeDiff)
  rleHypoxia <- rle(dataJ$hypoxicPeriods)
  # Create a vector to store the cumulative time for each run
  cumulativeTimeVector <- rep(NA, length(dataJ$hypoxicPeriods))
  nonCumulativeTimeVector <- rep(NA, length(dataJ$hypoxicPeriods))
  start <- 1
  dataJ$validHypoxia <- rep(0, length(dataJ$hypoxicPeriods))
  dataJ$nonHypoxia <- rep(0, length(dataJ$hypoxicPeriods))
  counter <- 1
  counter2 <- 1
  for (i in seq_along(rleHypoxia$lengths)) {
    if (rleHypoxia$values[i]) { 
      end <- start + rleHypoxia$lengths[i] - 1
      cumulativeTimeVector[start:end] <- cumsum(dataJ$timeDiff[start:end])
      if (cumulativeTimeVector[end] >= 900) {
        dataJ$validHypoxia[end] <- counter
        counter <- counter + 1
      } 
    }
    else { 
      end <- start + rleHypoxia$lengths[i] - 1
      nonCumulativeTimeVector[start:end] <- cumsum(dataJ$timeDiff[start:end])
      if (nonCumulativeTimeVector[end] >= 900) {
        dataJ$nonHypoxia[end] <- counter2
        counter2 <- counter2 + 1
      }
    }
    start <- start + rleHypoxia$lengths[i]
  }
  dataJ$validHypoxia <- as.factor(dataJ$validHypoxia)
  dataJ$nonHypoxia <- as.factor(dataJ$nonHypoxia)
  dataJ$cumulativeTimeVector <- cumulativeTimeVector
  dataJ$nonCumulativeTimeVector <- nonCumulativeTimeVector
  if (nrow(dataJ %>% filter(validHypoxia != 0)) != 0) {
    windows3[[j]] <- dataJ %>% filter(validHypoxia != 0) %>%
      group_by(validHypoxia) %>% 
      summarize(end = dev_reading_tm,
                duration = cumulativeTimeVector,
                pre_time = dev_reading_tm - as.difftime(30, units = "mins") - cumulativeTimeVector, 
                start = dev_reading_tm - cumulativeTimeVector) %>% mutate(person = j)
  }
  else {
    windows3[[j]] <- data.frame()
  }
  if (nrow(dataJ %>% filter(nonHypoxia != 0)) != 0) {
    windows3nonHypoxia[[j]] <- dataJ %>% filter(nonHypoxia != 0) %>% 
      group_by(nonHypoxia) %>% 
      summarize(end = dev_reading_tm,duration = nonCumulativeTimeVector,
                pre_time = dev_reading_tm - as.difftime(30, units = "mins") - nonCumulativeTimeVector, 
                start = dev_reading_tm - nonCumulativeTimeVector) %>% mutate(person = j)
  }
  else {
    windows3nonHypoxia[[j]] <- data.frame()
  }
  pVisi_data[[j]] <- dataJ
}

# Generating window features ----
functions <- list(mean, sd, median, max, min,
                    function(x) max(x) - min(x),
                    function(x) mean(abs(diff(x))),
                    function(x) median(abs(diff(x))),
                    IQR, 
                    function(x) sum(x > mean(x)),
                    function(x) sum(x > geometric.mean(x)), 
                    function(x) length(which(diff(sign(diff(x))) == 2)),
                    skewness, 
                    kurtosis, 
                    function(x) sum(x^2),
                    function(x) which.min(x),
                    function(x) which.max(x),
                    function(x) diff(c(which.max(x), which.min(x))))
# Names for each function
functionNames <- c("Mean", "Std", "Median", "Max", "Min", 
                     "MaxMinusMin", "MeanAbsoluteDiff", "MedianAbsoluteDiff",
                     "IQR", "ValsAboveMean", "ValsAboveGeoMean", "NumsPeaks",
                     "Skew", "Kurtosis", "Energy", "Argmin", "Argmax", "ArgDiff")
#windows 1
features <- list()
for (j in 1:length(pVisi_data)) {
    if (length(windows[[j]]) != 0) {
      dataX <- windows[[j]]
      features[[j]] <- data.frame(matrix(NA, nrow=nrow(windows[[j]]), ncol=length(functions)))
      colnames(features[[j]]) <- functionNames
      for (xc in 1:nrow(windows[[j]])) {
        x <- subset(pVisi_data[[j]], dev_reading_tm >= windows[[j]]$pre_time[xc] & dev_reading_tm <= windows[[j]]$start[xc])
        for (i in 1:length(functions)) {
          features[[j]][xc, i] <- functions[[i]](x$SPO2)
        }
      }
    }
  }
#windows 2
features2 <- list()
for (j in 1:length(pVisi_data)) {
    print(j)
    if (length(windows2[[j]]) != 0) {
      dataX <- windows2[[j]]
      features2[[j]] <- data.frame(matrix(NA, nrow=nrow(windows2[[j]]), ncol=length(functions)))
      colnames(features2[[j]]) <- functionNames
      for (xc in 1:nrow(windows2[[j]])) {
        x <- subset(pVisi_data[[j]], dev_reading_tm >= windows2[[j]]$pre_time[xc] & dev_reading_tm <= windows2[[j]]$start[xc])
        for (i in 1:length(functions)) {
          features2[[j]][xc, i] <- functions[[i]](x$SPO2)
        }
      }
    }
  }
#windows 3
features3 <- list()
  for (j in 1:length(pVisi_data)) {
    print(j)
    if (length(windows3[[j]]) != 0) {
      dataX <- windows3[[j]]
      features3[[j]] <- data.frame(matrix(NA, nrow=nrow(windows3[[j]]), ncol=length(functions)))
      colnames(features3[[j]]) <- functionNames
      for (xc in 1:nrow(windows3[[j]])) {
        x <- subset(pVisi_data[[j]], dev_reading_tm >= windows3[[j]]$pre_time[xc] & dev_reading_tm <= windows3[[j]]$start[xc])
        for (i in 1:length(functions)) {
          features3[[j]][xc, i] <- functions[[i]](x$SPO2)
        }
      }
    }
  }
# windows 3 non hypoxia ----
  features3NH <- list()
  features3 <- list()
  for (j in 1:length(pVisi_data)) {
    print(j)
    if ((length(windows3nonHypoxia[[j]]) && (length(windows3[[j]]) != 0)) != 0) {
      if (nrow(windows3[[j]]) <= nrow(windows3nonHypoxia[[j]])) {
        dataX <- windows3[[j]]
        features3[[j]] <- data.frame(matrix(NA, nrow=nrow(windows3[[j]]), ncol=length(functions)))
        colnames(features3[[j]]) <- functionNames
        for (xc in 1:nrow(windows3[[j]])) {
          x <- subset(pVisi_data[[j]], dev_reading_tm >= windows3[[j]]$pre_time[xc] & dev_reading_tm <= windows3[[j]]$start[xc])
          for (i in 1:length(functions)) {
            features3[[j]][xc, i] <- functions[[i]](x$SPO2)
          }
        }
        dataX <- windows3nonHypoxia[[j]]
        features3NH[[j]] <- data.frame(matrix(NA, nrow=nrow(windows3[[j]]), ncol=length(functions)))
        sampled_indices <- sample(1:nrow(windows3nonHypoxia[[j]]), nrow(windows3[[j]]))
        sampled_windows <- windows3nonHypoxia[[j]][sampled_indices, ]
        colnames(features3NH[[j]]) <- functionNames
        for (xc in 1:nrow(sampled_windows)) {
          x <- subset(pVisi_data[[j]], dev_reading_tm >= sampled_windows$pre_time[xc] & dev_reading_tm <= sampled_windows$start[xc])
          for (i in 1:length(functions)) {
            features3NH[[j]][xc, i] <- functions[[i]](x$SPO2)
          }
        }
      }
      else {
        dataX <- windows3[[j]]
        sampled_indices <- sample(1:nrow(windows3[[j]]), nrow(windows3nonHypoxia[[j]]))
        sampled_windows <- windows3[[j]][sampled_indices, ]
        features3[[j]] <- data.frame(matrix(NA, nrow=nrow(windows3nonHypoxia[[j]]), ncol=length(functions)))
        colnames(features3[[j]]) <- functionNames
        for (xc in 1:nrow(sampled_windows)) {
          x <- subset(pVisi_data[[j]], dev_reading_tm >= sampled_windows$pre_time[xc] & dev_reading_tm <= sampled_windows$start[xc])
          for (i in 1:length(functions)) {
            features3[[j]][xc, i] <- functions[[i]](x$SPO2)
          }
        }
        dataX <- windows3nonHypoxia[[j]]
        features3NH[[j]] <- data.frame(matrix(NA, nrow=nrow(windows3nonHypoxia[[j]]), ncol=length(functions)))
        colnames(features3NH[[j]]) <- functionNames
        for (xc in 1:nrow(windows3nonHypoxia[[j]])) {
          x <- subset(pVisi_data[[j]], dev_reading_tm >= windows3nonHypoxia[[j]]$pre_time[xc] & dev_reading_tm <= windows3nonHypoxia[[j]]$start[xc])
          for (i in 1:length(functions)) {
            features3NH[[j]][xc, i] <- functions[[i]](x$SPO2)
          }
        }
      }
    }
  }
  
# plotting ----
  features3_combined <- bind_rows(features3, .id = "Patient") %>%
    mutate(Windows = "Hypoxia")
  features3NH_combined <- bind_rows(features3NH, .id = "Patient") %>%
    mutate(Windows = "Non-Hypoxia")
  features_combined <- bind_rows(features3_combined, features3NH_combined)
  features_combined <- features_combined %>%
    pivot_longer(
      cols = -c(Patient, Windows),
      names_to = "Measure",
      values_to = "Value"
    )
  # Create the plot with 18 facets
  ggplot(features_combined, aes(x = Windows, y = Value, fill = Windows)) +
    geom_boxplot(outlier.alpha = 0.5) +
    facet_wrap(~ Measure, scales = "free", ncol = 6) +  # 6 columns for better visualization
    theme_minimal(base_size = 14) +
    theme(legend.position = "top", strip.text = element_text(size = 10, face = "bold")) +
    labs(title = "Comparative Measures for 30-Minute Periods Before Hypoxia and Normoxia Windows",
         x = "Windows",y = "Values") +
    scale_fill_manual(values = c("Hypoxia" = "lightcoral", "Non-Hypoxia" = "lightseagreen"), 
                     labels = c("Pre-hypoxia window", "Pre-normoxia window")) +
    stat_compare_means(aes(group = Windows), method = "t.test", label = "p.signif",
                       label.x.npc = 'center')


# windows 4 (rowan's code )----
#windows4 <- list()
#for (j in 1:length(pVisi_data)) {
    print(j)
    dataJ <- pVisi_data[[j]]
    dataJ <- dataJ[dataJ$SPO2 != "XX",] #removing all XXs in the SPO2
    dataJ$SPO2 <- as.numeric(dataJ$SPO2) #turning spo2 is numeric
    dataJ <- dataJ %>% drop_na(SPO2)
    dataJ$dev_reading_tm = as.POSIXct(dataJ$dev_reading_tm, format="%d%b%Y:%T", tz="EST")
    dataJ <- dataJ[order(dataJ$dev_reading_tm),]
    dataJ$timeDiff <- c(0, as.numeric(difftime(dataJ$dev_reading_tm[-1], dataJ$dev_reading_tm[-nrow(dataJ)], units = "secs")))
    dataJ$hypoxic <- dataJ$SPO2 < 90
    dataJ$hypoxicPeriods <- dataJ$SPO2 < 90
    #dataJ <- dataJ[dataJ$timeDiff < 90,]
    dataJ$cumulativeTime <- cumsum(dataJ$timeDiff)
    nr <- nrow(dataJ)
    windowSize <- 60*60/15
    windows4[[j]] = lapply(seq(1, nr - windowSize + 1, by = windowSize / 10), function(i) dataJ[i:(i + windowSize - 1), ])
  }  
# windows function
window_overlap <- function(pVisi_data, overlap_percent = 50) {
    results <- list()
    window_minutes <- 60
    interval_seconds <- 15
    window_size <- (window_minutes * 60) / interval_seconds
    step_size <- window_size * (1 - overlap_percent / 100)
    
    for (j in seq_along(pVisi_data)) {
      if (j == 461) next
      print(j)
      dataJ <- pVisi_data[[j]]
      
      # Clean and preprocess data
      dataJ <- dataJ[dataJ$SPO2 != "XX", ]
      dataJ$SPO2 <- as.numeric(dataJ$SPO2)
      dataJ <- dataJ %>% drop_na(SPO2)
      dataJ$dev_reading_tm <- as.POSIXct(dataJ$dev_reading_tm, format = "%d%b%Y:%T", tz = "EST")
      dataJ <- dataJ[order(dataJ$dev_reading_tm), ]
      dataJ$timeDiff <- c(0, diff(as.numeric(dataJ$dev_reading_tm)))
      dataJ$hypoxic <- dataJ$SPO2 < 90
      dataJ$hypoxicPeriods <- dataJ$SPO2 < 90
      dataJ$cumulativeTime <- cumsum(dataJ$timeDiff)
      
      nr <- nrow(dataJ)
      idx_seq <- seq(1, nr - window_size + 1, by = step_size)
      
      results[[j]] <- lapply(idx_seq, function(i) {
        dataJ[i:(i + window_size - 1), ]
      })
    }
    return(results)
  }
windows4_10percOverlap <- window_overlap(pVisi_data, overlap_percent = 10)
windows4_30percOverlap <- window_overlap(pVisi_data, overlap_percent = 30)
windows4_50percOverlap <- window_overlap(pVisi_data, overlap_percent = 50)
windows4_70percOverlap <- window_overlap(pVisi_data, overlap_percent = 70)
windows4_90percOverlap <- window_overlap(pVisi_data, overlap_percent = 90)

# Summary matrices (rowan's code) ----
summarize_patient_windows <- function(windows_list, functions, functionNames) {
  summaries <- list()
  get_longest_seq <- function(x, threshold) {
    hypoxicThreshold <- x < threshold
    hypoxic_sequence <- rle(hypoxicThreshold)
    longest_seq <- max(hypoxic_sequence$lengths[hypoxic_sequence$values == TRUE], na.rm = TRUE)
    return(ifelse(is.finite(longest_seq), longest_seq, 0))
  }
  
  for (pid in seq_along(windows_list)) {
    if (pid == 461) next
    print(pid)
    features1 <- list()
    features2 <- list()
    longest_sequences <- list()
    longest_sequences_88 <- list()
    longest_sequences_89 <- list()
    longest_sequences_91 <- list()
    longest_sequences_92 <- list()
    maximumGap <- list()
    
    for (i in seq_along(windows_list[[pid]])) {
      current_window <- windows_list[[pid]][[i]]
      first_half <- current_window[1:(nrow(current_window) %/% 2), ]
      second_half <- current_window[(nrow(current_window) %/% 2 + 1):nrow(current_window), ]
      
      # if (segment == "first10") {
      #   selected_segment <- first_half[first_half$cumtime <= 600, ]
      # } else if (segment == "last10") {
      #   selected_segment <- first_half[first_half$cumtime >= (total_time_first_half - 600), ]
      # } else if (segment == "middle10") {
      #   selected_segment <- first_half[first_half$cumtime >= (total_time_first_half / 2 - 300) &
      #                                    first_half$cumtime <= (total_time_first_half / 2 + 300), ]
      # } else if (segment == "first20") {
      #   selected_segment <- first_half[first_half$cumtime <= 1200, ]
      # } else if (segment == "last20") {
      #   selected_segment <- first_half[first_half$cumtime >= (total_time_first_half - 1200), ]
      # } else if (segment == "middle20") {
      #   selected_segment <- first_half[first_half$cumtime >= (total_time_first_half / 2 - 600) &
      #                                    first_half$cumtime <= (total_time_first_half / 2 + 600), ]
      # } else {
      #   selected_segment <- first_half  # fallback: full first half
      # }
      
      features1[[i]] <- sapply(functions, function(f) f(first_half$SPO2))
      features2[[i]] <- sapply(functions, function(f) f(second_half$SPO2))
      
      longest_sequences[[i]]     <- get_longest_seq(current_window$SPO2, 90)
      longest_sequences_88[[i]]  <- get_longest_seq(current_window$SPO2, 88)
      longest_sequences_89[[i]]  <- get_longest_seq(current_window$SPO2, 89)
      longest_sequences_91[[i]]  <- get_longest_seq(current_window$SPO2, 91)
      longest_sequences_92[[i]]  <- get_longest_seq(current_window$SPO2, 92)
      maximumGap[[i]]            <- max(current_window$timeDiff)
    }
    
    features1 <- do.call(rbind, features1)
    features2 <- do.call(rbind, features2)
    colnames(features1) <- paste0("First_", functionNames)
    colnames(features2) <- paste0("Second_", functionNames)
    
    windows_hypoxia <- sapply(windows_list[[pid]], function(i) sum(i$timeDiff * i$hypoxicPeriods))
    windows_mean    <- sapply(windows_list[[pid]], function(i) mean(i$SPO2))
    windows_median  <- sapply(windows_list[[pid]], function(i) median(i$SPO2))
    windows_continuous <- sapply(windows_list[[pid]], function(i) any(i$timeDiff > 60))
    
    windows_hyp_1st <- sapply(windows_list[[pid]], function(i) {
      first_half <- i[1:(nrow(i) %/% 2), ]
      sum(first_half$SPO2 < 90) / (0.5 * nrow(i))
    })
    windows_hyp_2nd <- sapply(windows_list[[pid]], function(i) {
      second_half <- i[(nrow(i) %/% 2 + 1):nrow(i), ]
      sum(second_half$SPO2 < 90) / (0.5 * nrow(i))
    })
    df <- data.frame(
      time = windows_hypoxia,
      disjoint = windows_continuous,
      oxygen_mean = windows_mean,
      oxygen_median = windows_median,
      pre_hyp = windows_hyp_1st,
      post_hyp = windows_hyp_2nd,
      patientID = pid
    )
    df <- cbind(df, features1, features2,
                longest_hypoxic_seq = unlist(longest_sequences),
                longest_hypoxic_seq_88 = unlist(longest_sequences_88),
                longest_hypoxic_seq_89 = unlist(longest_sequences_89),
                longest_hypoxic_seq_91 = unlist(longest_sequences_91),
                longest_hypoxic_seq_92 = unlist(longest_sequences_92),
                maximumGap = unlist(maximumGap))
    df$meandiff <- c(0, diff(df$oxygen_mean))
    df$mediandiff <- c(0, diff(df$oxygen_median))
    summaries[[pid]] <- df
  }
  return(summaries)
}
summarize_patient_windows <- function(windows_list, functions, functionNames) {
  summaries <- list()
  get_longest_seq <- function(x, threshold) {
    hypoxicThreshold <- x < threshold
    hypoxic_sequence <- rle(hypoxicThreshold)
    longest_seq <- max(hypoxic_sequence$lengths[hypoxic_sequence$values == TRUE], na.rm = TRUE)
    return(ifelse(is.finite(longest_seq), longest_seq, 0))
  }
  
  for (pid in seq_along(windows_list)) {
    if (pid == 461) next
    print(pid)
    features1 <- list()
    features2 <- list()
    longest_sequences <- list()
    longest_sequences_88 <- list()
    longest_sequences_89 <- list()
    longest_sequences_91 <- list()
    longest_sequences_92 <- list()
    maximumGap <- list()
    
    for (i in seq_along(windows_list[[pid]])) {
      current_window <- windows_list[[pid]][[i]]
      first_half <- current_window[1:(nrow(current_window) %/% 2), ]
      second_half <- current_window[(nrow(current_window) %/% 2 + 1):nrow(current_window), ]
      
      features1[[i]] <- sapply(functions, function(f) f(first_half$SPO2))
      features2[[i]] <- sapply(functions, function(f) f(second_half$SPO2))
      
      longest_sequences[[i]]     <- get_longest_seq(current_window$SPO2, 90)
      longest_sequences_88[[i]]  <- get_longest_seq(current_window$SPO2, 88)
      longest_sequences_89[[i]]  <- get_longest_seq(current_window$SPO2, 89)
      longest_sequences_91[[i]]  <- get_longest_seq(current_window$SPO2, 91)
      longest_sequences_92[[i]]  <- get_longest_seq(current_window$SPO2, 92)
      maximumGap[[i]]            <- max(current_window$timeDiff)
    }
    
    features1 <- do.call(rbind, features1)
    features2 <- do.call(rbind, features2)
    colnames(features1) <- paste0("First_", functionNames)
    colnames(features2) <- paste0("Second_", functionNames)
    
    windows_hypoxia <- sapply(windows_list[[pid]], function(i) sum(i$timeDiff * i$hypoxicPeriods))
    windows_mean    <- sapply(windows_list[[pid]], function(i) mean(i$SPO2))
    windows_median  <- sapply(windows_list[[pid]], function(i) median(i$SPO2))
    windows_continuous <- sapply(windows_list[[pid]], function(i) any(i$timeDiff > 60))
    
    windows_hyp_1st <- sapply(windows_list[[pid]], function(i) {
      first_half <- i[1:(nrow(i) %/% 2), ]
      sum(first_half$SPO2 < 90) / (0.5 * nrow(i))
    })
    windows_hyp_2nd <- sapply(windows_list[[pid]], function(i) {
      second_half <- i[(nrow(i) %/% 2 + 1):nrow(i), ]
      sum(second_half$SPO2 < 90) / (0.5 * nrow(i))
    })
    df <- data.frame(
      time = windows_hypoxia,
      disjoint = windows_continuous,
      oxygen_mean = windows_mean,
      oxygen_median = windows_median,
      pre_hyp = windows_hyp_1st,
      post_hyp = windows_hyp_2nd,
      patientID = pid
    )
    df <- cbind(df, features1, features2,
                longest_hypoxic_seq = unlist(longest_sequences),
                longest_hypoxic_seq_88 = unlist(longest_sequences_88),
                longest_hypoxic_seq_89 = unlist(longest_sequences_89),
                longest_hypoxic_seq_91 = unlist(longest_sequences_91),
                longest_hypoxic_seq_92 = unlist(longest_sequences_92),
                maximumGap = unlist(maximumGap))
    df$meandiff <- c(0, diff(df$oxygen_mean))
    df$mediandiff <- c(0, diff(df$oxygen_median))
    summaries[[pid]] <- df
  }
  return(summaries)
}

summaries_10percOverlap <- summarize_patient_windows(windows4_10percOverlap, functions, functionNames)
summaries_30percOverlap <- summarize_patient_windows(windows4_30percOverlap, functions, functionNames)
summaries_50percOverlap <- summarize_patient_windows(windows4_50percOverlap, functions, functionNames)
summaries_70percOverlap <- summarize_patient_windows(windows4_70percOverlap, functions, functionNames)
summaries_90percOverlap <- summarize_patient_windows(windows4_90percOverlap, functions, functionNames)
 
summaries_10percOverlap <- dplyr::bind_rows(summaries_10percOverlap)
  summaries_cont_10 <- summaries_10percOverlap %>% dplyr::filter(disjoint==FALSE)
  summaries_stable_10 <- summaries_cont_10[summaries_cont_10$pre_hyp < 0.9 & summaries_cont_10$post_hyp < 0.9, ]
  summaries_deter_10 <- summaries_cont_10[summaries_cont_10$pre_hyp < 0.9 & summaries_cont_10$post_hyp >= 0.9, ]
summaries_30percOverlap <- dplyr::bind_rows(summaries_30percOverlap)
  summaries_cont_30 <- summaries_30percOverlap %>% dplyr::filter(disjoint==FALSE)
  summaries_stable_30 <- summaries_cont_30[summaries_cont_30$pre_hyp < 0.9 & summaries_cont_30$post_hyp < 0.9, ]
  summaries_deter_30 <- summaries_cont_30[summaries_cont_30$pre_hyp < 0.9 & summaries_cont_30$post_hyp >= 0.9, ]
summaries_50percOverlap <- dplyr::bind_rows(summaries_50percOverlap)
  summaries_cont_50 <- summaries_50percOverlap %>% dplyr::filter(disjoint==FALSE)
  summaries_stable_50 <- summaries_cont_50[summaries_cont_50$pre_hyp < 0.9 & summaries_cont_50$post_hyp < 0.9, ]
  summaries_deter_50 <- summaries_cont_50[summaries_cont_50$pre_hyp < 0.9 & summaries_cont_50$post_hyp >= 0.9, ]
summaries_70percOverlap <- dplyr::bind_rows(summaries_70percOverlap)
  summaries_cont_70 <- summaries_70percOverlap %>% dplyr::filter(disjoint==FALSE)
  summaries_stable_70 <- summaries_cont_70[summaries_cont_70$pre_hyp < 0.9 & summaries_cont_70$post_hyp < 0.9, ]
  summaries_deter_70 <- summaries_cont_70[summaries_cont_70$pre_hyp < 0.9 & summaries_cont_70$post_hyp >= 0.9, ]
summaries_90percOverlap <- dplyr::bind_rows(summaries_90percOverlap)
  summaries_cont_90 <- summaries_90percOverlap %>% dplyr::filter(disjoint==FALSE)
  summaries_stable_90 <- summaries_cont_90[summaries_cont_90$pre_hyp < 0.9 & summaries_cont_90$post_hyp < 0.9, ]
  summaries_deter_90 <- summaries_cont_90[summaries_cont_90$pre_hyp < 0.9 & summaries_cont_90$post_hyp >= 0.9, ]

rm(windows4_10percOverlap, summaries_10percOverlap, summaries_cont_10, summaries_stable_10, summaries_deter_10)
rm(windows4_30percOverlap, summaries_30percOverlap, summaries_cont_30, summaries_stable_30, summaries_deter_30)
rm(windows4_50percOverlap, summaries_50percOverlap, summaries_cont_50, summaries_stable_50, summaries_deter_50)
rm(windows4_70percOverlap, summaries_70percOverlap, summaries_cont_70, summaries_stable_70, summaries_deter_70)
rm(windows4_90percOverlap, summaries_90percOverlap, summaries_cont_90, summaries_stable_90, summaries_deter_90)


#90% overlap
#for (pid in seq(1,length(windows4)+1,1)){
  print(pid)
  if (pid != 461) {
    features1 <- list()
    features2 <- list()
    longest_sequences <- list()
    longest_sequences_88 <- list()
    longest_sequences_89 <- list()
    longest_sequences_91 <- list()
    longest_sequences_92 <- list()
    maximumGap <- list()
    for (i in seq_along(windows4[[pid]])) {
      current_window <- windows4[[pid]][[i]]
      first_half <- current_window[1:(nrow(current_window) %/% 2), ]
      second_half <- current_window[(nrow(current_window) %/% 2 + 1):nrow(current_window), ]
      features1[[i]] <- sapply(functions, function(f) f(first_half$SPO2))
      features2[[i]] <- sapply(functions, function(f) f(second_half$SPO2))
      # hypoxic_sequence <- rle(current_window$hypoxicPeriods) # Run-length encoding
      # longest_seq <- max(hypoxic_sequence$lengths[hypoxic_sequence$values == TRUE], na.rm = TRUE)
      # longest_sequences[[i]] <- ifelse(is.finite(longest_seq), longest_seq, 0) # Handle NA cases
      # Function to get longest consecutive values < threshold
      get_longest_seq <- function(x, threshold) {
        hypoxicThreshold <- x < threshold
        hypoxic_sequence <- rle(hypoxicThreshold)
        longest_seq <- max(hypoxic_sequence$lengths[hypoxic_sequence$values == TRUE], na.rm = TRUE)
        return(ifelse(is.finite(longest_seq), longest_seq, 0))
      }
      longest_sequences[[i]] <- get_longest_seq(current_window$SPO2, 90)
      longest_sequences_88[[i]] <- get_longest_seq(current_window$SPO2, 88)
      longest_sequences_89[[i]] <- get_longest_seq(current_window$SPO2, 89)
      longest_sequences_91[[i]] <- get_longest_seq(current_window$SPO2, 91)
      longest_sequences_92[[i]] <- get_longest_seq(current_window$SPO2, 92)
      maximumGap[[i]] <- max(current_window$timeDiff)
    }
    features1 <- do.call(rbind, features1)
    features2 <- do.call(rbind, features2)
    colnames(features1) <- paste0("First_", functionNames)
    colnames(features2) <- paste0("Second_", functionNames)
    windows_hypoxia = unlist(lapply(windows4[[pid]], function(i) return(sum(i$timeDiff*i$hypoxicPeriods))))
    windows_mean = unlist(lapply(windows4[[pid]], function(i) return(mean(i$SPO2))))
    windows_median = unlist(lapply(windows4[[pid]], function(i) return(median(i$SPO2))))
    windows_continuous = unlist(lapply(windows4[[pid]], function(i) return(any(i$timeDiff>60))))
    windows_hyp_1st = unlist(lapply(windows4[[pid]], function(i){
      first_half = i[1:(nrow(i) %/% 2), ]
      second_half = i[(nrow(i) %/% 2 + 1):nrow(i), ]
      # Check conditions for first and second halves
      return(sum(first_half$SPO2 < 90)/(0.5*nrow(i)))
    }))
    windows_hyp_2nd = unlist(lapply(windows4[[pid]], function(i){
      first_half = i[1:(nrow(i) %/% 2), ]
      second_half = i[(nrow(i) %/% 2 + 1):nrow(i), ]
      # Check conditions for first and second halves
      return(sum(second_half$SPO2 < 90)/(0.5*nrow(i)))
    }))
    df = data.frame(time = windows_hypoxia, disjoint = windows_continuous, 
                    oxygen_mean = windows_mean, oxygen_median=windows_median, 
                    pre_hyp= windows_hyp_1st, post_hyp=windows_hyp_2nd, 
                    patientID = pid)
    #df = cbind(df, features1, features2, longest_hypoxic_seq = unlist(longest_sequences), maximumGap = unlist(maximumGap))
    df = cbind(df, features1, features2,
               longest_hypoxic_seq = unlist(longest_sequences),
               longest_hypoxic_seq_88 = unlist(longest_sequences_88),
               longest_hypoxic_seq_89 = unlist(longest_sequences_89),
               longest_hypoxic_seq_91 = unlist(longest_sequences_91),
               longest_hypoxic_seq_92 = unlist(longest_sequences_92),
               maximumGap = unlist(maximumGap))
    df$meandiff <- c(0, diff(df$oxygen_mean))
    df$mediandiff <- c(0, diff(df$oxygen_median))
    summaries[[pid]]=df
  }
}

# summary_windows = dplyr::bind_rows(summaries)
# summary_windows10 = dplyr::bind_rows(summaries_10percOverlap)
#summary_windows_cont10 <- summary_windows10 %>% dplyr::filter(disjoint==FALSE)
#   stable_summary_windows_10 <- summary_windows_cont10[summary_windows_cont10$pre_hyp < 0.9 & summary_windows_cont10$post_hyp < 0.9, ]
#   detoriating_summary_windows_10 <- summary_windows_cont10[summary_windows_cont10$pre_hyp < 0.9 & summary_windows_cont10$post_hyp >= 0.9, ]
# summary_windows50 = dplyr::bind_rows(summaries_50percOverlap)
# summary_windows_cont50 <- summary_windows50 %>% dplyr::filter(disjoint==FALSE)
#   stable_summary_windows_50 <- summary_windows_cont50[summary_windows_cont50$pre_hyp < 0.9 & summary_windows_cont50$post_hyp < 0.9, ]
#   detoriating_summary_windows_50 <- summary_windows_cont50[summary_windows_cont50$pre_hyp < 0.9 & summary_windows_cont50$post_hyp >= 0.9, ]
# summaries_50perc90minus = dplyr::bind_rows(summaries_50perc90minus)
# summary_windows_cont50_90minus <- summaries_50perc90minus %>% dplyr::filter(disjoint==FALSE)
#   stable_summary_windows50_90minus <- summary_windows_cont50_90minus[summary_windows_cont50_90minus$pre_hyp < 0.9 & summary_windows_cont50_90minus$post_hyp < 0.9, ]
#   detoriating_summary_windows50_90minus <- summary_windows_cont50_90minus[summary_windows_cont50_90minus$pre_hyp < 0.9 & summary_windows_cont50_90minus$post_hyp >= 0.9, ]
#   
# 
# length(which(summary_windows$disjoint == TRUE))
# length(unique(summary_windows$patientID[summary_windows$disjoint == TRUE]))
# #summary_windows_cont = summary_windows %>% dplyr::filter(disjoint==FALSE, pre_hyp<0.9)
# summary_windows_cont <- summary_windows %>% dplyr::filter(disjoint==FALSE)
# summary_windows_disjoint <-  summary_windows %>% dplyr::filter(disjoint==TRUE)

#plotting sample SPO2 traces ----
  #stable
  plot(windows4[[1]][[1]]$SPO2, type = "l", col = "black",  xlab = "Time", ylab = "SPO2", ylim = c(80, 100))
  abline(h = 90, col = "red", lty = 2)
  abline(v = 120, col = "blue", lwd = 1)
  #deteriorating
  plot(windows4[[2]][[317]]$SPO2, type = "l", col = "black",  xlab = "Time", ylab = "SPO2", ylim = c(80, 100))
  abline(h = 90, col = "red", lty = 2)
  abline(v = 120, col = "blue", lwd = 1)
  #extended
  plot(windows4[[2]][[289]]$SPO2, type = "l", col = "black",  xlab = "Time", ylab = "SPO2", ylim = c(80, 100))
  abline(h = 90, col = "red", lty = 2)
  abline(v = 120, col = "blue", lwd = 1)
  #improving
  plot(windows4[[1]][[124]]$SPO2, type = "l", col = "black",  xlab = "Time", ylab = "SPO2", ylim = c(80, 100))
  abline(h = 90, col = "red", lty = 2)
  abline(v = 120, col = "blue", lwd = 1)
# feature correlation + plotting ----
  selected <- c("longest_hypoxic_seq", "First_Min", "First_Mean", "First_Energy", 
                "First_Median", "First_Max")
  ggpairs(summary_windows_cont[,selected]) + labs(title = "Most important features in pre-hypoxia windows")
  ggpairs(summary_windows_cont[,c(8:25, 44)]) + labs(title = "All. pre-hypoxia window features")

  # plot 2
  functionNames <- c("Mean", "Std", "Median", "Max", "Min", 
                     "MaxMinusMin", "MeanAbsoluteDiff", "MedianAbsoluteDiff",
                     "IQR", "ValsAboveMean", "ValsAboveGeoMean", "NumsPeaks",
                     "Skew", "Kurtosis", "Energy", "Argmin", "Argmax", "ArgDiff")
  
  plot_data <- do.call(rbind, lapply(functionNames, function(feature) {
    data.frame(Feature = feature,
      First = summary_windows_cont[[paste0("First_", feature)]],
      Second = summary_windows_cont[[paste0("Second_", feature)]])}))
  # 1. Calculate R-squared for each feature and store it in a separate variable
  r_squared_values <- plot_data %>%
    group_by(Feature) %>%
    do({
      model <- lm(Second ~ First, data = .)  
      r_squared <- summary(model)$r.squared  
      data.frame(R_squared = r_squared)
    })
  # 2. Reorder 'Feature' based on R-squared values (to sort facets by R-squared)
  plot_data$Feature <- factor(plot_data$Feature, levels = r_squared_values$Feature[order(r_squared_values$R_squared, decreasing = TRUE)])
  ggplot(plot_data, aes(x = First, y = Second)) +
    geom_point(alpha = 0.7, color = "lightcoral",  size = 0.2) +
    geom_smooth(method = "lm", se = FALSE, color = "lightseagreen") +
    facet_wrap(~ Feature, scales = "free") + # Facet by feature with independent scales
    theme_minimal() + labs(x = "Pre-window feature value", y = "Target window feature value",
      caption = "Linear regression and R-squared values")
  for (feat in levels(plot_data$Feature)) {
    df <- plot_data %>% filter(Feature == feat)
    p <- ggplot(df, aes(x = First, y = Second)) +
      geom_point(alpha = 0.7, color = "lightcoral", size = 0.2) +
      geom_smooth(method = "lm", se = FALSE, color = "lightseagreen") +
      theme_minimal() +
      labs(title = feat, x = "Pre-window feature value", y = "Target window feature value",
           caption = paste0("R² = ", round(r_squared_values$R_squared[r_squared_values$Feature == feat], 3)))
    print(p)  
  }
# different windows ----
  stable_summary_windows <- summary_windows_cont[summary_windows_cont$pre_hyp < 0.9 & summary_windows_cont$post_hyp < 0.9, ]
  detoriating_summary_windows <- summary_windows_cont[summary_windows_cont$pre_hyp < 0.9 & summary_windows_cont$post_hyp >= 0.9, ]
  extened_summary_windows <- summary_windows_cont[summary_windows_cont$pre_hyp >= 0.9 & summary_windows_cont$post_hyp >= 0.9, ]
  improving_summary_windows <- summary_windows_cont[summary_windows_cont$pre_hyp >= 0.9 & summary_windows_cont$post_hyp < 0.9, ]
  stable_summary_windows$Group <- "Stable"
  detoriating_summary_windows$Group <- "Deteriorating"
  combined_data <- bind_rows(stable_summary_windows, detoriating_summary_windows)
  dim(stable_summary_windows)
  dim(extened_summary_windows)
  dim(detoriating_summary_windows)
  dim(improving_summary_windows)
  length(unique(stable_summary_windows$patientID))
  length(unique(extened_summary_windows$patientID))
  length(unique(detoriating_summary_windows$patientID))
  length(unique(improving_summary_windows$patientID))
  selected_features <- c("First_Mean", "First_Std", "First_Median", "First_Max", "First_Min", "First_MaxMinusMin",
                         "First_MeanAbsoluteDiff", "First_MedianAbsoluteDiff", "First_IQR", "First_ValsAboveMean",
                         "First_ValsAboveGeoMean", "First_NumsPeaks", "First_Skew", "First_Kurtosis", "First_Energy",
                         "First_Argmin", "First_Argmax", "First_ArgDiff")
  combined_data <- combined_data %>% dplyr::select(Group, all_of(selected_features)) %>%
    pivot_longer(cols = -Group, names_to = "Feature", values_to = "Value")
  ggplot(combined_data, aes(x = Group, y = Value, fill = Group)) +
    geom_boxplot(outlier.alpha = 0.5) + # Boxplot without outliers for cleaner appearance
    facet_wrap(~ Feature, scales = "free_y", ncol = 4) + # Separate plots for each feature
    theme_minimal() +
    labs(title = "Boxplots for Each Feature by Group",
         x = "Group",
         y = "Value") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(size = 8)) # Adjust facet text size

  df_cw = summary_windows %>% dplyr::filter(disjoint==FALSE, pre_hyp<0.9)
  sample <- sample(c(TRUE, FALSE), nrow(df_cw), replace=TRUE, prob=c(0.7,0.3))
  train  <- df_cw[sample, ]
  test   <- df_cw[!sample, ]
  #plotting features of stable vs summary vs
  stable_summary_windows$group <- "Stable"
  detoriating_summary_windows$group <- "Deteriorating"
  extened_summary_windows$group <- "Extended"
  combined_data_plot <- bind_rows(stable_summary_windows, detoriating_summary_windows)
  combined_data_plot <- combined_data_plot %>%
    dplyr::select(starts_with("First_"), group) %>%
    pivot_longer(cols = -group, names_to = "Feature", values_to = "Value") %>%
    mutate(Feature = gsub("^First_", "", Feature),  # Remove "First_" from Feature names
           group = fct_relevel(group, "Stable", "Deteriorating"))  # Reorder factor levels
  # Create the boxplot
  combined_data_plot <- combined_data_plot  %>% mutate(Feature = as.factor(Feature)) 
  axis_limits <- combined_data_plot %>% group_by(Feature) %>%
    summarise(ymin = quantile(Value, 0.05, na.rm = TRUE), ymax = quantile(Value, 0.95, na.rm = TRUE))
  combined_data_plot <- combined_data_plot %>% left_join(axis_limits, by = "Feature")
  ggplot(combined_data_plot, aes(x = group, y = Value, fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~ Feature, scales = "free", ncol = 6, nrow = 3) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),strip.background = element_blank()) +
    labs(title = "Comparison of Features Across Groups", x = "Group", y = "Features", fill = "Group") +
    geom_signif(comparisons = list(c("Stable", "Deteriorating"),  c("Deteriorating", "Extended"), c("Stable", "Extended")), 
      map_signif_level = TRUE, textsize = 3, step_increase = 0.05) 
  # window feature comparisons table
  comparisons_table <- combined_data_plot %>%
    group_by(Feature) %>%
    t_test(Value ~ group) %>%
    adjust_pvalue(method = "bonferroni") 
  # creating separate plots
   for (feat in unique(combined_data_plot$Feature)) {
     df <- combined_data_plot %>% filter(Feature == feat)
     p <- ggplot(df, aes(x = group, y = Value, fill = group)) +
       geom_boxplot() +
       theme_minimal() +
       labs(title = feat, x = "Group", y = "Value") +
       theme(legend.position = "none") +
       geom_signif(comparisons = list(c("Stable", "Deteriorating")), 
                   map_signif_level = TRUE)  # Add significance comparisons
     print(p)  # Print each plot separately
   }
   
