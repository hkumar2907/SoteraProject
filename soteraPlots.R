library(readr)
library(ggplot2)
visi_data <- read_csv("visi_data.csv")
plot(visi_data$CNIBP_MAP)
plot(visi_data$SPO2)
#assuming each point is 15 seconds
#Desaturation episodes: within 30 minutes, atleast 90% of the measurements are <90% oxygen saturation 
#Hypotensive episodes: MAP < 70 mm Hg for atleast 15 minutes
visi_data <- as.data.frame(visi_data)
visi_data$CNIBP_MAP <- as.numeric(visi_data$CNIBP_MAP)
visi_data$SPO2 <- as.numeric(visi_data$SPO2)
visi_data$DE <- 0
visi_data$HE <- 0
visi_data$indx <- rownames(visi_data)
DE_len <- 120
HE_len <- 60
DETracker <- visi_data$SPO2[1:DE_len]
HETracker <- visi_data$CNIBP_MAP[1:HE_len]
x <- 0
for (i in DE_len:nrow(visi_data)) {
  if ((length(na.omit(DETracker)) != 0) && !is.na(visi_data$SPO2[i])) {
    if (sort(DETracker, decreasing = TRUE)[length(na.omit(DETracker))*0.9] < 90) { #starting at the start
      visi_data$DE[i-DE_len] <- 1
    }
  }
  DETracker <- c(DETracker[-1], visi_data$SPO2[i])
  x <- x+1
}
for (i in (HE_len+1):nrow(visi_data)) {
  if (all(na.omit(HETracker) < 70) && !is.na(visi_data$CNIBP_MAP[i])) {
    #for this one want the whole window to be green
    visi_data$HE[(i-HE_len):i] <- 1
  }
  HETracker <- c(HETracker[-1], visi_data$CNIBP_MAP[i])
}
visi_data$indx <- as.numeric(visi_data$indx) * (15/3600) #changing indices to seconds then to hours
ggplot(visi_data, aes(x = as.numeric(indx), y = SPO2)) + geom_line(col = ifelse(visi_data$DE == 1,'red','black')) + geom_hline(yintercept=90,linetype=2) + ggtitle("Sample SPO2 Values", ) + xlab("Time(hours)") + ylab("SPO2")
ggplot(visi_data, aes(x = as.numeric(indx), y = CNIBP_MAP)) + geom_line(col = ifelse(visi_data$HE == 1,'green','black')) + geom_hline(yintercept=70,linetype=2) + ggtitle("Sample MAP Values") + xlab("Time(hours)") + ylab("MAP")
