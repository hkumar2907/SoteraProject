library(ggplot2)
library(gridExtra)
library(dplyr)
library(readr)
#DATA ----
set.seed(123)
x <- 1:100
y <- sin(x/10) + rnorm(100, mean = 0, sd = 0.1)
yNoisy <- y + rnorm(100, mean = 0, sd = 0.5)
yAmp <- y*5
yPerm <- c(y[1:24],NA, y[51:74],NA, y[26:49], NA, y[76:100])
y2 <- 0.5*sin(x/5 + 1) + rnorm(100, mean = 0, sd = 0.1) + 0.3
y3 <- 0.7*sin(x/8 + 0.5) + rnorm(100, mean = 0, sd = 0.1)
figDF <- as.data.frame(cbind(x,y,yNoisy, yAmp, yPerm, y2, y3))
#NOISE
fig1 <- ggplot(figDF, aes(x, y)) +geom_line(color = 'blue') + ggtitle("Original Data") +theme_minimal() + labs(x = "Time(seconds)", y = "x") 
fig2 <- ggplot(figDF, aes(x, yNoisy)) + geom_line(color = 'red') + ggtitle("Data with Noise") +theme_minimal() + labs(x = "Time(seconds)", y = "x") 
grid.arrange(fig1, fig2, ncol = 2)
#INCREASING AMPLITUDE
fig1 <- ggplot(figDF, aes(x, y)) +geom_line(color = 'blue') + ggtitle("Original Data") +theme_minimal() + labs(x = "Time(seconds)", y = "x")  + scale_y_continuous(limits = c(-5, 5))
fig2 <- ggplot(figDF, aes(x, yAmp)) + geom_line(color = 'red') + ggtitle("Amplitude Scaling") +theme_minimal() + labs(x = "Time(seconds)", y = "x")  + scale_y_continuous(limits = c(-5, 5))
grid.arrange(fig1, fig2, ncol = 2)
#PERMUTATION
fig1 <- ggplot(figDF, aes(x, y)) +geom_line(color = 'blue') + ggtitle("Original Data") +theme_minimal() + labs(x = "Time(seconds)", y = "x") 
fig2 <- ggplot(figDF) + geom_line(size = 0.5, color = "red", aes(x, yPerm)) + ggtitle("Permutation") +theme_minimal() + labs(x = "Time(seconds)", y = "x")  
grid.arrange(fig1, fig2, ncol = 2)
#SWAPPING
fig1 <- ggplot(figDF) + geom_line(aes(x,y),color = 'blue', linewidth = 1, linetype = "dotdash") + geom_line(aes(x,y2),color = 'red', linewidth = 0.5, linetype = "dashed") + geom_line(aes(x,y3),color = 'black', linewidth = 0.75)+ ggtitle("Original Data") +theme_minimal() + labs(x = "Time(seconds)", y = "x") 
fig2 <- ggplot(figDF) + geom_line(aes(x,y2),color = 'blue', linewidth = 1, linetype = "dotdash") + geom_line(aes(x,y),color = 'red', linewidth = 0.5, linetype = "dashed") + geom_line(aes(x,y3),color = 'black', linewidth = 0.75)+ ggtitle("Swapping") +theme_minimal() + labs(x = "Time(seconds)", y = "x") + theme(legend.position = "bottom")
grid.arrange(fig1, fig2, ncol = 2)
#SPHERICAL COORDINATES
rectToSph <- function(xyz) {
  xSqySq <- xyz[,1]^2 + xyz[,2]^2
  y <- sqrt(xSqySq + xyz[,3]^2)
  y2 <- atan2(sqrt(xSqySq), xyz[,3])
  y3 <- atan2(xyz[,2], xyz[,1])
  new <- cbind(y, y2, y3)
  return(new)
}
spherical = rectToSph(figDF[,c(2,6,7)])
fig1 <- ggplot(figDF) + geom_line(aes(x,y),color = 'blue', linewidth = 1) + geom_line(aes(x,y2),color = 'red', linewidth = 0.5) + geom_line(aes(x,y3),color = 'black')+ ggtitle("Original Data") +theme_minimal() + labs(x = "Time(seconds)", y = "x") 
fig2 <- ggplot(spherical) + geom_line(aes(x,y2),color = 'blue', linewidth = 1) + geom_line(aes(x,y),color = 'red', linewidth = 0.5) + geom_line(aes(x,y3),color = 'black')+ ggtitle("Spherical Data") +theme_minimal() + labs(x = "Time(seconds)", y = "x") + theme(legend.position = "bottom")
grid.arrange(fig1, fig2, ncol = 2)

#TDA ----
radiomicsdata <- read_csv("radiomicsdata.csv")
#noise
radiomicsdata$xNoisy <- radiomicsdata$x_axis + rnorm(49, mean = 0, sd = 1)
fig1 <- ggplot(radiomicsdata, aes(...1, x_axis)) +geom_line(color = 'blue') + ggtitle("Original Data") +theme_minimal() + labs(x = "Time(seconds)", y = "x") 
fig2 <- ggplot(radiomicsdata, aes(...1, xNoisy)) + geom_line(color = 'red') + ggtitle("Data with Noise") +theme_minimal() + labs(x = "Time(seconds)", y = "x") 
grid.arrange(fig1, fig2, ncol = 2)
#increased amplityde
radiomicsdata$xAmp <- radiomicsdata$x_axis * 3
fig1 <- ggplot(radiomicsdata, aes(...1, x_axis)) +geom_line(color = 'blue') + ggtitle("Original Data") +theme_minimal() + labs(x = "Time(seconds)", y = "x") + scale_y_continuous(limits = c(-9, 9))
fig2 <- ggplot(radiomicsdata, aes(...1, xAmp)) + geom_line(color = 'red') + ggtitle("Amplitude Scaling") +theme_minimal() + labs(x = "Time(seconds)", y = "x") + scale_y_continuous(limits = c(-9, 9))
grid.arrange(fig1, fig2, ncol = 2)
#permutation
radiomicsdata$xPerm <- c(radiomicsdata$x_axis[1:12],NA, radiomicsdata$x_axis[26:37],NA, radiomicsdata$x_axis[14:24], NA, radiomicsdata$x_axis[39:49])
fig1 <- ggplot(radiomicsdata, aes(...1, x_axis)) +geom_line(color = 'blue') + ggtitle("Original Data") +theme_minimal() + labs(x = "Time(seconds)", y = "x") 
fig2 <- ggplot(radiomicsdata, aes(...1, xPerm)) + geom_line(color = 'red') + ggtitle("Permutation") +theme_minimal() + labs(x = "Time(seconds)", y = "x") 
grid.arrange(fig1, fig2, ncol = 2)
#swapping
fig1 <- ggplot(radiomicsdata) + geom_line(aes(...1,x_axis),color = 'blue', linewidth = 1, linetype = "dotdash") + geom_line(aes(...1,y_axis),color = 'red', linewidth = 0.5, linetype = "dashed") + geom_line(aes(...1,z_axis),color = 'black', linewidth = 0.75)+ ggtitle("Original Data") +theme_minimal() + labs(x = "Time(seconds)", y = "x") 
fig2 <- ggplot(radiomicsdata) + geom_line(aes(...1,y_axis),color = 'blue', linewidth = 1, linetype = "dotdash") + geom_line(aes(...1,x_axis),color = 'red', linewidth = 0.5, linetype = "dashed") + geom_line(aes(...1,z_axis),color = 'black', linewidth = 0.75)+ ggtitle("Swapping") +theme_minimal() + labs(x = "Time(seconds)", y = "x") + theme(legend.position = "bottom") +
  grid.arrange(fig1, fig2, ncol = 2)
#spherical
spherical = as.data.frame(rectToSph(cbind(as.numeric(radiomicsdata$x_axis), as.numeric(radiomicsdata$y_axis), as.numeric(radiomicsdata$z_axis))))
spherical$t <- radiomicsdata$...1
ggplot(radiomicsdata) + geom_line(aes(...1,x_axis,color = "x_axis", linetype = "x_axis"), linewidth = 1) + geom_line(aes(...1,y_axis, color = "y_axis", linetype = "y_axis"), linewidth = 0.5) + geom_line(aes(...1,z_axis,  color = "z_axis", linetype = "z_axis"), linewidth = 0.75)+ ggtitle("Original Data") +theme_minimal() + labs(x = "Time(seconds)", y = "x") 
+  scale_color_manual(values = c("x_axis" = "red", "y_axis" = "blue", "z_axis" = "black")) + scale_linetype_manual(values = c("x_axis" = "dashed", "y_axis" = "dotdash", "z_axis" = "solid")) + theme(legend.position = "bottom") +
  guides(color = guide_legend("Variables"), linetype = guide_legend("Variables"))

fig2 <- ggplot(spherical) +  geom_line(aes(t,y),color = 'blue', linewidth = 1, linetype = "dotdash") + geom_line(aes(t,y2),color = 'red', linewidth = 0.5, linetype = "dashed") + geom_line(aes(t,y3),color = 'black', linewidth = 0.75)+ ggtitle("Spherical Data") +theme_minimal() + labs(x = "Time(seconds)", y = "x") 
grid.arrange(fig1, fig2, ncol = 2)
