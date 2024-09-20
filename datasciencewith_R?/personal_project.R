set.seed(1001)
x1=1:100+rnorm(100,mean=0,sd=15)
y1=1:100


plot(x1, y1, main = "Scatter Plot of x1 vs y1", xlab = "x1 Values", ylab = "y1 Values", col = "blue", pch = 19)

mtext(side = 3, text = "hi there")  # Text at the top margin
mtext(side = 2, text = "hi there")  # Text at the left margin


correlation <- cor(x1, y1)
mtext(side = 3, text = paste("Correlation: ", round(correlation, 2)))


plot(x1, y1, main = "Scatter Plot of x1 vs y1", xlab = "x1 Values", ylab = "y1 Values", col = "red", pch = 18)


hist(x1, col = "lightblue", xlab = "x1 Values", ylab = "Frequency", main = "Histogram of x1")


##Box plot
boxplot(y1, main = "Boxplot of y1")
boxplot(x1, y1, names=c("x1", "y1"), horizontal=TRUE)

##multiple plot
par(mfrow=c(2, 1))
boxplot(y1, main="Boxplot of y1")
hist(x1, col="lightgreen", xlab="x1 Values", ylab="Frequency", main="Histogram of x1")
