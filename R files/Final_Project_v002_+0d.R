#install.packages("forecast") #-- do this only once 
#Check the book: https://www.otexts.org/fpp2 and the blog: http://robjhyndman.com/hyndsight 
library(forecast)

covid_data<-read.csv(file.choose(), header=TRUE, sep=",")

FRA <- ts(covid_data$deaths_FRA,start=0, frequency=1)
cFRA <- ts(covid_data$Cdeaths_FRA,start=0, frequency=1)
cFRA2 <- ts(covid_data$Cdeaths_FRA[45:75],start=2020+52/365, frequency=365) #t=0 @ day75. +14d
cFRA3 <- ts(covid_data$Cdeaths_FRA[45:146],start=2020+52/365, frequency=365)
#FRA <- ts(covid_data$deaths_FRA,start=2020, frequency=365)
#cFRA <- ts(covid_data$Cdeaths_FRA,start=2020, frequency=365)

plot(cFRA)
cFRA <- FRA



fit <- decompose(cFRA, type="multiplicative") #decompose using "classical" method, multiplicative form
plot(fit)

fit <- decompose(cFRA, type="additive") #decompose using "classical" method, additive form
plot(fit)

#fit <- stl(cFRA, t.window=12, s.window="periodic") #decompose using STL (Season and trend using Loess)

cFRA2<-cFRA2+1

cFRA2_AAN <- ets(cFRA2, model="AAN", damped=FALSE)
cFRA2_AAZ <- ets(cFRA2, model="AAZ", damped=FALSE)
cFRA2_MMN <- ets(cFRA2, model="MMN", damped=FALSE)
cFRA2_MMZ <- ets(cFRA2, model="MMZ", damped=FALSE)

cFRA2_AAND <- ets(cFRA2, model="AAN", damped=TRUE)
cFRA2_AAZD <- ets(cFRA2, model="AAZ", damped=TRUE)
cFRA2_MMND <- ets(cFRA2, model="MMN", damped=TRUE)
cFRA2_MMZD <- ets(cFRA2, model="MMZ", damped=TRUE)


# Create their prediction "cones" for 360 months (30 years) into the future with quintile confidence intervals
cFRA2_AAN_pred <- forecast(cFRA2_AAN, h=30, level=c(0.8, 0.95))
cFRA2_AAZ_pred <- forecast(cFRA2_AAZ, h=30, level=c(0.8, 0.95))
cFRA2_MMN_pred <- forecast(cFRA2_MMN, h=30, level=c(0.8, 0.95))
cFRA2_MMZ_pred <- forecast(cFRA2_MMZ, h=30, level=c(0.8, 0.95))

cFRA2_AAND_pred <- forecast(cFRA2_AAND, h=30, level=c(0.8, 0.95))
cFRA2_AAZD_pred <- forecast(cFRA2_AAZD, h=30, level=c(0.8, 0.95))
cFRA2_MMND_pred <- forecast(cFRA2_MMND, h=30, level=c(0.8, 0.95))
cFRA2_MMZD_pred <- forecast(cFRA2_MMZD, h=30, level=c(0.8, 0.95))

par(mfrow=c(1,4)) # This command sets the plot window to show 1 row of 4 plots
plot(cFRA2_AAN_pred, xlab="Year", ylab="Acummulated deaths", ylim=c(0,100000))
plot(cFRA2_MMN_pred, xlab="Year", ylab="Acummulated deaths", ylim=c(0,100000))
plot(cFRA2_AAZ_pred, xlab="Year", ylab="Acummulated deaths", ylim=c(0,100000))
plot(cFRA2_MMZ_pred, xlab="Year", ylab="Acummulated deaths", ylim=c(0,100000))

par(mfrow=c(1,4)) # This command sets the plot window to show 1 row of 4 plots
plot(cFRA2_AAND_pred, xlab="Year", ylab="Acummulated deaths", ylim=c(0,100000))
plot(cFRA2_MMND_pred, xlab="Year", ylab="Acummulated deaths", ylim=c(0,100000))
plot(cFRA2_AAZD_pred, xlab="Year", ylab="Acummulated deaths", ylim=c(0,100000))
plot(cFRA2_MMZD_pred, xlab="Year", ylab="Acummulated deaths", ylim=c(0,100000))

par(mfrow=c(1,4)) # This command sets the plot window to show 1 row of 4 plots
plot(cFRA2_AAN_pred, xlab="Year", ylab="Acummulated deaths", ylim=c(0,100000))
lines(cFRA3, col="red")
plot(cFRA2_MMN_pred, xlab="Year", ylab="Acummulated deaths", ylim=c(0,100000))
lines(cFRA3, col="red")
plot(cFRA2_AAND_pred, xlab="Year", ylab="Acummulated deaths", ylim=c(0,100000))
lines(cFRA3, col="red")
plot(cFRA2_MMND_pred, xlab="Year", ylab="Acummulated deaths", ylim=c(0,100000))
lines(cFRA3, col="red")

cFRA2_tbats <- tbats(cFRA2)
cFRA2_tbats_pred <-forecast(cFRA2_tbats, h=30, level=c(0.8, 0.95))
par(mfrow=c(1,1))
plot(cFRA2_tbats_pred, xlab="Year", ylab="Predicted Electric Rate")




###
### Comparing models -- Time series Cross Validation (Rolling Horizon Holdout)
###

f_AANF  <- function(y, h) forecast(ets(y, model="AAN"), h = h)
errors_AANF <- tsCV(cFRA2, f_AANF, h=1, window=27)

f_MMNF<- function(y, h) forecast(ets(y, model="MMN"), h = h)
errors_MMNF <- tsCV(cFRA2, f_MMNF, h=1, window=27)

f_AAND  <- function(y, h) forecast(ets(y, model="AAN", damped=TRUE), h = h)
errors_AAND <- tsCV(cFRA2, f_AAND, h=1, window=27)

f_MMND  <- function(y, h) forecast(ets(y, model="MMN", damped=TRUE), h = h)
errors_MMND <- tsCV(cFRA2, f_MMND, h=1, window=27)

f_TBATS  <- function(y, h) forecast(tbats(y), h = h)
errors_TBATS <- tsCV(cFRA2, f_TBATS, h=1, window=27)

par(mfrow=c(1,1)) 
plot(errors_AANF, ylab='tsCV errors')
abline(0,0)
lines(errors_MMNF, col="red")
lines(errors_AAND, col="green")
lines(errors_MMND, col="blue")
lines(errors_TBATS, col="grey")
legend("left", legend=c("CV_error_AANF", "CV_error_MMNF","CV_error_AAND","CV_error_MMND", "CV_error_TBATS"), col=c("black", "red", "green", "blue", "grey"), lty=1:4)

mean(abs(errors_AANF/cFRA2), na.rm=TRUE)*100
mean(abs(errors_MMNF/cFRA2), na.rm=TRUE)*100
mean(abs(errors_AAND/cFRA2), na.rm=TRUE)*100
mean(abs(errors_MMND/cFRA2), na.rm=TRUE)*100
mean(abs(errors_TBATS/cFRA2), na.rm=TRUE)*100

# Print the mean and confidence intervals for the MMZ model
cFRA2_MMN_pred

# Export the results out
write.csv(cFRA2_tbats_pred, file = "FRA no lockdown+0d (BATS).csv") # export the selected model's predictions into a CSV file

par(mfrow=c(1,1)) # This command sets the plot window to show 1 row of 4 plots
plot(cFRA2_tbats_pred, xlab="Year", ylab="Acummulated deaths", ylim=c(0,80000))
lines(cFRA3, col="red")
abline(v=2020+75/365, col="green")
legend("left", legend=c("Best fitted model","Real death toll", "Start of Lockdown"), col=c("blue", "red","green"), lty=1)

library(forecast)
library(fpp)

cFRA21 <-cFRA2+10
View(cFRA21)
plot(cFRA21)
plot(log(cFRA21))
plot(log(log(cFRA21)))


par(mfrow=c(1,3))
plot(cFRA21, xlab="Year",
     ylab="A10 sales")
plot(log(cFRA21), xlab="Year",
     ylab="log A10 sales")
plot(diff(log(cFRA21),12), xlab="Year",
     ylab="Annual change in monthly log A10 sales")

# decomposition 
fit <- stl(log(cFRA21), t.window=12, s.window="periodic", robust=TRUE)
plot(fit)

#a10 dataset from fpp - sales of antidiabetic drug in Australia
par(mfrow=c(1,2))
Acf(diff(log(cFRA21),12)) # auto-correlation function
Pacf(diff(log(cFRA21),12)) # partial auto-correlation function


# fit ARIMA models

# non-seasonal first
fit <- auto.arima(cFRA21,seasonal=FALSE)
fit

#check residuals for autocorrelation
par(mfrow=c(1,1))
Acf(residuals(fit))

# and now seasonal
fit <- auto.arima(cFRA21,seasonal=TRUE)
fit

plot(forecast(fit,30)) #make a prediction for 5 years (60 months)

###
