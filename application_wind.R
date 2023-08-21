## ============================================================================================

# read the data
wind <- read.csv2("application_wind.csv")

# separating training and testing
training <- wind$average_wind_speed[1:62]
test <- wind$average_wind_speed[65:74] #obs 63 and 64 are NA

# descritive analysis
summary(training)

# some graphics  
acf(training,ylab="ACF",xlab="Lag",main="")
monthplot(training, ylab = "Average wind speed (m/s)", xlab = "Months")
plot.ts(training,ylab="Average wind speed (m/s)",xlab="Time (months)")

# covariates
x <- ts(training, start= c(2009, 12), frequency = 12)
saz <- decompose(x)
c <- saz$seasonal
xreg <- as.vector(c)
c_hat <- xreg[3:14]

#----------------------------------------------

source("charma_fit.r")
tau <- 0.5 # median

#----------------------------------------------
# fit the model CHARMA(3,2) with covariates
fit_training <- try(charma.fit(y=training, ar=c(1,2,3), ma=c(1,2), 
                               X=c, X_hat = c_hat, h1 = 12, 
                               graf=1), T)


# accuracy measures of the CHARMA model (in-sample)
dif_charma <- (training[4:62]-fit_training$fitted[4:62])
(mse_charma <- (sum(dif_charma^2))/length(dif_charma))
(mape_charma <- sum(abs(dif_charma)/abs(training[4:62]))/length(dif_charma))

forecast::accuracy(fit_training$fitted,training[1:62])

# evaluating the prediction of the CHARMA model (out-of-sample)
predict_charma <- fit_training$forecast[3:12]

dif_charmap <- (test-predict_charma)
(mse_charmap <- (sum(dif_charmap^2))/length(dif_charmap))
(mape_charmap <- sum(abs(dif_charmap)/abs(test))/length(dif_charmap))


# fit the model SARMA(2,0,3)(2,0,0) 
sarma_training <- arima(x = training,order = c(2,0,3),
                        seasonal = list(order = c(2, 0, 0), period = 12))

# accuracy measures of the SARMA model (in-sample)
difsarma <- (sarma_training$residuals)
(mse_sarma <- (sum(difsarma^2))/length(difsarma))
(mape_sarma <- sum(abs(difsarma)/abs(training))/length(difsarma))

# evaluating the prediction of the SARMA model (out-of-sample)
predict_sarma <- predict(sarma_training, n.ahead=12)

dif_out_sarma <- (test-predict_sarma$pred[3:12])
(mse_sarmap <- (sum(dif_out_sarma^2))/length(dif_out_sarma))
(mape_sarmap <- sum(abs(dif_out_sarma)/abs(test))/length(dif_out_sarma))

prev_sarma <- training-sarma_training$residuals


M <- cbind(c(mse_charma, mape_charma),
           c(mse_sarma ,mape_sarma ))
colnames(M) <- c("CHARMA", "SARMA")
rownames(M) <- c("EQM", "MAPE(%)")
M

Mprev <- cbind(c(mse_charmap, mape_charmap*100),
               c(mse_sarmap,mape_sarmap*100))
colnames(Mprev)<-c("CHARMA", "SARMA")
rownames(Mprev)<-c("EQM", "MAPE(%)")
Mprev

#-------------------------
# graphic
plot(wind$average_wind_speed, col="black", type="l", 
     ylim=c(2, (max(c(training, test))+0.75)) , axes =T, main="", xlab="Time (months)", ylab="Average wind speed (m/s)")
lines(c(fit_training$fitted, fit_training$forecast), lty = 2, lwd = 1, col="red")
lines(c(prev_sarma, predict_sarma$pred), lty = 3, lwd = 1, col="blue4")
legend("top",legend=c( "Observations", "CHARMA", "SARMA"),
       pt.bg="white", lty=c(1,2,3), col=c("black","red","blue4"), bty="n" )
fim <- (end(wind$average_wind_speed)[1]+end(wind$average_wind_speed)[2]/12)-12
abline(v=fim,lty=2)
