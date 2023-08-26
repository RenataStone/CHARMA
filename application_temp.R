## ============================================================================================

# read the data
temp <- read.csv2("temperature.csv")

# separating training and testing
training <- temp$temperature[1:62]
test <- temp$temperature[64:71] #obs 63 is NA

#descritive analysis
summary(training)

# some graphics  
acf(training,ylab="ACF",xlab="Lag",main="")
monthplot(training, ylab = "Maximum temperature (°C)", xlab = "Months")
plot.ts(training,ylab="Maximum temperature (°C)",xlab="Time (months)")

# covariates
x <- ts(training, start= c(2010, 2), frequency = 12)
saz <- decompose(x)
c <- saz$seasonal
xreg <- as.vector(c)
c_hat <- xreg[51:60]

#----------------------------------------------

source("charma_fit.r")
tau <- 0.5 # median

#----------------------------------------------
# fit the model CHARMA(3,0)* with covariates
fit_training <- try(charma.fit(y=training, ar=c(1,3), ma=c(NA), 
                               X=c, X_hat = c_hat, h1 = 9, 
                               graf=1), T)


# accuracy measures of the CHARMA model (in-sample)
dif_charma <- (training[4:62]-fit_training$fitted[4:62])
(mse_charma <- (sum(dif_charma^2))/length(dif_charma))
(mape_charma <- sum(abs(dif_charma)/abs(training[4:62]))/length(dif_charma))

forecast::accuracy(fit_training$fitted,training[1:62])

# evaluating the prediction of the CHARMA model (out-of-sample)
predict_charma <- fit_training$forecast[2:9]

dif_charmap <- (test-predict_charma)
(mse_charmap <- (sum(dif_charmap^2))/length(dif_charmap))
(mape_charmap <- sum(abs(dif_charmap)/abs(test))/length(dif_charmap))


# fit the model ARMA(2,0,1) with covariates
arma_training <- forecast::auto.arima(training)


# accuracy measures of the ARMA model (in-sample)
difarma <- (arma_training$residuals)
(mse_arma <- (sum(difarma^2))/length(difarma))
(mape_arma <- sum(abs(difarma)/abs(training))/length(difarma))

# evaluating the prediction of the ARMA model (out-of-sample)
predict_arma <- predict(arma_training, n.ahead=9)
dif_out_arma <- (test-predict_arma$pred[2:9])

(mse_armap <- (sum(dif_out_arma^2))/length(dif_out_arma))
(mape_armap <- sum(abs(dif_out_arma)/abs(test))/length(dif_out_arma))


prev_arma <- training-arma_training$residuals


M <- cbind(c(mse_charma, mape_charma),
           c(mse_arma ,mape_arma ))
colnames(M) <- c("CHARMA", "ARMA")
rownames(M) <- c("EQM", "MAPE(%)")
M

Mprev <- cbind(c(mse_charmap, mape_charmap*100),
               c(mse_armap,mape_armap*100))
colnames(Mprev)<-c("CHARMA", "ARMA")
rownames(Mprev)<-c("EQM", "MAPE(%)")
Mprev

#-------------------------
# graphic
plot(temp$temperature, col="black", type="l", 
     ylim=c(30, (max(c(training, test))+3)) , axes =T, main="", xlab="Time (months)", ylab="Maximum temperature (°C)")
lines(c(fit_training$fitted, fit_training$forecast), lty = 2, lwd = 1, col="red")
lines(c(prev_arma, predict_arma$pred), lty = 3, lwd = 1, col="blue4")
legend("top",legend=c( "Observations", "CHARMA", "ARMA"),
       pt.bg="white", lty=c(1,2,3), col=c("black","red","blue4"), bty="n" )
fim <- (end(temp$temperature)[1]+end(temp$temperature)[2]/9)-9
abline(v=fim,lty=2)
