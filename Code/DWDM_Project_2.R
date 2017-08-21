#Load reauired libraries
library(quantmod)
library(fArma)
library(fGarch)
library(rugarch)

#Load Dataset
dataset <- read.csv("../Data/FreePhysicalMemory.csv",TRUE)
dataset = dataset[, 2]
dataset.return = diff(dataset)

#Split test and train
no_of_samples = length(dataset.return)
no_of_train_samples = as.integer(no_of_samples * 0.75)
no_of_test_samples = no_of_samples - no_of_train_samples
dataset.return.train = head(dataset.return, no_of_train_samples)
dataset.return.test = tail(dataset.return, no_of_test_samples)
#dataset.return.test = dataset.return[no_of_train_samples + 1 : no_of_samples, 1]

#ACF PCF Analysis on the data
acf(dataset.return.train, lag.max = 100, main= "ACF")
pacf(dataset.return.train, lag.max = 100 , main = "PACF")

#Convert to Time Series
dataset.return.train = as.ts(dataset.return.train)
dataset.return.test = as.ts(dataset.return.test)

#Fidning ARMA Lag order using ACF
best_p = 0
best_q = 0
best_model = armaFit(formula = ~arma(0, 0), data = dataset.return.train)
best_aic = best_model@fit$aic
#cat(sprintf("ARMA(%d, %d)\t: %f\n", best_p, best_q, best_aic))
for(p in 0 : 5){
    for(q in 0 : 5){
        formula = as.formula( paste( sep="", "~arma(", p, ",", q, ")" ) )
        arma_model = armaFit(formula = formula, data = dataset.return.train)
        aic = arma_model@fit$aic
        cat(sprintf("ARMA(%d, %d)\t: %f\n", p, q, aic))
        if(aic < best_aic){
            best_model = arma_model
            best_p = p
            best_q = q
            best_aic = aic
        }
    }
}

cat(sprintf("Best ARMA Model - ARMA(%d, %d)\t: %f\n", best_p, best_q, best_aic))

formula = as.formula( paste( sep="", "~arma(", best_p, ",", best_q, ")" ) )
best_arma_model = armaFit(formula = formula, data = dataset.return.train)
Box.test(residuals(best_arma_model) , lag = 12 , type ="Ljung-Box")
Box.test(residuals(best_arma_model) , lag = 24 , type ="Ljung-Box")
Box.test(residuals(best_arma_model) , lag = 36 , type ="Ljung-Box")

best_aic1=999
for(p in 1 : 5){
    for(q in 1 : 5){
                    arch.spec.std = ugarchspec(variance.model=list(garchOrder=c(p,q)),mean.model= list(armaOrder=c(best_p,best_q)))
                    arch.fit.std = ugarchfit(spec=arch.spec.std, data=dataset.return.train,solver.control=list(trace = 1))
                    aic= infocriteria(arch.fit.std)[1]
                    if(aic < best_aic1){
                            best_model = arma_model
                            best_p1 = p
                            best_q1 = q
                            best_aic1 = aic
                                    }
                    }
        }
        
cat(sprintf("Best GARCH Model : ARMA(%d, %d)-GARCH(%d, %d)\t: %f\n", best_p, best_q, best_p1, best_q1, best_aic1))

arch.spec = ugarchspec(variance.model=list(garchOrder=c(best_p1,best_q1)),mean.model= list(armaOrder=c(best_p,best_q)))
arch.fit = ugarchfit(spec=arch.spec, data=dataset.return.train,solver.control=list(trace = 0))
#sqrt( mean( (dataset.return.test-arch.fit.std.pred_values)^2 , na.rm = TRUE ) )

#forecast
model_forecast=ugarchforecast(arch.fit, data = NULL, n.ahead = no_of_test_samples, n.roll= 0, out.sample = 0)

#rolling forecast
rollmodel=ugarchspec (
variance.model = list(model = "sGARCH", garchOrder = c(best_p1, best_q1)),
mean.model = list(armaOrder = c(best_p,best_q)),distribution.model = "norm")
modelfit=ugarchfit(rollmodel,data=dataset.return.train)
model_forecast_rolling=ugarchforecast(modelfit, data = NULL, n.ahead = no_of_test_samples, n.roll= 0)   
#plot(modelfor,dataset.return.test) 


#Calculating RMSE and MAE
residuals = dataset.return.test - fitted(model_forecast)
squared_residuals = residuals^2
root_mean_squared_ploterror = sqrt(mean(squared_residuals))
mean_absolute_error = mean(abs(residuals))


#testing for first 50 values
root_mean_squared_error = sqrt(mean(head(squared_residuals, 40)))
mean_absolute_error = mean(head(squared_residuals, 40))
