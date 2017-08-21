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

#arma_model = armaFit(formula = ~arma(2,2), data = dataset.return.train)


#Finding GARCH Lag order using ACF
#garch_model = garchFit(~arma(0, 2) + garch(1, 1), data = dataset.return.train)
#garch_model@fit$ics[["AIC"]]

best_r = 1
best_s = 1
formula = as.formula( paste( sep="", "~arma(", best_p, ",", best_q, ") + garch(", best_r, ",", best_s,")" ) )
best_model = garchFit(formula = formula, data = dataset.return.train, trace = FALSE)
best_aic = best_model@fit$ics[["AIC"]]
cat(sprintf("GARCH(%d, %d)\t: %f\n", best_r, best_s, best_aic))
for(r in 1 : 5){
	for(s in 1 : 5){
		formula = as.formula( paste( sep="", "~arma(", best_p, ",", best_q, ") + garch(", r, ",", s,")" ) )
		garch_model = garchFit(formula = formula, data = dataset.return.train, trace = FALSE)
		aic = garch_model@fit$ics[["AIC"]]
		cat(sprintf("ARMA(%d, %d)-GARCH(%d, %d)\t: %f\n", best_p, best_q, r, s, aic))
		if(aic < best_aic){
			best_model = garch_model
			best_r = r
			best_s = s
			best_aic = aic
		}
	}
}

cat(sprintf("Best GARCH Model : ARMA(%d, %d)-GARCH(%d, %d)\t: %f\n", best_p, best_q, best_r, best_s, best_aic))

#USE best_p, best_q, best_r, best_s to build model

#Building model using fGarch
formula = as.formula( paste( sep="", "~arma(", best_p, ",", best_q, ") + garch(", best_r, ",", best_s,")" ) )
best_model = garchFit(formula = formula, data = dataset.return.train, trace = TRUE)

#Building model using rugarch
#arch.spec.norm = ugarchspec(variance.model=list(garchOrder=c(best_r, best_s)),mean.model= list(armaOrder=c(best_p, best_q)),distribution.model = "norm")
#arch.fit = ugarchfit(spec=arch.spec.norm, data=dataset.return.train, solver.control=list(trace = 1))

#Use any of the above two models for prediction. Preferably use fgarch because we've used that for finding the best model




