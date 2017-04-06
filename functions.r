library(forecast)  # Required for tslm(), forecast(), Acf(), and Pacf()
library(nortest)   # For Kolmogorov-Smirnov Test
library(qpcR)      # For AICc
library(GeneCycle) # For Fisher test

# ==== Data Processing Tools ====
# Helper function: Modulus
mod<-function(x,m)
{
  t1<-floor(x/m)
  return(x-t1*m)
}

# Helper function: Create quarterly data from monthly
quarterize <- function(series) {
  # Set quarter value to value of last month in quarter
  output <- c()
  month <- 1
  
  for (i in seq(1, length(series))) {
    if (mod(month, 3) == 0) {
      output <- append(output, series[month])
    }
    month <- month + 1
  }
  
  return(output)
}

quarterize.avg <- function(series) {
  # Set quarter value to average of monthly values
  output <- c()
  month <- 1
  quarter.sum <- 0
  
  for (i in seq(1, length(series))) {
    quarter.sum <- quarter.sum + series[month]
    
    if (mod(month, 3) == 0) {
      quarter.avg <- quarter.sum/3
      quarter.sum <- 0
      output <- append(output, quarter.avg)
    }
    month <- month + 1
  }
  
  return(output)
}

# ==== Linear Regression Model Functions ====
lm.id <- function(x, y, log=TRUE, poly=TRUE, recip=TRUE) {
  # Used to help identify an appropriate linear regression model
  # x:    Predictor
  # y:    Response
  
  # Unused arguments (for now)
  # log:  If TRUE, regress against the log transformed predictor
  # poly: If TRUE, regress against polynomial predictors
  # reicp: ...
  
  lm.slr <- lm(y ~ x)
  lm.log <- lm(y ~ log(x))
  lm.recip <- lm(y ~ I(1/x))
  lm.poly <- lm(y ~ x + I(x^2) + I(x^3))
  models <- list(lm.slr, lm.log, lm.recip, lm.poly)
  
  par(mfrow=c(2, 2))
  
  # Plot a regression model
  plot.reg <- function(model, main) {
    # main: Main title for plot
    values.xlim = range(x)
    values.ylim = range(y, fitted(model))
    
    plot(x, y, xlim=values.xlim, ylim=values.ylim, main=main)
    par(new=TRUE)
    plot(x, fitted(lm.slr), xlim=values.xlim, ylim=values.ylim,
         type="l", xlab="", ylab="")
  }
  
  plot.reg(lm.slr, main="Simple Linear Regression Model")
  plot.reg(lm.poly, main="Polynomial Regression Model")
  plot.reg(lm.log, main="Log Transformed Predictor")
  plot.reg(lm.recip, main="Reciprocal Transformed Predictor")
  
  # Print summaries
  lapply(models, summary)
  
  # Return lm() function calls in a list for further analysis (if necessary)
  return(models)
}

# ======== ARIMA Model Functions ========
## ==== Helper Functions ====
arma.params <- function(model) {
  # Goal: Return well-formatted string of model parameters
  
  # Grab ARIMA parameters, assuming using arima() or auto.arima()
  # Result: (AR, MA, SAR, SMA, period, non-seasonal diff, seasonal diff)
  params <- model$arma
  
  formatted <- paste(
    "ARIMA (", paste(params[1], params[6], params[2], sep=","), ")(",
    paste(params[3], params[7], params[4], sep=","), ")",
    "[", params[5], "]", sep="")
  
  return(formatted)
}

# Goal: Create well formatted residual plots. Intended to be used with sapply().
resid.plot <- function(model, type, main=NULL, show.params=TRUE,
                       zline=FALSE, ...) {
  # model:        Fitted ARIMA model
  # type:         Plotting function to call, e.g. plot(), Acf(), Pacf(), ...
  # main:         Main title text
  # zline:        Plot abline at a=0, b=0
  # show.params:  Show ARIMA model parameters in the subtitle
  # ...:          Additional plotting parameters
  
  # If title not specified, use name of type of graph
  if (is.null(main)) { main=paste(toupper(type), "of Residuals") }
  
  # Add titles to list of graphing parameters
  if (show.params) { subtitle <- paste(arma.params(model), "on", model$series) }
  else { subtitle <- "" }
  
  params <- append(alist(...), alist(main=main, sub=subtitle))
  
  # Call plotting function
  do.call(what=type, args=append(alist(residuals(model)), params))
  
  # Draw lines
  if (zline) { abline(a=0, b=0) }
}

# Create well-formatted Normal quantile-quantile plots
resid.norm <- function(model, main="Normal Probability Plot of Residuals", ...) {
  do.call(what="resid.plot", args=append(
    alist(model, type="qqnorm", main=main), alist(...)))
  qqline(residuals(model))
}

resid.hist <- function(model, main="Histogram of Residuals", xlab="Residuals", ...) {
  do.call(what="resid.plot", args=append(
    alist(model, type="hist", main=main, xlab=xlab), alist(...)))
}

## ==== Main Functions =====
arma.id <- function(data, plot.data=FALSE, arch=TRUE, d=0, d2=0, lag.max=48) {
  # Used to help identify an appropriate ARMA model
  # plot.data:  Plot undifferenced data
  # arch:       Analyze data for ARCH/GARCH effects
  # d:          First-order difference at lag d
  # d2:         Second-order difference at lag d2
  
  # Plot undifferenced data
  if (plot.data) {
    plot(data)
    abline(a=0, b=0)
  }
  
  # Apply differencing
  if (d > 0) {
    if (d2 > 0) { data <- diff(diff(data, d), d2) }
    data <- diff(data, d) 
  }
  
  # ==== ACF/PACF of Observations ====
  par(mfrow=c(1,2))
  
  Acf(data, lag.max=lag.max)
  Pacf(data, lag.max=lag.max)
  
  # Identify ARCH/GARCH Effects
  if (arch==TRUE) {
    Acf(data^2, lag.max=lag.max)
    Pacf(data^2, lag.max=lag.max)
  }
}

arma.diff <- function(data, lag, title=NULL, ylab="Data") {
  # Difference the data at specified lag, return the results, and plot graphs to help check for trend
  
  if (is.null(title)) { title = paste("Data after Lag-", lag, " Difference", sep="") }
  
  diff.data <- diff(data, lag)
  diff.tslm <- tslm(diff.data ~ trend)
  
  plot(diff.data, main=title, ylab=ylab)
  abline(a=0, b=0)          # Plot zero line
  lines(fitted(diff.tslm), lwd="5", col=rgb(0, 0, 1, 0.5))  # Plot trend line
  
  cat("Variance of the original data: ", var(data), "\n",
      "Variance of the lag-", lag, " differenced data: ", var(diff.data),
      sep="")
  
  return(diff.data)
}

arma.diag <- function(model, lag.max=48) {
  # ==== Plot Residuals vs. Time ====
  resid.plot(model, type="plot", main="Residual Values vs. Time",
             xlab="", ylab="Residuals", zline=TRUE)
  
  # ==== Check for Serial Correlation in Residuals ====
  par(mfrow=c(2,2))
  Acf(residuals(model), lag.max=lag.max)
  Pacf(residuals(model), lag.max=lag.max)
  
  # Normality Checks
  resid.hist(model, breaks=50, show.params=FALSE)
  resid.norm(model, show.params=FALSE)
}

# Formal Statistical Tests
arma.test <- function(model) {
  test.lag <- sqrt(length(fitted(model)))  # Take h = sqrt(number of observations)
  
  # ==== Normality Tests ====
  # Shapiro-Wilks
  nortest1 <- shapiro.test(residuals(model))

  # Anderson-Darling Test
  nortest2 <- ad.test(residuals(model))
  
  # ==== Portmanteau Tests ====
  
  # Grab ARIMA parameters, assuming using arima() or auto.arima()
  # Result: (AR, MA, SAR, SMA, period, non-seasonal diff, seasonal diff)
  params <- model$arma
  
  test.df <- sum(params[1:4])  # Take fitdf = p + q + P + Q
  
  # White Noise Hypothesis with h - p - q degrees of freedom
  box <- Box.test(residuals(model), type=c("Ljung-Box"), lag=test.lag, fitdf=test.df)
  
  # McLeod-Li Test with h degrees of freedom
  ml <- Box.test(residuals(model)^2, type=c("Ljung-Box"), lag=test.lag)
  
  #cat("Box-Pierce on h - p - q = ", test.lag - test.df,"df", "\n")
  #cat("h = sqrt(n): ", test.lag, "\n")
  #cat("AR/MA Coefficients: ", test.df, "\n")
  return(list(nortest1, nortest2, box, ml))
}

arma.forecast <- function(data, forecast, zoom=0, zoom.ylim=NULL) {
  # data:       Original data to plot forecast against
  # forecast:   A list of arguments to be passed to the forecast function
  # zoom:       Zoom in an x-year period around forecast
  # zoom.ylim:  Range for y-axis when zooming
   # TODO: Make this automatic
  
  predicted_values <- do.call(what="forecast", args=forecast)
  
  # ==== Automatically identify start/end of fitted values of forecast ====
  model.start <- min(time(predicted_values$fitted))
  model.freq <- predicted_values$model$arma[5]
  model.end <- max(time(predicted_values$fitted)) + forecast$h * (1/model.freq)
  
  # ==== Plotting Code ====
  if (zoom > 0) {
    main=c("Zoomed-In Forecast")
    
    if (!is.null(zoom.ylim)) {
      plot(predicted_values, main=main,
           xlim=range(model.end - zoom, end=model.end), ylim=range(zoom.ylim))
    }
    else {
      plot(predicted_values, main=main,
           xlim=range(model.end - zoom, end=model.end))
    }
  } else {
    plot(predicted_values)
  }
  
  # Plot actual data
  lines(window(data, start=model.start, end=model.end),
        lwd=5, col=rgb(0, 1, 0, 0.25))
  
  #return(predicted_values)
}

resid.compare <- function(model1, model2, share.axis=TRUE, lag.max=48, bins=50) {
  # Compare residuals of models
  # share.axis:  Make residuals vs. time plots share same y-axis range
  # lag.max:     Lag at which to display ACF/PACF graphs
  # bins:        Number of bars for the Normality histograms
  
  models <- list(model1, model2)
  resids <- lapply(models, FUN=residuals)
  
  # ==== Residuals vs. Time ====
  par(mfrow=c(1,2))
  
  if (share.axis) {
    # Identify range for residuals
    xlim=range(time(residuals(model1)), time(residuals(model2)))
    ylim=range(residuals(model1), residuals(model2))
    sapply(models, "resid.plot", type="plot", xlim=xlim, ylim=ylim, zline=TRUE,
           main="Residual Values vs. Time", xlab="", ylab="Residuals")
  } else {
    sapply(models, "resid.plot", type="plot", zline=TRUE,
           main="Residual Values vs. Time", xlab="", ylab="Residuals")
  }
  
  # ==== PACF, ACF ====
  par(mfrow=c(2,2))
  sapply(models, "resid.plot", type="Acf", lag.max=lag.max)
  sapply(models, "resid.plot", type="Pacf", lag.max=lag.max)
  
  # ==== Normality ====
  par(mfrow=c(2,2))
  sapply(models, "resid.hist", breaks=bins) # Histogram
  sapply(models, "resid.norm")              # Normal QQ-Plots
}

# Compare results of formal statistical tests
test.compare <- function(...) {
  # http://stackoverflow.com/questions/5754367/using-substitute-to-get-argument-name-with
  get.args <- function(a, ...) {
    arg <- deparse(substitute(a))
    dots <- substitute(list(...))[-1]
    c(arg, sapply(dots, deparse))
  }
  
  model.names.temp <- get.args(...)
  
  # ...: Fitted ARIMA models to run residual analysis on
  models <- list(...)
  
  # Ordering error
  # results <- lapply(models, arma.test)
  results <- list()
  
  for (i in seq(1, length(models))) {
    results <- append(results, list(arma.test(models[[i]])))
  }
  
  # Store test results
  statistics <- c()
  parameters <- c()  # Degrees of freedom
  p.values <- c()
  methods <- c()     # Test used
  # data.names <- c()  # Name of the time series
  model.names <- c() # Name of the model
  
  # Loop through models
  for (i in seq(1, length(results))) {
    
    # Loop through tests
    for (j in seq(1, length(results[[i]]))) {
      statistics <- append(statistics, round(results[[i]][[j]]$statistic, 2))
      parameters <- append(parameters, round(as.numeric(results[[i]][[j]]$parameter), 2))
      p.values <- append(p.values, round(results[[i]][[j]]$p.value, 2))
      
      # Add test names
      if (j == 1) { test.name <- "Shapiro-Wilk" }
      if (j == 2) { test.name <- "Anderson-Darling" }
      if (j == 3) { test.name <- "Box-Ljung" }
      if (j == 4) { test.name <- "McLeod-Li" }
      
      methods <- append(methods, test.name)
      
      # methods <- append(methods, results[[i]][[j]]$method)
      # data.names <- append(data.names, results[[i]][[j]]$data.name)
      model.names <- append(model.names, model.names.temp[i])
    }
  }
    
  results.df <- data.frame("Model"=model.names, "Test"=methods, "Statistic"=statistics,
                           "df"=parameters, "p-value"=p.values)
  return(results.df)
}

# ==== Functions for Comparing Fit ====
# Compare log-likelihood, AIC, AICc, etc.
fit.compare <- function(...) {
  get.args <- function(a, ...) {
    arg <- deparse(substitute(a))
    dots <- substitute(list(...))[-1]
    c(arg, sapply(dots, deparse))
  }
  
  model.names <- get.args(...)
  
  # ...: Fitted ARIMA models to analyze
  models <- list(...)
  
  n.coefs <- c()      # Number of total ARMA coefficients
  log.liks <- c()
  aics <- c()
  aiccs <- c()
  sigma2s <- c()
  codes <- c()    # Convergence code: 0 if converged, some other number o/w
  
  # Loop through models
  for (i in seq(1, length(models))) {
    n.coefs <- append(n.coefs, length(models[[i]]$coef))
    log.liks <- append(log.liks, round(models[[i]]$loglik, 2))
    aics <- append(aics, round(models[[i]]$aic, 2))
    aiccs <- append(aiccs, round(AICc(models[[i]]), 2))
    sigma2s <- append(sigma2s, round(models[[i]]$sigma2, 2))
    if (models[[i]]$code == 0) {
      codes <- append(codes, "Yes")
    } else {
      codes <- append(codes, "Error")
    }
  }
  
  results.df <- data.frame(
    "Model"=model.names, "Order"=n.coefs, "Log Lik"=log.liks,
    "AIC"=aics, "AICc"=aiccs, "sigma2"=sigma2s, "Convergence?"=codes)
  
  return(results.df)
}