# Define functions
# AIC compensated for small sample sizes
AICc <- function(k, n, lpd) {
  # k - number of coefficients
  # n - number of data points
  # lpd - log likelihood
  f <-  2 * k * n / (n - k - 1) - 2 * lpd
  return(f)
}

# Load Libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(rstan)
  library(BH)
  library(RcppEigen)
  library(Rcpp)
  library(inline)
  library(fitdistrplus)
  library(modeest)
})

# Select Model and Load Data
fitModel <- "gibbs_model.stan"
# fitModel <- "maxwell_model.stan"
#
# Select the appropriate table either "11_agi", "11_ti", or "21"
table <- "11_agi"
binMin <- readRDS(paste("data/processed/irs/binMin_", table, ".rds", sep = ""))
binReturns <- readRDS(
  paste("data/processed/irs/binReturns_", table, ".rds", sep = ""))

# Configure data for processing
years <- names(binReturns)
binVect <- names(binMin)

# Stan's integrate_1d doesn't like the [0,1] endpoint and chokes
if (grepl(table, "11_agi") | grepl(table, "11_ti")) {
  for (i in 1:length(years)) {
    binReturns[[years[i]]][2] <-
      binReturns[[years[i]]][2] + binReturns[[years[i]]][1]
    binReturns[[years[i]]] <- binReturns[[years[i]]][-1]
  }
  for (i in 1:length(binVect)) {
    binMin[[binVect[i]]] <- binMin[[binVect[i]]][-2]
  }
}

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

message("Compiling the Stan code modules... ", appendLF = FALSE)
suppressMessages({
  expose_stan_functions(paste("src", fitModel, sep = "/"),
                        show_compiler_warnings = FALSE)
})
message("COMPLETE")

# Configure the environment
nChains <- 4
nCores <- parallel::detectCores()



# Initialization data from prior distributions
initf <- function() {
  list(
    "sigma" = rgamma(1, 2, 8e-6),
    "T" = rgamma(1, 4, 8e-5),
    "r0" = rgamma(1, 4, 4e-5),
    "alpha" = rgamma(1, 2, 1))
}

# Output prototypes
stanFit <- list()
fitEval <- data.frame(sse = double(),
                      r2 = double(),
                      aic = double())
regressors <- data.frame("C" = double(),
                         "T" = double(),
                         "r0" = double(),
                         "alpha" = double())
posteriorDistributions <- data.frame("sigma.shape" = double(),
                                     "sigma.rate" = double(),
                                     "T.shape" = double(),
                                     "T.rate" = double(),
                                     "r0.shape" = double(),
                                     "r0.rate" = double(),
                                     "alpha.shape" = double(),
                                     "alpha.rate" = double())
params <- c("sigma", "T", "r0", "alpha", "lpd") # output parameters

# Loop through each of the years
for (i in 1:length(years)) {
  # The IRS data bins vary from year to year
  for (j in length(binVect):1) {
    if ( years[i] >= binVect[j]) {
      bins <- length(binReturns[[i]])
      # Adjust the bins based on the marginal utility of money
      x <-  binMin[[j]] #/ 1000
    }
  }
  y_hat <- vector(mode = "numeric", length = bins)
  y <- binReturns[[i]]
  y <- y/sum(y)
  incomeData <- list('N' = bins,
                     'x' = x,
                     'y' = y)
  # print(incomeData)

  # Process the data
  message(paste("Processing tax data from:",years[i]))
  tryCatch(
    suppressMessages({
      stanFit[[i]] <- stan(
        file = paste("src", fitModel, sep = "/"),
        data = incomeData,
        chains = nChains,
        init = initf,
        pars = params,
        cores = nCores
      )
    }),
    error = function(e) {
      print(paste(years[i]," unable to find a fit."))
    },
    finally = {
      # Condition the results for output
      draws <- extract(stanFit[[i]], pars = params)
      reg <- data.frame(row.names = years[i],
                        "T" = mean(draws$T),
                        "r0" = mean(draws$r0),
                        "alpha" = mean(draws$alpha))
      for (j in 1:bins) {
        if (j < bins) {
          y_hat[j] <- model_int(c(x[j], x[j + 1]), as.numeric(reg[1, ]), double())
        } else {
          y_hat[j] <- model_int(c(x[j], Inf), as.numeric(reg[1, ]), double())
        }
      }
      reg$C <- sum(y_hat)
      y_hat <- y_hat/reg$C
      fit <- data.frame(row.names = years[i],
                        sse = sum((y - y_hat)^2))
      fit$r2 <- 1 - fit$sse / sum((y - mean(y))^2)
      fit$aic <-
        AICc(3, bins, mlv(draws$lpd, method = "parzen", kernel = "gaussian"))

      # Determine the posterior-prior information distribution
      prior <- list(sigma = fitdist(as.vector(draws$sigma),"gamma"),
                    T = fitdist(as.vector(draws$T),"gamma", method = "mme"),
                    r0 = fitdist(as.vector(draws$r0),"gamma", method = "mme"),
                    alpha = fitdist(as.vector(draws$alpha),"gamma"))

      # export data
      posteriorDistributions <- posteriorDistributions %>%
        rbind(data.frame(row.names = years[i],
                         "sigma.shape" = prior$sigma$estimate[1],
                         "sigma.rate" = prior$sigma$estimate[2],
                         "T.shape" = prior$T$estimate[1],
                         "T.rate" = prior$T$estimate[2],
                         "r0.shape" = prior$r0$estimate[1],
                         "r0.rate" = prior$r0$estimate[2],
                         "alpha.shape" = prior$alpha$estimate[1],
                         "alpha.rate" = prior$alpha$estimate[2]))
      regressors <- rbind(regressors, reg)
      fitEval <- rbind(fitEval, fit)
    }
  )
}

# Export the desired results
regressors %>% saveRDS(paste("data/processed/stan/",fitModel,"/modelRegressors",
                       table, ".rds", sep = ""))
fitEval %>% saveRDS(paste("data/processed/stan/",fitModel,"/modelEvaluation",
                             table, ".rds", sep = ""))
posteriorDistributions %>% saveRDS(paste("data/processed/stan/",fitModel,"/modelPosteriors",
                          table, ".rds", sep = ""))

stanFit %>% saveRDS(paste("data/processed/stan/",fitModel,"/modelStanOutputs",
                            table, ".rds", sep = ""))

# Clean up the workspace
rm(prior, fit, y_hat, reg, draws, initf, y, years, nChains, params,
   i, j, bins, AICc, model_dist, model_int,
   nCores, binVect, x, incomeData, table, binReturns, binMin)
