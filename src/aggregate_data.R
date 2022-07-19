library(tidyverse)


integrand <- function(x,T,r0,alpha,c) {
  #Gibbs
  v = 1 / (1 + (x / r0)^2)^((1 + alpha) / 2) * exp(-r0/T * atan(x / r0));
  #Maxwell
  # v = ((1 - x/r0 + (x/r0)^2) / (1 + x/r0)^2)^(r0 / (T * 6)) *
  #   x / (1 + (x/r0)^3)^((2 + alpha)/3) *
  #   exp(r0 / (sqrt(3) * T) * atan((1 - 2 * x / r0) / sqrt(3)))
  return(-log(v / c) * v / c)
}

# Configure the analysis
fitModel <- "gibbs_model.stan"
# fitModel <- "maxwell_model.stan"

# Select the appropriate table either "11_agi", "11_ti", or "21"
table <- "11_agi"

# Load the data
raw <- readRDS("data/raw/EIA-FRED.rds")
binIncome <- readRDS(
  paste("data/processed/irs/binIncome_", table, ".rds", sep = ""))
binReturns <- readRDS(
  paste("data/processed/irs/binReturns_", table, ".rds", sep = ""))
regressors <- readRDS(
  paste("data/processed/stan/",fitModel,"/modelRegressors",
        table, ".rds", sep = ""))

# rstan::rstan_options(auto_write = TRUE)
# options(mc.cores = parallel::detectCores())
#
# message("Compiling the Stan code modules... ", appendLF = FALSE)
# suppressMessages({
#   rstan::expose_stan_functions(paste("src", fitModel, sep = "/"),
#                                show_compiler_warnings = FALSE)
# })
# message("COMPLETE")

years <- names(binReturns)

econParameters <- tibble(year = years, raw[raw$year %in% years, 10]) %>%
  column_to_rownames(var = "year")
names(econParameters) <- "lambda"
econParameters$M <- double(length(years))
econParameters$mu <- double(length(years))
econParameters$N <- double(length(years))
econParameters$T <- double(length(years))
econParameters$S <- double(length(years))
econParameters$U <- raw$UTILITY[raw$year %in% years] * 1e6

for (i in 1:length(years)) {
  econParameters$M[i] <- sum(binIncome[[i]]) * 1000
  econParameters$N[i] <- sum(binReturns[[i]])
  econParameters$S[i] <- integrate(integrand, 0.01, 1e9, T = regressors[i,1], r0 = regressors[i,2], alpha = regressors[i,3], c = regressors[i,4])

}

econParameters$S <- as.numeric(econParameters$S) * as.numeric(econParameters$N)
econParameters$T <- regressors$T * econParameters$lambda
econParameters$mu <- (econParameters$U - econParameters$S * econParameters$T + econParameters$lambda*econParameters$M)/econParameters$N
econParameters$F <- econParameters$U - econParameters$T*econParameters$S
econParameters$H <- econParameters$U + econParameters$lambda*econParameters$M
econParameters$G <- econParameters$H - econParameters$T*econParameters$S
econParameters$R <- econParameters$U / econParameters$T / econParameters$N

econParameters$U <- econParameters$lambda * econParameters$M
data_input <- econParameters %>% subset(select = -c(F, G, R, mu, H))
# Fit the specific money supply and and temperature to specific entropy
data_input$S <- data_input$S / data_input$N
data_input$M <- data_input$M / data_input$N
data_input$T <- data_input$T / data_input$T[1]
data_input$M <- data_input$M / data_input$M[1]
data_input[, c("U", "M", "N", "T", "lambda")] <- log(data_input[, c("U", "M", "N", "T","lambda")])
split <- caTools::sample.split(data_input$S, SplitRatio = 0.8)

train <- subset(data_input, split == "TRUE")
test <- subset(data_input, split == "FALSE")
fit <- lm(S ~ T + M, train) # (3.38) from Callen
pred <- predict(fit, test, se.fit = TRUE)
res <- resid(fit)
fit_norm <- fitdistrplus::fitdist(res, distr = "norm", method = "mle")
# Evaluate the fit
plot(fitted(fit), res)
plot(fit_norm)
print(summary(fit))
#Evaluate the model agains the testing dataset
plot(density(pred$fit - test$S))
pred_r2 <- 1 - sum((pred$fit - test$S)^2) / sum((test$S - mean(test$S))^2)
# pred_t <- (mean(pred$fit) - mean(test$S))/(sqrt(var(pred$fit - test$S)*length(test$S)))
# pred_p <- pt(q = pred_t, df = length(test$S) - 1)

specEcon <- econParameters
specEcon$M <- specEcon$M / specEcon$N
specEcon$U <- specEcon$U / specEcon$N
specEcon$S <- specEcon$S / specEcon$N
specEcon$F <- specEcon$F / specEcon$N
specEcon$H <- specEcon$H / specEcon$N
specEcon$G <- specEcon$G / specEcon$N


Ts_plot <- econParameters %>%
  ggplot(aes(x = S/N, y = T)) +
  geom_line(size = 1) +
  labs(x = "s",
       y = "T [GJ]",
       title = "T-s Diagram US Economy",
       subtitle = "1996-2019")

lambda_m_plot <- econParameters %>%
  ggplot(aes(x = M/N, y = lambda)) +
  geom_line(size = 1) +
  labs(x = "m [$]",
       y = "lambda [GJ/$]",
       title = "lambda-m Diagram US Economy",
       subtitle = "1996-2019")
print(Ts_plot)
print(lambda_m_plot)
