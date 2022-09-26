library(tidyverse)
library(latex2exp)
# library(rayshader)
# library(magick)
# library(hrbrthemes)


integrand <- function(x,T,r0,alpha,c) {
  #Gibbs
  v = 1 / (1 + (x / r0)^2)^((1 + alpha) / 2) * exp(-r0/T * atan(x / r0));
  #Maxwell
  # v = ((1 - x/r0 + (x/r0)^2) / (1 + x/r0)^2)^(r0 / (T * 6)) *
  #   x / (1 + (x/r0)^3)^((2 + alpha)/3) *
  #   exp(r0 / (sqrt(3) * T) * atan((1 - 2 * x / r0) / sqrt(3)))
  return(-log(v / c) * v / c)
}

utility <- function(u, m, c, R, s_0, u_0, m_0) {
  return(s_0 + R*(c*log(u / u_0) + log(m / m_0)))
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

years <- names(binReturns)

econParameters <- tibble(year = years, raw[raw$year %in% years, 10]) %>%
  column_to_rownames(var = "year")
names(econParameters) <- "lambda"
econParameters$M <- double(length(years))
econParameters$N <- double(length(years))
econParameters$T <- regressors$T
econParameters$S <- double(length(years))
econParameters$U <- raw$UTILITY[raw$year %in% years] * 1e6

for (i in 1:length(years)) {
  econParameters$M[i] <- sum(binIncome[[i]]) * 1000
  econParameters$N[i] <- sum(binReturns[[i]])
  econParameters$S[i] <- integrate(integrand, 0.01, 1e9,
                                   T = regressors[i,1],
                                   r0 = regressors[i,2],
                                   alpha = regressors[i,3],
                                   c = regressors[i,4])[1] %>% as.numeric()
}

data_input <- econParameters

# Fit the specific money supply and and temperature to specific entropy
data_input$M <- data_input$M / data_input$N
data_input$U <- data_input$U / data_input$N
fit_m_T <- lm(M ~0 + T, data_input)
print(summary(fit_m_T))

#Normalize the specific extensive parameters to an initial value
data_input$U <- data_input$U / data_input$U[1]
data_input$M <- data_input$M / data_input$M[1]
data_input[, c("U", "M")] <- log(data_input[, c("U", "M")])

# Remove the contribution of money to entropy
data_input$S <- data_input$S - fit_m_T$coefficients[1] * data_input$M

# split <- caTools::sample.split(data_input$S, SplitRatio = 0.8)
# train <- subset(data_input, split == "TRUE")
# test <- subset(data_input, split == "FALSE")

fit <- lm(S ~ U , data_input) # (3.34) from Callen

res <- resid(fit)
fit_norm <- fitdistrplus::fitdist(res, distr = "norm", method = "mle")
# Evaluate the fit
plot(fitted(fit), res)
plot(fit_norm)
print(summary(fit))

#Evaluate the model against the testing data set
# pred <- predict(fit, test, se.fit = TRUE)
# plot(density(pred$fit - test$S))
# pred_r2 <- 1 - sum((pred$fit - test$S)^2) / sum((test$S - mean(test$S))^2)
# pred_t <- (mean(pred$fit) - mean(test$S))/(sqrt(var(pred$fit -
#                                                       test$S)*length(test$S)))
# pred_p <- pt(q = pred_t, df = length(test$S) - 1)

# Generate remaining intensive parameters and collect everything into one
# dataframe
R <- fit_m_T$coefficients[1]
names(R) <- "R"
c <- fit$coefficients[2]/fit_m_T$coefficients[1]
names(c) <- "c"
s_0 <- fit$coefficients[1]
names(s_0) <- "s_0"
specEcon <- econParameters
specEcon$M <- specEcon$M / specEcon$N / 1000
specEcon$U <- specEcon$U / specEcon$N
specEcon$T <- specEcon$U / fit$coefficients[2]
specEcon$lambda <- specEcon$T / econParameters$T * 1000
specEcon$mu <- specEcon$T * ((c + 1) * R - specEcon$S)

# Compare economic deflators
deflators <- readRDS("data/processed/Deflators.rds")

# Normalize to 1996
deflator_names <- unique(deflators$Deflator)
for (i in 1:length(deflator_names)) {
  ind <- deflators$Deflator %in% deflator_names[i]
  ind_96 <- ind & (deflators$Year %in% 1996)
  deflators$Value[ind] <- deflators$Value[ind] / deflators$Value[ind_96]
}
deflators <- rbind(deflators,
                   data.frame("Year" = row.names(specEcon),
                              "Deflator" = "PI",
                              "Value" = specEcon$lambda[1] / specEcon$lambda))
deflators$Year <- as.Date(paste(deflators$Year,"1","1",sep = "-"))

plot_index <- deflators %>%
  ggplot(aes(x = Year, y = Value, color = Deflator, group = Deflator)) +
  geom_line( size = 1) +
  scale_x_date("Year",
               breaks = scales::breaks_width("10 years", offset = "3 years"),
               labels = scales::label_date_short()) +
  labs(title = "Comparison of Various Economic Deflators",
       caption = "Data Source: M2, CPI, and GDP taken from FRED API.\
       EPI is derived from EIA data for total energy consumed.")

ggsave("plots/deflators.pdf",
       width = 10, height = 5.625, units = "in")
ggsave("plots/deflators.jpg",
       width = 10, height = 5.625, units = "in")

# Polytropic Analysis
gamma <- 1 + R / c
names(gamma) <- "gamma"
data_poly <- specEcon %>% subset(select = c(M, N, lambda))
data_poly$M <- log(data_poly$M)
data_poly$lambda <- log(data_poly$lambda)

fit_poly <- lm(lambda ~ M, data_poly)
print(summary(fit_poly))

C <- exp(fit_poly$coefficients[1])
names(C) <- "C"
n <- -fit_poly$coefficients[2]
names(n) <- "n"
K <- (n - gamma) / (1 - gamma)
names(K) <- "K"
specEcon$lambda.fit <- exp(fit_poly$fitted.values)
Const <- specEcon$lambda.fit[1]*specEcon$M[1]^(-fit_poly$coefficients[2])
del_w <- Const / (1 + fit_poly$coefficients[2]) *
  (specEcon$M[24]^(1 + fit_poly$coefficients[2]) -
     specEcon$M[1]^(1 + fit_poly$coefficients[2]))

# Generate summary plots
Ts_plot <- specEcon %>%
  ggplot(aes(x = S, y = T)) +
  geom_line(size = 1) +
  labs(x = TeX("$\\bar{s} [{person}^{-1}]$"),
       y = "T [GJ]",
       title = TeX("$T-\\bar{s}$ Diagram US Economy"),
       subtitle = "1996-2019")
ggsave("plots/T-s.pdf",
       width = 6, height = 4, units = "in")
ggsave("plots/T-s.jpg",
       width = 6, height = 4, units = "in")

lambda_m_plot <- specEcon %>%
  ggplot(aes(x = M, y = lambda)) +
  geom_line(size = 1) +
  geom_line(data = specEcon, aes(x = M, y = lambda.fit), color = 'red') +
  labs(x = TeX("$\\bar{m}$ [k\\$/person]"),
       y = TeX("$P$ [GJ/k\\$]"),
       title = TeX("$P-\\bar{m}$ Diagram US Economy"),
       subtitle = "1996-2019")
ggsave("plots/lambda-m.pdf",
       width = 6, height = 4, units = "in")
ggsave("plots/lambda-m.jpg",
       width = 6, height = 4, units = "in")

mu_N_plot <- specEcon %>%
  ggplot(aes(x = N, y = mu)) +
  geom_line(size = 1) +
  labs(x = TeX("$N$"),
       y = TeX("$\\mu$ [GJ/person]"),
       title = TeX("$\\mu-N$ Diagram US Economy"),
       subtitle = "1996-2019")
ggsave("plots/mu_N.pdf",
       width = 6, height = 4, units = "in")
ggsave("plots/mu_N.jpg",
       width = 6, height = 4, units = "in")


T_m_data <- fit_m_T$model / 1000
T_m_data$Fit <- fit_m_T$fitted.values / 1000

T_m_plot <- T_m_data %>%
  ggplot(aes(x = T, y = M)) +
  geom_point(size = 1) +
  geom_line(data = T_m_data, aes(x = T, y = Fit), color = 'red') +
  labs(x = TeX("$T^*$ [k\\$]"),
       y = TeX("$\\bar{m}$ [k\\$/person]"),
       title = TeX("$\\bar{m}-T^*$ Plot US Economy"),
       subtitle = "1996-2019")
ggsave("plots/m-T*.pdf",
       width = 6, height = 4, units = "in")
ggsave("plots/m-T*.jpg",
       width = 6, height = 4, units = "in")

s_u_data <- fit$model
s_u_data$Fit <- fit$fitted.values
s_u_plot <- s_u_data %>%
  ggplot(aes(x = U, y = S)) +
  geom_point(size = 1) +
  geom_line(data = s_u_data, aes(x = U, y = Fit), color = 'red') +
  labs(x = TeX("$\\log \\left(\\bar{e}/\\bar{e}_0\\right)$"),
       y = TeX("$\\bar{s}{\\prime} [{person}^{-1}]$"),
       title = TeX("$\\bar{s}\\prime-\\log\\left(\\frac{\\bar{e}}{\\bar{e}_0}
                   \\right)$ Plot US Economy"),
       subtitle = "1996-2019")
ggsave("plots/s*-u.pdf",
       width = 6, height = 4, units = "in")
ggsave("plots/s*-u.jpg",
       width = 6, height = 4, units = "in")


# Utility Surface Plot
# len <- 100
# x <- rep(0, len)
# y <- rep(0, len)
# x[1] <- round(min(specEcon$U) / 1.1)
# x[len] <- round(max(specEcon$U) * 1.1)
# del_x <- (x[len] - x[1]) / (len - 1)
# y[1] <- round(min(specEcon$M) / 1.1)
# y[len] <- round(max(specEcon$M) * 1.1)
# del_y <- (y[len] - y[1]) / (len - 1)
# for (i in 2:len) {
#   x[i] <- x[i - 1] + del_x
#   y[i] <- y[i - 1] + del_y
# }

# data <- expand.grid(U = x, M = y)
# data$S <- utility(data$U, data$M, c, R, s_0, specEcon$U[1], specEcon$M[1])
#
# s_surface <- data %>% ggplot(aes(x = U, y = M, fill = S)) +
#   geom_tile() +
#   geom_line(specEcon, mapping = aes(U, M, color = S), size = 3) +
#   scale_color_continuous(limits = c(min(data$S), max(data$S)))
#
#
# plot_gg(s_surface, width = 5, height = 5, multicore = TRUE, scale = 250,
#         zoom = 0.7, theta = 10, phi = 30, windowsize = c(800, 800))
# Sys.sleep(0.2)
#
# render_snapshot(clear = TRUE)


# Print output plots
print(Ts_plot)
print(mu_N_plot)
print(lambda_m_plot)
print(T_m_plot)
print(s_u_plot)
print(plot_index)

# Move plots to paper directory
flist <- list.files("plots", "^.+[.]pdf$", full.names = TRUE)
file.copy(flist,
          "../../Papers/Quantum\ Foundations\ of\ Utility/images",
          overwrite = TRUE)
