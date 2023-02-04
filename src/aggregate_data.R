library(tidyverse)
library(latex2exp)


# Configure the analysis
fitModel <- "gibbs_model.stan"
# fitModel <- "maxwell_model.stan"

# # Select the appropriate table either "11_agi", "11_ti", or "21"
table <- "11_agi"


# Function Declarations
integrand <- function(x,T,r0,alpha,c) {
  if (fitModel == "gibbs_model.stan") {
    #Gibbs
    v = 1 / (1 + (x / r0)^2)^((1 + alpha) / 2) * exp(-r0/T * atan(x / r0));
  } else {
    #Maxwell
    v = ((1 - x/r0 + (x/r0)^2) / (1 + x/r0)^2)^(r0 / (T * 6)) *
      x / (1 + (x/r0)^3)^((2 + alpha)/3) *
      exp(r0 / (sqrt(3) * T) * atan((1 - 2 * x / r0) / sqrt(3)))
  }
  return(-log(v / c) * v / c)
}

utility <- function(u, m, c, R, s_0, u_0, m_0) {
  return(s_0 + R*(c*log(u / u_0) + log(m / m_0)))
}


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
names(econParameters) <- "P"
econParameters$M <- double(length(years))
econParameters$Tstar <- regressors$T
econParameters$T <- double(length(years))
econParameters$S <- double(length(years))
econParameters$mu <- double(length(years))
econParameters$N <- double(length(years))
econParameters$E <- raw$VALUE[raw$year %in% years] * 1e6
econParameters$F <- double(length(years))
econParameters$H <- double(length(years))


for (i in 1:length(years)) {
  econParameters$M[i] <- sum(binIncome[[i]]) * 1000
  econParameters$N[i] <- sum(binReturns[[i]])
  econParameters$S[i] <- integrate(integrand, 0.01, 1e9,
                                   T = regressors[i,1],
                                   r0 = regressors[i,2],
                                   alpha = regressors[i,3],
                                   c = regressors[i,4])[1] %>% as.numeric()
}

specEcon <- econParameters # Specific values (per person) for extensive vars

# Fit for the ideal money constant
specEcon$M <- specEcon$M / specEcon$N
specEcon$Tstar <- specEcon$M / (specEcon$S - log(regressors$C))
fit_R <- lm(M ~0 + Tstar, specEcon)
print(summary(fit_R))

params <- list()
params$R <- fit_R$coefficients[1] # Ideal Money Constant
params$R <- unname(params$R, force = TRUE)

# Fit for the specific heat/value capacity
specEcon$E <- specEcon$E / specEcon$N
specEcon$T <- specEcon$E / (specEcon$S - log(regressors$C))
fit_c <- lm(E ~ 0 + T , specEcon) # (3.34) from Callen
print(summary(fit_c))

params$c <- fit_c$coefficients[1]/params$R #Specific Heat Capacity
params$c <- unname(params$c, force = TRUE)


# Remaining Economic Parameters
specEcon$F <- -specEcon$T * log(regressors$C) # [GJ/person]
specEcon$P <- specEcon$T / specEcon$Tstar * 1000 # [GJ/k$]
specEcon$M <- specEcon$M / 1000 # [k$/person]
specEcon$H <- specEcon$E + specEcon$P * specEcon$M # [GJ/person]
specEcon$mu <- specEcon$T * ((params$c + 1) * params$R - specEcon$S) # [GJ/person]
specEcon$year <- as.numeric(years)


# Generate the composite equation of state
model <- params$c * params$R * log(specEcon$E / specEcon$E[1]) +
  params$R * log(specEcon$M / specEcon$M[1])
res <- specEcon$S - model
fit_norm <- fitdistrplus::fitdist(res, distr = "norm", method = "mle")
params$s_0 <- fit_norm$estimate[1]
params$s_0 <- unname(params$s_0, force = TRUE)
model <- params$s_0 + model
res <- res - params$s_0


# Evaluate the fit
plot(fit_norm)
summary(fit_norm)
R_squared <- 1 - sum(res^2) / sum((specEcon$S - mean(specEcon$S))^2)
t_test <- list(
  t.test(model - params$s_0, mu = mean(specEcon$S)),
  t.test(model - params$c * params$R * log(specEcon$E / specEcon$E[1]),
         mu = mean(specEcon$S)),
  t.test(model - params$R * log(specEcon$M / specEcon$M[1]) ,
         mu = mean(specEcon$S))
)
names(t_test) <- c("s_0","c", "R")

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
                              "Value" = specEcon$P[1] / specEcon$P))
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
params$gamma <- 1 + params$R / params$c
params$gamma <- unname(params$gamma, force = TRUE)

data_poly <- specEcon %>% subset(select = c(M, P))
data_poly$M <- log(data_poly$M)
data_poly$P <- log(data_poly$P)

fit_poly <- lm(P ~ M, data_poly)
print(summary(fit_poly))

params_poly <- list()
params_poly$C <- exp(fit_poly$coefficients[1])
params_poly$C <- unname(params_poly$C, force = TRUE)
params_poly$n <- -fit_poly$coefficients[2]
params_poly$n <- unname(params_poly$n, force = TRUE)
params_poly$K <- (params_poly$n - params$gamma) / (1 - params$gamma)
params_poly$K <- unname(params_poly$K, force = TRUE)
specEcon$P.fit <- exp(fit_poly$fitted.values)

# Work extracted from society
params_poly$del_w <- params_poly$C / (1 - params_poly$n) *
  (specEcon$M[24]^(1 - params_poly$n) -
     specEcon$M[1]^(1 - params_poly$n))

# Generate summary plots
Ts_plot <- specEcon %>%
  ggplot(aes(x = S, y = T)) +
  geom_point(size = 1) +
  geom_text(aes(label = year)) +
  labs(x = TeX("$\\bar{s} [{person}^{-1}]$"),
       y = "T [GJ]",
       title = TeX("$T-\\bar{s}$ Diagram US Economy"),
       subtitle = "1996-2019")
ggsave("plots/T-s.pdf",
       width = 6, height = 4, units = "in")
ggsave("plots/T-s.jpg",
       width = 6, height = 4, units = "in")

P_m_plot <- specEcon %>%
  ggplot(aes(x = M, y = P)) +
  geom_point(size = 1) +
  geom_text(aes(label = year)) +
  geom_line(data = specEcon, aes(x = M, y = P.fit), color = 'red') +
  labs(x = TeX("$\\bar{m}$ [k\\$/person]"),
       y = TeX("$P$ [GJ/k\\$]"),
       title = TeX("$P-\\bar{m}$ Diagram US Economy"),
       subtitle = "1996-2019")
ggsave("plots/P-m.pdf",
       width = 6, height = 4, units = "in")
ggsave("plots/P-m.jpg",
       width = 6, height = 4, units = "in")

P_H_plot <- specEcon %>%
  ggplot(aes(x = H, y = P)) +
  geom_point(size = 1) +
  geom_text(aes(label = year)) +
  labs(x = TeX("$\\bar{h}$ [GJ/person]"),
       y = TeX("$P$ [GJ/k\\$]"),
       title = TeX("$P-\\bar{h}$ Diagram US Economy"),
       subtitle = "1996-2019")
ggsave("plots/P-H.pdf",
       width = 6, height = 4, units = "in")
ggsave("plots/P-H.jpg",
       width = 6, height = 4, units = "in")

mu_N_plot <- specEcon %>%
  ggplot(aes(x = N, y = mu)) +
  geom_point(size = 1) +
  geom_text(aes(label = year)) +
  labs(x = TeX("$N$"),
       y = TeX("$\\mu$ [GJ/person]"),
       title = TeX("$\\mu-N$ Diagram US Economy"),
       subtitle = "1996-2019")
ggsave("plots/mu_N.pdf",
       width = 6, height = 4, units = "in")
ggsave("plots/mu_N.jpg",
       width = 6, height = 4, units = "in")


T_m_data <- fit_R$model / 1000
T_m_data$Fit <- fit_R$fitted.values / 1000

T_m_plot <- T_m_data %>%
  ggplot(aes(x = Tstar, y = M)) +
  geom_point(size = 1) +
  geom_line(data = T_m_data, aes(x = Tstar, y = Fit), color = 'red') +
  labs(x = TeX("$T^*$ [k\\$]"),
       y = TeX("$\\bar{m}$ [k\\$/person]"),
       title = TeX("$\\bar{m}-T^*$ Plot US Economy"),
       subtitle = "1996-2019")
ggsave("plots/m-T*.pdf",
       width = 6, height = 4, units = "in")
ggsave("plots/m-T*.jpg",
       width = 6, height = 4, units = "in")

e_T_data <- fit_c$model
e_T_data$Fit <- fit_c$fitted.values
e_T_plot <- e_T_data %>%
  ggplot(aes(x = T, y = E)) +
  geom_point(size = 1) +
  geom_line(data = e_T_data, aes(x = T, y = Fit), color = 'red') +
  labs(x = TeX("$T [GJ]$"),
       y = TeX("$\\bar{e} [GJ\\cdot{person}^{-1}]$"),
       title = TeX("$\\bar{e}-T$ Plot US Economy"),
       subtitle = "1996-2019")
ggsave("plots/e-T.pdf",
       width = 6, height = 4, units = "in")
ggsave("plots/e-T.jpg",
       width = 6, height = 4, units = "in")


# Print output plots
print(Ts_plot)
print(mu_N_plot)
print(P_m_plot)
print(P_H_plot)
print(T_m_plot)
print(e_T_plot)
print(plot_index)

# Move plots to paper directory
flist <- list.files("plots", "^.+[.]pdf$", full.names = TRUE)
# file.copy(flist,
#           "../../Papers/Quantum\ Foundations\ of\ Utility/images",
#           overwrite = TRUE)

rm(i, deflator_names, flist, ind, ind_96, res, years)
