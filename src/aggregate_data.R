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
  econParameters$M[i] <- sum(binIncome[[i]])
  econParameters$N[i] <- sum(binReturns[[i]])
  econParameters$S[i] <- integrate(integrand, 0.01, 1e9,
                                   T = regressors[i,1],
                                   r0 = regressors[i,2],
                                   alpha = regressors[i,3],
                                   c = regressors[i,4])[1] %>% as.numeric()
}

# Fit for the ideal money constant
specEcon <- econParameters # Specific values (per person) for extensive vars
specEcon$M <- specEcon$M / specEcon$N # [k$/person]
specEcon$E <- specEcon$E / specEcon$N # [GJ/person]
specEcon$Tstar <- specEcon$M / (specEcon$S - log(regressors$C)) #[k$]
fit_R <- lm(M ~ Tstar-1, specEcon)
print(summary(fit_R))

params <- list()
params$R <- fit_R$coefficients[[1]] # Ideal Money Constant

# Fit for the specific heat/value capacity
data.c <- data.frame(s = specEcon$S - specEcon$S[1],
                     m = params$R * log(specEcon$M / specEcon$M[1]),
                     e = params$R * log(specEcon$E / specEcon$E[1]))
fit_c <- lm(s - m ~ e-1, data.c)
print(summary(fit_c))
params$c <- fit_c$coefficients[[1]]
params$s0 <- specEcon$S[1]
params$e0 <- specEcon$E[1]
params$m0 <- specEcon$M[1]


# Remaining Economic Parameters
specEcon$T <- specEcon$E / (params$c * params$R) # [GJ]
specEcon$F <- -specEcon$T * log(regressors$C) # [GJ/person]
specEcon$P <- specEcon$T / specEcon$Tstar  # [GJ/k$]
specEcon$H <- specEcon$E + specEcon$P * specEcon$M # [GJ/person]
specEcon$mu <- specEcon$T * ((params$c + 1) * params$R - specEcon$S) # [GJ/person]
specEcon$year <- years

# Generate the composite equation of state
model <- params$s0 + params$R * (params$c * log(specEcon$E / params$e0) +
                                   log(specEcon$M / params$m0))
res <- specEcon$S - model

# Evaluate the fit
R_squared <- 1 - sum(res^2) / sum((specEcon$S - mean(specEcon$S))^2)
t_test <- list(
  #t.test(model - params$s_0, mu = mean(specEcon$S)),
  t.test(model - params$c * params$R * log(specEcon$E / params$e0),
         mu = mean(specEcon$S)),
  t.test(model - params$R * log(specEcon$M / params$m0) ,
         mu = mean(specEcon$S))
)
names(t_test) <- c("c", "R")

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
  geom_line( linewidth = 1) +
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
# params$gamma <- unname(params$gamma, force = TRUE)

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


T_m_data <- fit_R$model
T_m_data$Fit <- fit_R$fitted.values

T_m_plot <- T_m_data %>%
  ggplot(aes(x = Tstar, y = M)) +
  geom_point(size = 1) +
  geom_line(data = T_m_data, aes(x = Tstar, y = Fit), color = 'red') +
  labs(x = TeX("$T_m$ [k\\$]"),
       y = TeX("$\\bar{m}$ [k\\$/person]"),
       title = TeX("$\\bar{m}-T_m$ Plot US Economy"),
       subtitle = "1996-2019")
ggsave("plots/m-T*.pdf",
       width = 6, height = 4, units = "in")
ggsave("plots/m-T*.jpg",
       width = 6, height = 4, units = "in")

e_T_data <- fit_c$model
e_T_data$T <- specEcon$T
e_T_data$Fit <- fit_c$fitted.values
e_T_plot <- e_T_data %>%
  ggplot(aes(x = T, y = e)) +
  geom_point(size = 1) +
  geom_line(data = e_T_data, aes(x = T, y = Fit), color = 'red') +
  labs(x = TeX("$T [GJ]$"),
       y = TeX("$\\bar{s} - R \\log\\left[\\bar{m}/m_0\\right]$"),
       title = TeX("$\\bar{e}-T$ Plot US Economy"),
       subtitle = "1996-2019")
ggsave("plots/e-T.pdf",
       width = 6, height = 4, units = "in")
ggsave("plots/e-T.jpg",
       width = 6, height = 4, units = "in")


EPI_plot <- raw %>% subset(!is.na(MARGINAL.VALUE)) %>%
  ggplot(aes(x = year, y = MARGINAL.VALUE * 1000)) +
  geom_line(color = 'red') +
  scale_y_continuous(trans = 'log10') +
  labs(x = TeX("$year$"),
       y = TeX("Marginal Value [MJ/\\$]$"),
       title = TeX("Marginal Value of the US dollar"),
       subtitle = "1970-2020")

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
file.copy(flist,
          "../../Papers/Resolving\ the\ Allais\ Paradox/images",
          overwrite = TRUE)

rm(i, deflator_names, flist, ind, ind_96, res, years)
