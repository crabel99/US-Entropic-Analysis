library(tidyverse)

source("src/api_key_control.R")
service_eia <- "https://api.eia.gov/"
service_fred <- "https://api.stlouisfed.org/fred/"

# You will need to execute createKey(service_name) for each of the services
# You will get a prompt asking for a password. This is the key that you will
# have from the EIA or from FRED. The EIA does not require an account for the
# API, just an email address. FRED requires an account to be able to get a key.

# Uncomment the next two lines the first time running the code:
# createKey(service_eia)
# createKey(service_fred)

eia::eia_set_key(accessKey(service_eia))
fredr::fredr_set_key(accessKey(service_fred))

# eia::eia_cats(711245)
# series_id       name                                                            f     units               updated
# <chr>           <chr>                                                           <chr> <chr>               <chr>
#   1 TOTAL.CDEGRUS.A Total Energy CO2 Emissions per Real Dollar of GDP, Annual       A     Metric Tons Carbon… 24-MAR…
# 2 TOTAL.CDTPRUS.A Total Energy CO2 Emissions per Capita, Annual                   A     Metric Tons Carbon… 24-MAR…
# 3 TOTAL.TEGDSUS.A Energy Expenditures as Share of GDP, Annual                     A     Percent             24-MAR…
# 4 TOTAL.TEGOSUS.A Energy Expenditures as Share of Gross Output, Annual            A     Percent             24-MAR…
# 5 TOTAL.TETCBUS.A Total Primary Energy Consumption, Annual                        A     Trillion Btu        21-APR…
# 6 TOTAL.TETCEUS.A Total Energy CO2 Emissions, Annual                              A     Million Metric Ton… 21-APR…
# 7 TOTAL.TETCHUS.A Energy Expenditures per Capita, Annual                          A     Nominal Dollars     24-MAR…
# 8 TOTAL.TETCVUS.A Energy Expenditures, Annual                                     A     Million Nominal Do… 24-MAR…
# 9 TOTAL.TETGRUS.A Total Primary Energy Consumption per Real Dollar of GDP, Annual A     Thousand Btu per C… 24-MAR…
# 10 TOTAL.TETPRUS.A Total Primary Energy Consumption per Capita, Annual             A     Million Btu         24-MAR…

# Energy Expenditure Nominal
data_raw <- eia::eia_series("TOTAL.TETCVUS.A")
data_ts <- data_raw$data[[1]] %>%
  relocate(value, .after = year) %>%
  # subset(select = -date) %>%
  rename(TOTAL.TETCVUS.A = value)

series_id <- c("TOTAL.TEACBUS.A", "TOTAL.TECCBUS.A", "TOTAL.TEICBUS.A",
               "TOTAL.TERCBUS.A", "TOTAL.TETCBUS.A")
for (i in 1:length(series_id)) {
  data_raw <- eia::eia_series(series_id[i])
  data_ts[series_id[i]] <- data_raw$data[[1]]$value
}
data_ts <- data_ts %>% arrange(-desc(date))

# Utility
# Sector efficiency taken from https://flowcharts.llnl.gov/commodities/energy
data_ts$VALUE <- (data_ts$TOTAL.TEACBUS.A * 0.21 +
  data_ts$TOTAL.TECCBUS.A * 0.65 +
  data_ts$TOTAL.TEICBUS.A * 0.49 +
  data_ts$TOTAL.TERCBUS.A * 0.65) * 1055.056 / 1000 # PJ


# Marginal Utility
data_ts$MARGINAL.VALUE <- (
  data_ts$VALUE / data_ts$TOTAL.TETCVUS.A  # GJ/$
)

# Nominal GDP
data_raw <- fredr::fredr(
  series_id = "GDPA",
  observation_start = as.Date(data_ts$date[1]),
  observation_end = as.Date(tail(data_ts$date,1))
)
data_ts$GDPA <- data_raw$value

# GDP Deflator
data_raw <- fredr::fredr(
  series_id = "A191RD3A086NBEA",
  observation_start = as.Date(data_ts$date[1]),
  observation_end = as.Date(tail(data_ts$date,1))
)
data_ts$GDP.DEFLATOR <- data_raw$value

# CPI
data_raw <- fredr::fredr(
  series_id = "USACPIALLAINMEI",
  observation_start = as.Date(data_ts$date[1]),
  observation_end = as.Date(tail(data_ts$date,1))
)
data_raw <- c(rep(NA, length(data_ts$year) - length(data_raw$value)),
              data_raw$value)
data_ts$CPI <- data_raw

# M2
data_raw <- fredr::fredr(
  series_id = "M2NS",
  observation_start = as.Date(data_ts$date[1]),
  observation_end = as.Date(tail(data_ts$date,1))
)
data_raw$year <- format(data_raw$date, format = "%Y")
data_raw <- aggregate(value~year, data_raw, mean)
data_raw <- c(rep(NA, length(data_ts$year) - length(data_raw$value)),
              data_raw$value)
data_ts$M2NS <- data_raw

# Energy Price Index
data_ts$EPI <- (rep(data_ts$MARGINAL.VALUE[22], length(data_ts$year)) /
                  data_ts$MARGINAL.VALUE)

# GDPA              Billions of dollars
# TOTAL.TEGDSUS.A   Percent
# TOTAL.TETCBUS.A   Trillion Btu
# VALUE             PJ
# MARGINAL.VALUE    GJ/$

saveRDS(data_ts, file = "data/raw/EIA-FRED.rds")
rm(data_raw, service_eia, service_fred, createKey, updateKey, accessKey, i)

# Select the price indexes, normal them, and plot them
data_index <- data_ts[!is.na(data_ts$EPI), ]
data_index$GDP.DEFLATOR <-
  data_index$GDP.DEFLATOR/rep(data_index$GDP.DEFLATOR[1],
                              length(data_index$year))
data_index$CPI <- data_index$CPI/rep(data_index$CPI[1],length(data_index$year))
data_index$M2NS <- data_index$M2NS/rep(data_index$M2NS[1],length(data_index$year))

data_index <- data_index %>%
  subset(select = c("year", "M2NS", "EPI","CPI","GDP.DEFLATOR")) %>%
  rename(GDP = GDP.DEFLATOR) %>%
  rename(M2 = M2NS) %>%
  as.data.frame() %>%
  reshape::melt(id.vars = "year") %>%
  rename(Year = year, Deflator = variable, Value = value)

saveRDS(data_index, file = "data/processed/Deflators.rds")

plot_index <- data_index %>%
  ggplot(aes(x = Year, y = Value, color = Deflator)) +
  geom_line( size = 1) +
  labs(title = "Comparison of Various Economic Deflators",
       caption = "Data Source: M2, CPI, and GDP taken from FRED API.\
       EPI is derived from EIA data for total energy consumed.")

saveRDS(plot_index, file = "data/processed/Deflator Plot.rds")
print(plot_index)

plot_marginal_utlity <- data_ts[!is.na(data_ts$MARGINAL.VALUE),] %>%
  ggplot(aes(x = year, y = MARGINAL.VALUE * 1000)) +
  geom_line(size = 1) +
  scale_y_continuous(trans = 'log10') +
  labs(x = "Year",
       y = "Marginal Utility [GJ/$]",
       title = "Marginal Utility of the Dollar [MJ/$]",
       caption = "Derived from EIA primary energy data.")

saveRDS(plot_marginal_utlity, file = "data/processed/Marginal Utility Plot.rds")
print(plot_marginal_utlity)

rm(data_ts, data_index, plot_index, plot_marginal_utlity)
