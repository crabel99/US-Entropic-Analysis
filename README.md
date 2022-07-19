# United States Macroeconomic Analysis

This `R` project looks at the US economy from a statistical mechanical framework.
The project constructs a statistical ensemble using the income distribution of the United States.
The data set comprises income data from 1996-2019.

## Model Data Sources

The model data comes from three different sources.

1. Internal Revenue Service [Income Tax Data ](https://www.irs.gov/statistics/soi-tax-stats-individual-statistical-tables-by-size-of-adjusted-gross-income)
2. Energy Information Agency [Open Data API](https://api.eia.gov/)
3. Federal Reserve Bank of St. Louis [FRED API](https://api.stlouisfed.org/fred/)

To avoid having to obtain the API keys, the data has already been downloaded and resides in the [`data/raw/`](data/raw/) directory.
To avoid having to process
the IRS data, the downloaded `*.xls` files have already been parsed into suitable `R` data structures in [`data/processed/irs/`](data/processed/irs/).

## Using the Model

### Setup
The model is developed using [`R-Studio`](https://www.rstudio.com/products/rstudio/download/) and can be installed from that link.
Once `R-Studio` is installed you will need to install the necessary `CRAN` libraries.

To install the libraries, from the R console run:

``` R
install.packages("renv")
```

Once `renv` is installed execute the following from the R console:

``` R
renv::restore()
```

### Estimating the Income Distribution

The income distribution model regression uses HMC( Hamiltonian Monte Carlo) with NUTS (No U-Turn Sampling) in a Bayesian framework to estimate the model's hyper-parameters.
This is done using `rstan` which is an implementation of [`Stan](https://mc-stan.org/).

There are two different distributions that can be selected:

1. [`gibbs_model.stan`](src/gibbs_model.stan)
2. [`maxwell_model.stan`](src/maxwell_model.stan)

The Gibbs model is taken directly from Banerjee and Yakovenko [(2010)](https://arxiv.org/abs/0912.4898) derivation.
Their derivation uses the Gibbs distribution as the basis of the thermal portion of the income distribution and the Pareto distribution for the epithermal portion,

$$
 \begin{align}
   \tag{1}
   f(u;T,u_0,\alpha) = \frac{e^{- \frac{u_0}{T} \arctan{\frac{u}{u_0}}}}{Z 
   \left(1 + \left( \frac{u}{u_0} \right)^2 \right)^\alpha}. \label{eq:1}
 \end{align}
$$

The Maxwell model is derived from the [Maxwell distribution](https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution#Distribution_for_the_energy) in the same manner that the Gibbs distribution was using the stationary Fokker-Plank equation to include the Pareto portion of the distribution.

The Gibbs distribution  provides the better fit of the Adjusted Gross Income (AGI), while the Maxwell distribution provides the better fit of the Taxable Income (TI).

By default, [`income_model_regression.R`](src/income_model_regression.R) is configured to follow the same method and dataset of Banerjee and Yakoveno (2010).
The regression can be configured to use the other distribution or other datasets, by modifying the two variables, `fitModel` and `fitData`,

``` R
# Select Model and Load Data
fitModel <- "gibbs_model.stan"
# fitModel <- "maxwell_model.stan"
#
# Select the appropriate table either "11_agi", "11_ti", or "21"
table <- "11_agi"
```

Once configured, the R script can be run from the console using,

``` R
source("src/income_model_regression.R")
```

This will generate an output file in `data/processed/stan/[gibbs_model.stan | maxwell_model.stan]` with a name based on the ensemble's data set depending on the configuration.

By default, the regression will run 4 simultaneous chains that are used to test for convergence.
This can be adjusted by setting `nChains` to the desired number of chains.

To examine the Stan output, you will need to install `shinystan`

### Estimating Remaining Ensemble Parameters

The remaining ensemble parameters are estimated using [`aggregate_data.R`](src/aggregate_data.R).
This file is configured similarly to `income_model_regression.R`,

``` R
# Configure the analysis
fitModel <- "gibbs_model.stan"
# fitModel <- "maxwell_model.stan"

# Select the appropriate table either "11_agi", "11_ti", or "21"
table <- "11_agi"
```

Once configured, the R script can be run from the console using,

``` R
source("src/aggregate_data.R")
```

This model takes the marginal utility of money, $\lambda$, and converts the units of the economic temperature, $T_{fit}$, equation $\eqref{eq:1}$, from dollars to energy, $GJ$.

$$
\begin{align}
  \tag{2}
  T = \lambda T_{fit}. \label{eq:2}
\end{align}
$$

Next the model computes the mean income,

$$
\begin{align}
  \tag{3}
  \langle m \rangle = \frac{M}{N} = \frac{\sum_j M_j}{\sum_j N_j}. \label{eq:3}
\end{align}
$$

Where $M$ is the total income of the ensemble, and $N$ is the total number of people in the ensemble.

The script then uses the built in `R` linear model `lm` to fit the following equation,
$$
\begin{align}
  \tag{4}
  s = s_0 + cR \log\left(\frac{T}{T_0}\right) + R \log\left(\frac{\langle m \rangle}{\langle m_0 \rangle}\right). \label{eq:4}
\end{align}
$$
This model is a slight adaptation of Callen's (1985, p. 68) equation 3.38 which is the equation of state for an ideal gas.

Where the ensemble's modeled entropy, $s$, is computed from the appropriate distribution, e.g. equation $\eqref{eq:1}$, $T_0 = 2\,815.8 \left[GJ\right]$, and $\langle m_0 \rangle = 37\,689 \left[\$/person\right]$.

When the package is run, it will provide an output something like

``` bash
Call:
lm(formula = S ~ T + M, data = train)

Residuals:
      Min        1Q    Median        3Q       Max 
-0.051839 -0.018045 -0.002761  0.014310  0.075569 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)    
(Intercept) 11.52009    0.01575 731.298  < 2e-16 ***
T            0.13336    0.03577   3.728  0.00183 ** 
M            0.85757    0.04188  20.474 6.66e-13 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.0314 on 16 degrees of freedom
Multiple R-squared:  0.9679,	Adjusted R-squared:  0.9639 
F-statistic: 241.1 on 2 and 16 DF,  p-value: 1.131e-12
```

The parameters $R$ and $c$ can be recovered by running,

``` R
R <- fit$coefficients[3]
c <- fit$coefficients[2] / R
print(R)
print(c)

# 0.8575689
# 0.1555147
```

## References
* Banerjee, A., & Yakovenko, V. M. (2010). Universal patterns of inequality. New Journal of Physics, 12, 1-25. doi:10.1088/1367-2630/12/7/075032
* Energy Information Agency (2022). Open Data. https://www.eia.gov/opendata/index.php.
* Federal Reserve Bank of St. Louis. (2022). FRED Economic Data API. https://fred.stlouisfed.org/docs/api/fred/. 
* Internal Revenue Service (2022). SOI Tax Stats - Individual Statistical Tables by Size of Adjusted Gross Income. https://www.irs.gov/statistics/soi-tax-stats-individual-statistical-tables-by-size-of-adjusted-gross-income.
* Lawrence Livermore National Laboratory (2022). Estimated U.S. Energy Consumption in 2021. [LLNL-MI-410527](https://flowcharts.llnl.gov/sites/flowcharts/files/2022-04/Energy_2021_United-States_0.png).
