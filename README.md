# composite-quantile-factor
The source code for simulations and real data analysis included in the paper 'A Data-Adaptive Factor Model Using Composite Quantile Approach' is available. For the real data analysis, two economic datasets are used: the CRSP monthly stock return data and the FRED-MD data.

## Overview
- `DAFM algorithm.R`: It is a code for functions that estimate data-adaptive factor model. 
- `RobPCA, QFM algorithms.R`: It is a code that implements algorithms of robust PCA proposed by He et al.(2022) and quantile factor model proposed by Chen et al.(2021).
- Simulations
  - This folder includes codes that derive the simulation results for the location-shift model and the location-scale-shift model, which are presented in Tables 1 and 2.
- Real data analysis
  - Stock return 
    - `Stock return.R` : We used CRSP stock return data, which is unfortunately, not freely available. However, many institutions have access to WRDS which provides CRSP data. The code assumes that the monthly and daily stock return data, from January 2004 to December 2023, are saved as 'CRSP rawdata(monthly).csv' and 'CRSP rawdata(daily).csv', resepctively. The code derives all results in Tables 3-5.
  - FRED-MD 
    - `2020-03.csv` : Raw FRED-MD data until 2020-03.
    - `Forecast.R` : The main code that forecasts the unemployment rate and derives results in Table 7.
    - `Forecasting functions.R` : It includes functions needed in `Forecast.R`.
    - `Factor numbers` : This code estimates the number of factors from all factor models used in this analysis.
    - `Get xxxx.R` : The factors from different models of each rolling window data are estimated and saved as .rds file. This file presents the code that derives the factors saved as `xxxx.rds`. xxxx can be one of the six factor models used in this analysis.
