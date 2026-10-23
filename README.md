# Time-varying transmission dynamics in seasonal influenza

Code accompanying the paper *"Bayesian inference for time-varying transmission in seasonal influenza: An application to Irish hospitalisation data"*. The model treat transmission rate as a Gaussian random walk and inference is perform using fixed-lag particle filter within particle marginal Metropolis–Hastings.

## Repository structure

```
.
├── data/
│   ├── simulated/          # cached synthetic datasets 
│   └── real/                # Laboratory-confirmed inpatient cases in Ireland (2022/2023-2025/2026)
│   
├── figures/
│   ├── real_data/            # figures generated from the real-data analysis
│   └── sim_data/             # figures generated from the simulation study
├── scripts/                   # entry-point scripts that call src/ on data/
├── src/                      # all model / inference code 

└── README.md
```



## `scripts/` (entry points)

| Script | Fits | Data |
|---|---|---|
| `run_multi_simulation.R` | age-stratified EpiSSM | synthetic data, cached to `data/simulated/` |
| `run_ireland_analysis.R` | age-stratified EpiSSM | real age-stratified cases, with 14-day forecast |
|`run_pmmh_ireland.R`      | age-stratified EpiSSM| inference via PMMH |
| `run_1age_covid.R` | homogeneous EpiSSM | real total daily cases (`data/real/Covid_Data_Ireland.csv`) |
| `run_epinow2.R` | Epinow2 Model| real total daily cases (`data/real/Covid_Data_Ireland.csv`) |


