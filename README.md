# SIR model for coronavirus outbreak dynamics driven by social stress

A model of the COVID-19 pandemic is proposed, combining the dynamics of social stress described by the tools of sociophysics with classical epidemic models.
The model parameters have been successfully fitted to best match the statistical observations of epidemics in different countries of the world.

The script illustrates the paper about COVID-19 outbreak spread model prepared by authors: Kastalskiy, I. A., Pankratova, E. V., Mirkes, E. M., Kazantsev, V. B., and Gorban, A. N.

Requires MATLAB R2016b or later



## Getting Started

Before runing simulations download following files from this repository: '*SIR_with_social_stress.m*' (code) and '*owid-covid-data.xlsx*' (most current COVID-19 data).

Preferred use of the COVID-19 data file format, as in the project "Our World in Data": https://ourworldindata.org/coronavirus

Direct link to Excel file: https://covid.ourworldindata.org/data/owid-covid-data.xlsx

Data repository of Our World in Data project: https://github.com/owid/covid-19-data/tree/master/public/data



## Setting minor parameters of the SIR_SS model

dt - time step for simulations (you can set 1, 0.5, 0.25, 0.2, 0.125, 0.1, ..., 0.01 or even smaller) decrease in the time step improves the accuracy of calculations, but increases their duration

thr_cases - threshold for total confirmed cases (TCC), which marks the onset of an outbreak
          (usually 100 cases)

N - number of points in parameter space (q, a, K2) along each dimension

N2 - number of points along the parameter Io

The preferred choice for N or N2 is (10n+1), where n is an integer.

The computation time grows explosively with increasing N and grows linearly with increasing N2, but the more correctly the optimization algorithm works.



## Major model parameters

q - stress response rate (see model description [Kastalskiy et al., 2021]). For European countries, q usually takes a value in the range from several 10k to several 100k: (20..300)*1000

a - infection rate, usually takes a value in the range 0.1..0.4

b - recovery rate, we take the value 0.1 as the reciprocal of the characteristic recovery period and as a reference for a

K2 - exhaustion rate (see model description [Kastalskiy et al., 2021]). For European countries, K2 usually takes a value in the range (4..8)*10e-3

K3 - relaxation rate to ignorant mode (slow), we take the value 0.01 as a reference for K2

Io - the initial fraction of infected people in the population



## Run

To run demo simulation execute script '*SIR_with_social_stress.m*' in MATLAB environment.

Please note that depending on the ranges for parameters (q,a,K2,Io), the fit result may differ slightly from that in the paper.
To get the exact result, as in **Table 1** in the paper use precise ranges and desired precisions for each country individually. Also, the result is influenced by possible data updates in file '*owid-covid-data.xlsx*' caused by changes in the population of countries over time.



## Authors

* **Innokentiy A. Kastalskiy** - *Model implementation, writing code* - [iakastalskiy](https://github.com/iakastalskiy)
* **Evgeniya V. Pankratova** - *Model simulation*
* **Evgeny M. Mirkes** - *Model formulation*
* **Victor B. Kazantsev** - *Model simulation*
* **Alexander N. Gorban** - *Model formulation*



## Cite

**Social stress drives the multi-wave dynamics of COVID-19 outbreaks**

Innokentiy A. Kastalskiy, Evgeniya V. Pankratova, Evgeny M. Mirkes, Victor B. Kazantsev, and Alexander N. Gorban
(2021)


## Acknowledgements

We acknowledge Russian Foundation for Basic Research (project No. 20-04-60078).
