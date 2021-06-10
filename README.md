# SIR model for coronavirus outbreak dynamics driven by social stress

MATLAB R2016b required

The script illustrates the paper about COVID-19 outbreak spread model by authors: Kastalskiy, I. A., Pankratova, E. V., Mirkes, E. M., Kazantsev, V. B., and Gorban, A. N.

A model of the COVID-19 pandemic is proposed, combining the dynamics of social stress described by the tools of sociophysics with classical epidemic models.
The model parameters have been successfully fitted to best match the statistical observations of epidemics in different countries of the world.

Please note that depending on the ranges for parameters (q,a,K2,Io), the fit result may differ slightly from that in the paper.
Use the exact ranges for each country individually. Also, the result depends on the change in population.

Innokentiy Kastalskiy, kastalskiy@neuro.nnov.ru
Copyright 2021



# Loading COVID-19 data and selecting a country or territory

Preferred use of the data file format, as in the project "Our World in Data":
https://ourworldindata.org/coronavirus

Direct link to Excel file:
https://covid.ourworldindata.org/data/owid-covid-data.xlsx

Data repository of Our World in Data project:
https://github.com/owid/covid-19-data/tree/master/public/data



# Setting the parameters of the SIR_SS model

dt - time step for simulations (you can set 1, 0.5, 0.25, 0.2, 0.125, 0.1, ..., 0.01 or even smaller) decrease in the time step improves the accuracy of calculations, but increases their duration
thr_cases - threshold for total confirmed cases (TCC), which marks the onset of an outbreak
          (usually 100 cases)

N - number of points in parameter space (q, a, K2) along each dimension
N2 - number of points along the parameter Io
The preferred choice for N or N2 is (10n+1), where n is an integer.
The computation time grows explosively with increasing N and grows linearly with increasing N2, but the more correctly the optimization algorithm works.



# Main parameters of the model

q - stress response rate (see model description [Kastalskiy et al., 2021]). For European countries, q usually takes a value in the range from several 10k to several 100k: (20..300)*1000
a - infection rate, usually takes a value in the range 0.1..0.4
b - recovery rate, we take the value 0.1 as the reciprocal of the characteristic recovery period and as a reference for a
K2 - exhaustion rate (see model description [Kastalskiy et al., 2021]). For European countries, K2 usually takes a value in the range (4..8)*10e-3
K3 - relaxation rate to ignorant mode (slow), we take the value 0.01 as a reference for K2
Io - the initial fraction of infected people in the population
