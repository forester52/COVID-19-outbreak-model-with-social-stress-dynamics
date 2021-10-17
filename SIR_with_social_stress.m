%% SIR model for coronavirus outbreak dynamics driven by social stress
% 
% Requires MATLAB R2016b
% (comment out the lines with "colororder" commands for MATLAB 2019a or earlier)
% 
% The script illustrates the paper about COVID-19 outbreak spread model
% Authors: Kastalskiy, I. A., Pankratova, E. V., Mirkes, E. M.,
%          Kazantsev, V. B., and Gorban, A. N.
% 
% A model of the COVID-19 pandemic is proposed, combining the dynamics of social
% stress described by the tools of sociophysics with classical epidemic models.
% The model parameters have been successfully fitted to best match the statistical
% observations of epidemics in different countries of the world.
% 
% Please note that depending on the ranges for parameters (q,a,K2,Io), the fit result
% may differ slightly from that in the paper. To get the exact result, as in Table 1
% in the paper use precise ranges and desired precisions for each country individually.
% Also, the result is influenced by possible data updates in file 'owid-covid-data.xlsx'
% caused by changes in the population of countries over time.
% 
% Innokentiy Kastalskiy, kastalskiy@neuro.nnov.ru
% Copyright 2021
% 
% 
%% Loading COVID-19 data and selecting a country or territory
% 
% Preferred use of the data file format, as in the project "Our World in Data":
% https://ourworldindata.org/coronavirus
% 
% Direct link to Excel file:
% https://covid.ourworldindata.org/data/owid-covid-data.xlsx
% 
% Data repository of Our World in Data project:
% https://github.com/owid/covid-19-data/tree/master/public/data


% ----------------------------------------------------------------------------------------------------
% Enter ISO country codes (three capital letters)
% and corresponding Numbers of days for analysis and simulation

% for batch analysis, uncomment the two following lines:
country_iso_code = {'CHN';'BRA';'COL';'FRA';'DEU';'IND';'IRN';'ISR';'ITA';'RUS';'ESP';'GBR';'USA'};
          N_days = [ 200;  300;  300;  200;  200;  400;  100;  100;  200;  200;  200;  200;  150];

% for single country analysis, uncomment the two following lines:
% country_iso_code = {'RUS'};
% N_days = 200;

% ----------------------------------------------------------------------------------------------------


if ~exist( 'table_data', 'var' )
    clearvars -except country_iso_code N_days
    [ filename, pathname ] = uigetfile( '*.xls*', 'Select COVID-19 data file', 'owid-covid-data.xlsx' );
    cd( pathname );
    table_data = readtable( filename );
end

for ii = 1:length(country_iso_code)

raw_data = table_data( strcmp( table_data.iso_code, country_iso_code{ii} ), : );
countryname = raw_data.location{1};
pop = raw_data.population(1);

%% Setting the parameters of the SIR_SS model

T0 = tic;
disp( [ 'Processing of the coronavirus outbreak in ' countryname ' started at:  ' datestr(fix(clock)) ] )


% -------------------------------------------------------------------------------------------------------
% dt - time step for simulations (you can set 1, 0.5, 0.25, 0.2, 0.125, 0.1, ..., 0.01 or even smaller)
%      decrease in the time step improves the accuracy of calculations, but increases their duration
% thr_cases - threshold for total confirmed cases (TCC), which marks the onset of an outbreak
%           (usually 100 cases)

dt = 1;
thr_cases = 100;

% -------------------------------------------------------------------------------------------------------


arrsize = round((N_days(ii)-1)/dt)+1;
TCC = [ 0; raw_data.total_cases ];
ix = find( TCC >= thr_cases, 1 );
Ndays = min( N_days(ii), length( TCC( ix:end ) ) );
data = TCC( ix:ix+Ndays-1 );
datanorm = norm( data-mean(data) );


% -----------------------------------------------------------------------------------------------------
% N - number of points in parameter space (q, a, K2) along each dimension
% N2 - number of points along the parameter Io
% The preferred choice for N or N2 is (10n+1), where n is an integer
% The computation time grows explosively with increasing N and grows linearly with increasing N2,
% but the more correctly the optimization algorithm works

N = 11;
N2 = 11;


% -----------------------------------------------------------------------------------------------------
% Main parameters of the model
% q - stress response rate (see model description [Kastalskiy et al., 2021]). For European countries,
%     q usually takes a value in the range from several 10k to several 100k: (20..300)*1000
% a - infection rate, usually takes a value in the range 0.1..0.4
% b - recovery rate, we take the value 0.1 as the reciprocal of the characteristic recovery period
%     and as a reference for a
% K2 - exhaustion rate (see model description [Kastalskiy et al., 2021]). For European countries,
%      K2 usually takes a value in the range (4..8)*10e-3
% K3 - relaxation rate to ignorant mode (slow), we take the value 0.01 as a reference for K2
% Io - the initial fraction of infected people in the population


% you can refine these 4 ranges for single country analysis:
q_var = linspace(0,500,N)*1000;
a_var = linspace(0.1,0.4,N);
K2_var = linspace(0.003,0.03,N);
Io_var = linspace(1.0e-6,50.0e-6,N2);

b  = 0.1;
K3 = 0.01;

prec_q = +2;                    % desired accuracy in parameter q (10^...)
prec_a = -4;                    % desired accuracy in parameter a (10^...)
prec_K2 = -5;                   % desired accuracy in parameter K2 (10^...)
prec_Io = -7;                   % desired accuracy in parameter Io (10^...)

% -----------------------------------------------------------------------------------------------------


q_range = log10( (q_var(end)-q_var(1)) / (N-1) );
a_range = log10( (a_var(end)-a_var(1)) / (N-1) );
K2_range = log10( (K2_var(end)-K2_var(1)) / (N-1) );
Io_range = log10( (Io_var(end)-Io_var(1)) / (N2-1) );

Sign = zeros(arrsize,1);
Signx = Sign; Sres = Sign; Sexh = Sign; I = Sign; R = Sign;
Sresx = Signx; Sexhx = Signx; Ix = Signx; Rx = Signx;


%% Performing calculations

h = waitbar( 0, ['Coronavirus outbreak in ' countryname ' is calculated...' newline 'parameters: q: ' num2str(q_var(end)/1000) ...
             'k, a: ' num2str(a_var(1)) ', K2: ' num2str(K2_var(1)) ', Io: ' num2str(Io_var(1)) newline 'precision: q: ' num2str(q_range) ...
             ', a: ' num2str(a_range) ', K2: ' num2str(K2_range) ', Io: ' num2str(Io_range) newline 'R^2 = -Inf'] );

flag2 = true;
x = NaN(1,9);
xc = 0;


% -----------------------------------------------------------------------------------------------------

flag0 = true;       % if 'true', then the optimization process is initiated in the parameter space
                    % (q,a,K2,Io) in order to maximize R^2, otherwise the while loop is executed once

% -----------------------------------------------------------------------------------------------------


flag1 = true;
magn = 4;
c = 0;
while flag2 && ( 0.01*round(100*q_range) > prec_q || 0.01*round(100*a_range) > prec_a || 0.01*round(100*K2_range) > prec_K2 ...
      || 0.01*round(100*Io_range) > prec_Io || ( 0.01*round(100*q_range) == prec_q && 0.01*round(100*a_range) == prec_a ...
      && 0.01*round(100*K2_range) == prec_K2 && 0.01*round(100*Io_range) == prec_Io && flag1 ) )
    
    flag2 = flag0;
    flag1 = false;
    c = c+1;
    
    if c > 1
        [ q_var, flagq ] = localization_for_SIR( magn, q, q_var, prec_q, N );
        [ a_var, flaga ] = localization_for_SIR( magn, a, a_var, prec_a, N );
        [ K2_var, flagK2 ] = localization_for_SIR( magn, K2, K2_var, prec_K2, N );
        [ Io_var, flagIo ] = localization_for_SIR( magn, I(1), Io_var, prec_Io, N2 );
        flag1 = flag1 | flagq | flaga | flagK2 | flagIo;
        
        q_range = log10( (q_var(end)-q_var(1)) / (N-1) );
        a_range = log10( (a_var(end)-a_var(1)) / (N-1) );
        K2_range = log10( (K2_var(end)-K2_var(1)) / (N-1) );
        Io_range = log10( (Io_var(end)-Io_var(1)) / (N2-1) );
    end
    
    R2 = -Inf;
    for n = 1:N2
        if Io_var(n) >= 0
            Ix(1) = Io_var(n);
            Signx(1) = 1 - Ix(1);
            Sresx(1) = 1 - Signx(1) - Ix(1);
            Sexhx(1) = 0;
            Rx(1) = 0;
            for m = 1:N
                if q_var(m) >=0
                    q0 = q_var(m);
                    for k = 1:N
                        if a_var(k) >= 0
                            a0 = a_var(k);
                            for j = 1:N
                                if K2_var(j) >= 0
                                    K20 = K2_var(j);
                                    for i = 2:arrsize
                                        v1 = Signx(i-1);
                                        v2 = Sresx(i-1);
                                        v3 = Sexhx(i-1);
                                        v4 = Ix(i-1);
                                        v5 = Rx(i-1);
                                        Signx(i) = v1 + dt*( -q0*v1*v4^2 - a0*v1*v4 + K3*v3 );
                                        Sresx(i) = v2 + dt*( -K20*v2 + q0*v1*v4^2 );
                                        Sexhx(i) = v3 + dt*( -K3*v3 - a0*v3*v4 + K20*v2 );
                                        Ix(i) = v4 + dt*( -b*v4 + a0*v1*v4 + a0*v3*v4 );
                                        Rx(i) = v5 + dt*b*v4;
                                    end
%                                     cmltv = (Ix+Rx)*pop;
%                                     datafit = cmltv(1:round(1/dt):round((Ndays-1)/dt)+1);
                                    R2curr = 1 - ( norm(data-(Ix(1:round(1/dt):round((Ndays-1)/dt)+1)+Rx(1:round(1/dt):round((Ndays-1)/dt)+1))*pop)/datanorm )^2;
%                                     R2curr = 1 - ( norm(data-(Ix+Rx)*pop)/datanorm )^2;
                                    if R2curr > R2
                                        xc = xc + 1;
                                        R2 = R2curr;
                                        K2 = K20; a = a0; q = q0;
                                        I = Ix; Sign = Signx; Sres = Sresx; Sexh = Sexhx; R = Rx;
                                    end
                                end
                            end
                        end
                    end
                end
            end
            waitbar( n/N2, h, ['Coronavirus outbreak in ' countryname ' is calculated...' newline 'parameters: q: ' num2str(q/1000) ...
                'k, a: ' num2str(a) ', K2: ' num2str(K2) ', Io: ' num2str(I(1)) newline 'precision: q: ' num2str(q_range) ...
                ', a: ' num2str(a_range) ', K2: ' num2str(K2_range) ', Io: ' num2str(Io_range) newline 'R^2 = ' num2str(R2)] ) ;
        end
    end

    flageq = false;
    for n = 1:c-1
        if sum( x(n,1:8) - [q a K2 I(1) q_range a_range K2_range Io_range] ) == 0
            flageq = true;
        end
    end
    if sum( x(end,1:4) - [q a K2 I(1)] ) == 0 && q_range == prec_q && a_range == prec_a && K2_range == prec_K2 && Io_range == prec_Io || flageq
        flag2 = false;
    end
    x(c,:) = [ q a K2 I(1) q_range a_range K2_range Io_range R2 ];
end
close(h)

eval( [ 'x_' country_iso_code{ii} '=x;' ] );
eval( [ 'Sign_' country_iso_code{ii} '=Sign;' ] );
eval( [ 'Sres_' country_iso_code{ii} '=Sres;' ] );
eval( [ 'Sexh_' country_iso_code{ii} '=Sexh;' ] );
eval( [ 'I_' country_iso_code{ii} '=I;' ] );
eval( [ 'R_' country_iso_code{ii} '=R;' ] );

toc(T0)

%% Plotting figures

figure(2*(ii-1)+1)
clf
colororder({'k','#0072BD'})
yyaxis left
set( gca, 'FontSize',14 )
plot( 1:dt:dt*(arrsize-1)+1, [Sign'; Sres'; Sexh'], 'LineWidth',2 )%, 'Color','#7E2F8E' )
axis([ 0.5 round((arrsize-1)*dt)+1.5 0 1 ])
grid on
xlabel(['Days since ' num2str(thr_cases) ' cases confirmed'], 'FontSize',14 )
ylabel('Normalized cases for S', 'FontSize',14 )
title(['Popul.~' num2str(round(pop/1e6)) 'm, Io=' num2str(I(1)) ', q=' num2str(q/1000) 'k, a=' num2str(a) ...
       ', K2=' num2str(K2)], 'FontWeight','normal' , 'FontSize',11 )

yyaxis right
plot( 1:dt:dt*(arrsize-1)+1, I, 'LineWidth',2 )
hold on
plot( 1:dt:dt*(arrsize-1)+1, R+I, '-', 'LineWidth',2 , 'Color',[0.6350 0.0780 0.1840]*1.2 )
ylim([ 0 R(end)+I(end) ])
ylabel( 'Normalized cases for I and CC', 'FontSize',14 )
legend( 'S_i_g_n', 'S_r_e_s', 'S_e_x_h', 'I', 'CC', 'Location','west' , 'FontSize',11 )

saveas( gcf, [countryname '_1.fig'] )


figure(2*ii)
clf
colororder([0 0.4470 0.7410; [0.6350 0.0780 0.1840]*1.2])
yyaxis right
set( gca, 'FontSize',14 )
plot( 2-ix:length(TCC)-ix+1, TCC, '.', 'Color',[1 0.25 0.25] , 'MarkerSize',8 )  
hold on
plot( 1:dt:round(dt*(arrsize-1))+1, (I+R)*pop, '-', 'Color',[0.6350 0.0780 0.1840]*1.2 , 'LineWidth',2 )
axis([ 0.5 round((arrsize-1)*dt)+1.5 0 max(TCC(round((arrsize-1)*dt)+ix),(I(end)+R(end))*pop) ])
xlabel(['Days since ' num2str(thr_cases) ' cases confirmed'], 'FontSize',14 )
ylabel('Total confirmed cases', 'FontSize',14 )
title(['COVID outbreak in ' countryname ',  R^2 = ' num2str(R2)], 'FontWeight','normal' , 'FontSize',11 )

yyaxis left
stem( 3-ix:length(TCC)-ix+1, diff(TCC), 'Marker','none' , 'Color','#4DBEEE' )
hold on
plot( 1:dt:round(dt*(arrsize-1))+1-dt, diff(I+R)*pop/dt, '-', 'LineWidth',2 )
ylim([ 0 max(max(diff(TCC(1:round((arrsize-1)*dt)+ix+1))),max(diff(I+R))*pop/dt) ])
ylabel('Daily new cases', 'FontSize',14 )
legend( 'DNC', 'CC''', 'TCC', 'CC', 'Location','east' , 'FontSize',11 )

saveas( gcf, [countryname '_2.fig'] )

end

save( [ 'SIR_SS_model_results_' num2str(dt) '_' num2str(N) '.mat' ] )


function [ q_var, flagq ] = localization_for_SIR( magn, q, q_var, precq, N )
    flagq = false;
    if q == q_var(1)
            q_var = linspace( q-magn*(q_var(end)-q)/2, q+magn*(q_var(end)-q)/2, N );
            flagq = true;
    else
        if q == q_var(end)
            q_var = linspace( q-magn*(q-q_var(1))/2, q+magn*(q-q_var(1))/2, N );
            flagq = true;
        else
            if 0.01*round(100*log10((q_var(end)-q_var(1))/(N-1))) > precq
                m = (N-1)*(10^(ceil(0.01*round(100*log10((q_var(end)-q_var(1))/(N-1))))))/20;
                q_var = linspace( round(q/10^(precq+1))*10^(precq+1)-m, round(q/10^(precq+1))*10^(precq+1)+m, N );
            end
        end
    end
end