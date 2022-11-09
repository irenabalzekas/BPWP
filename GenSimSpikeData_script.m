%% Generate simulated data: One year spike rate timeseries 

% DEFINE BASIC PARAMETERS + PATHS
dest = [];
years = 1 ;
n = 3*24*365*years; % 3 samples per hour for one year 
totalsec = years*60*60*365*24;
dt = round( totalsec / n); %year in seconds/n - % confirm this is indeed 20 min***
subject = 'S1';
maxsigamp = [60,20,10,10,40,25,20];
poi = [1, 7,15, 21, 30, 50, 100]; %  DEFINE PERIODS IN SIM SIGNAL periods (in days) of interest to include in simulated signal
polycoef = [-.005 1*10^-50 10^-100  0 0]; % for first, second, third, fourth, fifth, 
noiseparams = [16,10,10,10,6,4,4]; % MAKE THIS MATCH WITH NUMBER OF POI

% SIMULATE SIGNAL WITH NOISE
Fs = 1/dt; % sampling period and rate are inverse 
f = Fs*(0:n-1)/n;
f_adj = 1./(f*24*60*60); 
foi = [];
for i = 1: length(poi)
    [val, idx] = min(abs(f_adj-poi(i)));
    foi(i) = f(idx);
end 

t = linspace(1,totalsec, n);
basesig = zeros(1,n);
components = [];
noisecomponents = [];
for i= 1:length(foi)
    component = maxsigamp(i).*sin(2*pi*foi(i)*t);
    addednoise = noiseparams(i).*randn(1, n);
    basesig = basesig + component + addednoise;
    components(i,:) = component;
    noisecomponents(i,:) = addednoise;
end
t4pol = [1:n];
polynomialcomponent = polycoef(1)*t4pol + polycoef(2)*t4pol.^2 + polycoef(3)*t4pol.^3 +polycoef(4)*t4pol.^4 +polycoef(5)*t4pol.^5;
simsig = basesig + polynomialcomponent; % adding noise
simsig = simsig + abs(min(simsig));% shift up so all values above zero and keep some around zero

% CWT 
[meanpowerCWTL2norm,perioddays_cwt,WT, fwt] = wavelet_decomp_L2norm(simsig, 1/dt);
meanpowerCWTreg = nanmean(abs(WT).^2,2);

% CREATE STRUCT
simdata.spikes = simsig;
simdata.time = t;
simdata.dt = dt;
simdata.N = n;
simdata.time_days = t ./ (60*60*24);
simdata.periods_days = poi;
simdata.polycoef = polycoef;
simdata.years = years;

% PLOT GENERATED TIMESERIES 
figure
subplot(1,4,[1,2,3]) % timeseries
plot(simdata.time_days,simdata.spikes,'Color','black')
xlabel('Time (days)')
ylabel('Spike rate')
title('Simulated spike rate (spikes per hours) timeseries')
subplot(1,4,[4]) % CWT
plot(perioddays_cwt, meanpowerCWTreg)
xlabel('Period (days)')
ylabel('Average power')
title('CWT spectrum')

savename = strcat(dest,subject,'_simdata.mat');
save(savename, 'simdata')