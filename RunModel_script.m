%% DESCRIPTION: Run BPDN on simulated data with optimal delta parameter, test output for significance, plot overall outputs

subjects = {'S1'};
loc = [];
dest = [];
percent2drop = [];
numrepeats = [];
numsamps = [];
N = [];
sampperdaydesired = 5;

for i = 1: length(subjects)
    
    sub = subjects{i};
    % Load reference data from parameter sweep (N, dt, data, etc.)
    filename = strcat(loc,subject,'/deltaselection/MSE_2D_PercentDrop',num2str(percent2drop),'_repeats',num2str(numrepeats),'_N',num2str(N),'_numsamp',num2str(numsamps),'.mat');
    load(filename);

    % Define basics
    delta = deltaselection.delta;
    maxdegree = deltaselection.maxdegree;
    sampletype = deltaselection.sampletype;
    signal = deltaselection.signal;
    time = deltaselection.time; 
    dt_ogdata = time(2)-time(1);
    dt = deltaselection.dt;
    duration = deltaselection.monthsdata;
    N = deltaselection.N;
    t_full = deltaselection.t_full;
    f_k = deltaselection.f_k;
    desiredperioddays = deltaselection.desiredperioddays;
    dct_basis = deltaselection.dct_basis;
    scalefactor = deltaselection.scalefactor;
    poly_basis = deltaselection.poly_basis;
    % Select data based on desired num per day
    ind = find(deltaselection.sampperdays == sampperdaydesired);

    measureddata_wds = deltaselection.measuredata_init{ind}(:);
    time_wds =  deltaselection.t_init{ind}(:);
    phi = round((time_wds-min(time))/dt) + 1; % indices for when each time data were sampled
    numsamp = deltaselection.numsamps(ind);

    % First, BPDN
    [x_wds, xreshuffled, xpercentiles, z_wds, zreshuffled, zpercentiles, reshuff_optval] = BPDN_wReshuffling_short(delta, measureddata_wds, phi, dct_basis, poly_basis);
    sighits = find(xpercentiles > 99);
    [reconsig_wds] = BPDN_reconsig(f_k, x_wds, scalefactor, maxdegree, z_wds, t_full);

    % Rescale power by frequency (% Change from l1 to l2 norm for cwt.  l2 norm reduces amplitude of high freq data, while l1 does not. (scale is inversely related to frequency))
    x2 = [];
    for i=1:length(x_wds) 
        x2(i) = x_wds(i).*1/sqrt(f_k(i));  
    end
    x_power_rescale = abs(x2);

    % Nick CWT - with L2 adjustment
    [meanpowerCWTL2norm,perioddays_cwt,WT, fwt] = wavelet_decomp_L2norm(signal, 1/dt_ogdata);
    [perioddays_cwt,meanpowerCWTreg,pks,locs,reconsigCWT] = BPDN_cwt2findcycles(signal, dt_ogdata); % redundant, but just for peak finding

    % Modify time scale for plots
    time_full_days = (time - min(time))./(60*60*24);
    t_days = (time_wds - min(time_wds))./ (60*60*24);
    t_recon_days = (t_full*dt) ./ (60*60*24);

    % Saving data for fig
    fig.sampperday = sampperdaydesired;
    fig.x = x_wds;
    fig.xreshuffled = xreshuffled;
    fig.xpercentiles = xpercentiles;
    fig.z = z_wds;
    fig.zpercentiles = zpercentiles;
    fig.reshuff_optval = reshuff_optval;
    fig.measureddata = measureddata_wds;
    fig.time = time_wds;
    fig.sighits = sighits;
    fig.reconsig = reconsig_wds;
    fig.x_power_rescale = x_power_rescale;
    fig.meanpowerCWTreg = meanpowerCWTreg;
    fig.perioddays_cwt = perioddays_cwt;
    fig.WT = WT;
    fig.fwt = fwt;
    fig.reconsigCWT = reconsigCWT;
    fig.pks = pks;
    fig.locs = locs;
    fig.time_full_days = time_full_days;
    fig.t_days = t_days;
    fig.t_recon_days = t_recon_days;
    deltaselection.fig = fig;

    % Append and save
    savename = strcat(dest,subject,'_Outputfor_N',num2str(N),'_sampperday',num2str(sampperdaydesired),'_numsamp',num2str(numsamp),'.mat');
    save(savename,'deltaselection');
end 


%% Plot overall model outputs - spectra, reconstructions, etc. 

% INCLUDED IN PLOT 
% original data
% spectrogram of original data
% wavelet spec of original data
% subsampling of og data samples
% BPDN-based reconstruction
% cwt over BPDN spectral comparison

subject = 'S1';
N = 3000;
sampperdaydesired = [];
numsamp = [];
loc = [];

% filename = strcat(destfolder,'Outputfor',num2str(duration),'_N',num2str(N),'_numsamp',num2str(numsamp),'_forfig.mat');
filename = strcat(loc,subject,'_Outputfor_N',num2str(N),'_sampperday',num2str(sampperdaydesired),'_numsamp',num2str(numsamp),'.mat');
load(filename);

% Initialzie variables
signal = deltaselection.signal;
time_full_days = deltaselection.fig.time_full_days;
t_days = deltaselection.fig.t_days;
measureddata_wds = deltaselection.fig.measureddata;
t_recon_days = deltaselection.fig.t_recon_days;
reconsig_wds = deltaselection.fig.reconsig;
desiredperioddays = deltaselection.desiredperioddays;
x = deltaselection.fig.x;
sighits = deltaselection.fig.sighits;
perioddays_cwt = deltaselection.fig.perioddays_cwt;
meanpowerCWTreg = deltaselection.fig.meanpowerCWTreg;
power = zscore(abs(deltaselection.fig.WT)); 
wave_fam = 1./deltaselection.fig.fwt;
fs = 1 / deltaselection.dt;


fig=figure('Position',[10 10 1500 2000]);
fsize = 10;

% PART 1: BASIC RAW DATA AND CWT SPECTROGRAM 

% original signal 
subplot(4,6,[1,2,3,4])
plot(time_full_days,signal, 'Color','black','DisplayName','raw data')
ylabel('Spike rate')
xlim([0 max(time_full_days)])
colorbar()
ax=gca; ax.YAxis.Exponent = 3;
set(gca,'FontSize',fsize)

% spectrogram for cwt
subplot(4,6,[7,8,9,10])
contourf(time_full_days, log(perioddays_cwt), power,'edgecolor','none')
xlabel('Time (days)')
ylabel('Period (days)')
title({'';'Spectrogram, CWT of raw spike rate'})
xlim([0 max(time_full_days)])
logged = log(perioddays_cwt);
inds = [49,98];
yticks([logged(inds)])
yticklabels({num2str(round(perioddays_cwt(inds)))})
colormap(parula);
colorbar()
set(gca,'FontSize',fsize)

% Delta selection figure 
subplot(4,6,[5,6,11,12])
colors = {'red','black','blue','cyan','green','magenta'};
for nn = 1: length(deltaselection.sampperdays)

    x = deltaselection.deltas;
    dat = deltaselection.MSE{nn}'; % Number of ‘Experiments’ In Data Set
    yMean = mean(dat,1);  % Mean Of All Experiments At Each Value Of ‘x’
    ySEM = std(dat,0,1)/sqrt(size(dat,1));   % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
    CI95 = tinv([0.025 0.975], (size(dat,1)-1)); % Calculate 95% Probability Intervals Of t-Distribution
    yCI95 = bsxfun(@times, ySEM, CI95(:)); % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’

    plot(x, yCI95+yMean, 'Color','blue') % Plot 95% Confidence Intervals Of All Experiments
    CIs = yCI95+yMean;
    curve1 = CIs(1,:);
    curve2 = CIs(2,:);
    xs = [x, fliplr(x)];
    inbetween = [curve1,fliplr(curve2)];
    fill(xs,inbetween, colors{nn},'FaceAlpha',0.5,'DisplayName','','EdgeColor','none','HandleVisibility','off')
    hold on
    labelname = strcat(num2str(deltaselection.sampperdays(nn)),{' '},'samples per day');
    plot(x, yMean, 'Color',colors{nn}, 'LineWidth',2,'DisplayName',labelname{:}) % Plot Mean Of All Experiments
    hold on
    grid
end 
xlabel('Delta')
ylabel('MSE')
title(strcat(subject,{' '},'75/25 cross validation'))
legend
set(gca,'FontSize',fsize)

% PART 2: BPDN SAMPLES AND MODEL OUTPUT

% original signal with samples
subplot(4,6,[13,14,15,16])
plot(time_full_days,signal, 'Color','black', 'DisplayName','raw data')
hold on
plot(t_days, measureddata_wds,'.','Color','black', 'DisplayName','sparse samples','MarkerSize',2.5)
xlabel('Time (days)')
ylabel('Spike rate')
title({'';''})
xlim([0 max(time_full_days)])
ax=gca; ax.YAxis.Exponent = 3;
colorbar()
ylim([0 max(signal) + 2000])
legend()
set(gca,'FontSize',fsize)

% recon overlayed on original signal (or recon alone)
subplot(4,6,[19,20,21,22])
plot(time_full_days,signal, 'Color','black','DisplayName','raw data')
hold on
plot(t_recon_days, reconsig_wds,'Color',[0.8500 0.3250 0.0980],'DisplayName','BPDN recon')
xlabel('Time (days)')
ylabel('Spike rate')
xlim([0 max(time_full_days)])
title({'';'Reconstructed signal, BPDN'})
ax=gca; ax.YAxis.Exponent = 3;
xlim([0 max(time_full_days)])
ylim([0 max(reconsig_wds) + 1000])
legend()
colorbar()
set(gca,'FontSize',fsize)

% overlay BPDN and CWT spectra
subplot(4,6,[17,18,23,24])

% spectrum for BPDN
yyaxis right

% no rescaling
plot(desiredperioddays, x.^2, 'Color',[0.8500 0.3250 0.0980], 'LineWidth',1);
hold on
plot(desiredperioddays(sighits),x(sighits).^2,'*','Color','black','MarkerSize',4)
xlabel('Period (days)')
ylabel('Power, BPDN')
set(gca,'FontSize',fsize)

% spectrum for cwt
yyaxis left
plot(perioddays_cwt, meanpowerCWTreg, 'Color',[0 0 0 0.3], 'LineWidth',1.5); % blue color with 0.6 alpha
hold on
plot(perioddays_cwt(locs),pks,'Color','cyan','o', 'MarkerSize',2)
xlabel('Period (days)')
ylabel({'';'';'Power,CWT'})
title('BPDN and CWT spectra')
ax=gca; ax.YAxis(1).Exponent = 5; ax.YAxis(1).Color = 'black';
xlim([-10 max(perioddays_cwt)])
set(gca,'FontSize',fsize)
