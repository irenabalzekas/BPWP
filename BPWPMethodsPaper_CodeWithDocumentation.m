%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Irena Balzekas

% Code to accompany manuscript titled:
% "Method for cycle detection in sparse, irregularly sampled, 
% long-term neuro-behavioral timeseries: Basis pursuit denoising with 
% polynomial detrending of long-term, inter-ictal epileptiform activity"

% This code is published under creative commons license. Please cite our
% publication accordingly. 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter sweeps

% SCRIPT NAME: BPDN_crossval_script_parallel_2D_5050
% WHAT SCRIPT DOES:
% - loads data, creates data struct that will be foundation moving forward
% - Patient-specific crossval to find ideal delta parameter for subject

% DEFINE FIXED PARAMS PER SUBJECT
N = 3000; 
monthsdata = 12; 
numrepeats = 10; % number of times to recalculate MSE
percent2drop = 10; % percent of samples excluded from initial fit
maxdegree = 3; % highest degree polynomial to fit long-term trend
maxperiod = 120; % (days) upper end of period range of interest - what frequencies want well represented in the dct basis
minperiod = 10; % (days) lower end of period range of interest
sampletype = 'densesubrange'; % parameter in dct basis design
subjects = {'M4','M5','M1'};
sampperdays = [1,3,5]; % average number of samples per day 
deltas = 0:50000:1000000; % range of error (measurement and model) delta values to test - parameter sweep

for p = 1:length(subjects)  
    
    subject = subjects{p};

    % LOAD DATA
%     filename = strcat('/mnt/eplab/Personal/Irena/RCS/fromdashboard/',subject,'/periodicity/',subject,'_periodicities_spikeseries_left_nozscore.csv');
    filename = strcat('P:\Personal\Irena\RCS\fromdashboard\',subject,'\periodicity\',subject,'_periodicities_spikeseries_left_nozscore.csv');

    data = readtable(filename);  
%     if strcmp(subject,'M5') % get rid of M5's early drop
%         data = data(82:end,:);
%     end
    signal = data.spikes_interp';
    time = (data.start_uutc./10^6)'; % raw data timestamps were in milliseconds
    % cut time down to specified duration
    indmax = round((monthsdata * 30 *24 *60 *60) / (time(2) - time(1))); 
    signal = signal(1:indmax);
    time = time(1:indmax);
    dt = round((max(time) - min(time))) / (N-1);
    t_full = 1:N;
    
    % Building the dct basis can be slow if it is large, good to do outside of loop
    [f_k, desiredperioddays] = frequency_sampling(N,maxperiod, minperiod, dt, sampletype);
    [dct_basis, scalefactor] = DCT2_basis(N, f_k);
    [poly_basis] = polynomial_basis(N, maxdegree);

    MSEpersamp = [];
    Xspersamp = [];
    zspersamp = [];
    measureddatapersamp = [];
    tpersamp = [];
    phipersamp = [];
    numsamps = [];
    training_md_persamp = [];
    training_t_persamp = [];
    training_phi_persamp = []; 
    testing_md_persamp = [];
    testing_t_persamp = [];
    testing_phi_persamp = [];
         
    for s = 1:length(sampperdays)
        
        % Define numsamp needed 
        numsamp = monthsdata * 30 * sampperdays(s);
        numsamps(s) = numsamp;

        [measureddata_init, t_init, phi_init, training_md, training_t, training_phi, testing_md, testing_t, testing_phi] = BPDN_samplerealdata_traintest(signal,time,numsamp, dt, percent2drop, numrepeats);
        % demo fig for visualizing/debugging
        % figure;plot(phi_init, measureddata_init);hold on;plot(training_phi(1,:), training_md(1,:),'*');hold on;plot(training_phi(5,:), training_md(5,:),'o')
        
        % Initialize variables 
        MSE = zeros(length(deltas), numrepeats);
        xs = zeros(length(deltas),numrepeats,N);
        zs = zeros(length(deltas),numrepeats, maxdegree+1);
        
        parfor dd = 1: length(deltas)

            delta = deltas(dd);

            for nn = 1: numrepeats

                % Training data
                measureddata = training_md(nn,:);
                t = training_t(nn,:);
                phi = training_phi(nn,:);

                % Testing data
                measured_missed = testing_md(nn,:);
                t_missed = testing_t(nn,:);
                phi_missed = testing_phi(nn,:);

                %to visualize for debugging
                %figure;plot(phi_init,measureddata_init,'o');hold on;plot(phi,
                %measureddata,'*');hold on;plot(phi_missed,measured_missed,'x') %
                
                % Run BPDN
                [x, z, ~, ~, ~, ~] = BPDN(delta, measureddata', phi, dct_basis, poly_basis);
                [reconsig] = BPDN_reconsig(f_k, x,scalefactor, maxdegree, z, t_full);
                xs(dd,nn,:) = x;
                zs(dd,nn,:) = z;
                
                % Get missed samples from reconstruction
                measured_missed_recon = reconsig(phi_missed); 

                % Calculate MSE
                MSE(dd,nn) = immse(measured_missed, measured_missed_recon); 
       
            end

        end 
        
        MSEpersamp{s} = MSE;
        Xspersamp{s} = xs;
        zspersamp{s} = zs;
        measureddatapersamp{s} = measureddata_init;
        tpersamp{s} = t_init;
        phipersamp{s} = phi_init;
        training_md_persamp{s} = training_md;
        training_t_persamp{s} = training_t;
        training_phi_persamp{s} = training_phi; 
        testing_md_persamp{s} = testing_md;
        testing_t_persamp{s} = testing_t;
        testing_phi_persamp{s} = testing_phi;
        
    end

    deltaselection.subject = subject;
    deltaselection.signal = signal;
    deltaselection.time = time;
    deltaselection.deltas = deltas;
    deltaselection.sampperdays = sampperdays;
    deltaselection.numsamps = numsamps;
    deltaselection.note = 'Repeated cross validation for each combo of delta and sample density (number of samples). Each cell corresponds to cross validation MSE vs delta for a particular density';
    deltaselection.MSE = MSEpersamp;
    deltaselection.N = N;
    deltaselection.dt = dt;
    deltaselection.measuredata_init = measureddatapersamp;
    deltaselection.t_init = tpersamp;
    deltaselection.phi_init = phipersamp;
    deltaselection.training_md_persamp = training_md_persamp;
    deltaselection.training_md_persamp= training_md_persamp;
    deltaselection.training_phi_persamp = training_phi_persamp; 
    deltaselection.testing_md_persamp = testing_md_persamp;
    deltaselection.testing_t_persamp = testing_t_persamp;
    deltaselection.testing_phi_persamp = testing_phi_persamp;
    deltaselection.dct_basis = dct_basis;
    deltaselection.maxdegree = maxdegree;
    deltaselection.minperiod = minperiod;
    deltaselection.maxperiod = maxperiod;
    deltaselection.sampletype = sampletype;
    deltaselection.poly_basis = poly_basis;
    deltaselection.f_k = f_k;
    deltaselection.desiredperioddays = desiredperioddays;
    deltaselection.x = Xspersamp;
    deltaselection.z = zspersamp;
    deltaselection.scalefactor = scalefactor;
    deltaselection.t_full = t_full;
    deltaselection.monthsdata = monthsdata;
    deltaselection.percent2drop = percent2drop;

%     savename = strcat('/mnt/eplab/Personal/Irena/BPDN_methods_paper/realdata/',subject,'/deltaselection/MSE_2D_DataDuration',num2str(monthsdata),'_PercentDrop',num2str(percent2drop),'_repeats',num2str(numrepeats),'_N',num2str(N),'_numsamp',num2str(numsamps),'stabletraintestdata.mat');
    savename = strcat('P:\Personal\Irena\BPDN_methods_paper\realdata\',subject,'\deltaselection\MSE_2D_DataDuration',num2str(monthsdata),'_PercentDrop',num2str(percent2drop),'_repeats',num2str(numrepeats),'_N',num2str(N),'_numsamp',num2str(numsamps),'stabletraintestdata.mat');

    save(savename,'deltaselection');
    
end 

%% Plot parameter sweeps visualize/ select optimal delta and sampling density 
% Use to identify delta value that minimizes MSE for patient and sample
% density

subjects = {'M1','M4','M5'};

colors = {'red','black','blue','cyan','green','magenta'};
fsize = 15;
fig = figure('Position',[10,10,2000,1000]);
for i = 1:length(subjects)
    
    subplot(1,length(subjects),i)
    subject = subjects{i};
   
    filename = strcat('P:\Personal\Irena\BPDN_methods_paper\realdata\',subject,'\deltaselection\MSE_2D_DataDuration12_PercentDrop25_repeats10_N3000_numsamp360  1080  1800stabletraintestdata.mat');
    load(filename);
    
    for nn = 1: length(deltaselection.sampperdays)
    
        x = deltaselection.deltas;
        dat = deltaselection.MSE{nn}'; % Number of ‘Experiments’ In Data Set
        yMean = mean(dat,1);  % Mean Of All Experiments At Each Value Of ‘x’
        ySEM = std(dat,0,1)/sqrt(size(dat,1));   % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
        CI95 = tinv([0.025 0.975], (size(dat,1)-1)); % Calculate 95% Probability Intervals Of t-Distribution
        yCI95 = bsxfun(@times, ySEM, CI95(:)); % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’

        %plot(x, yCI95+yMean, 'Color','blue') % Plot 95% Confidence Intervals Of All Experiments
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
    title(strcat(subject,{' '},'cross validation'))
    legend
    set(gca,'FontSize',fsize)
    
end 
savename = 'H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\deltaparametersweeps.tiff';
exportgraphics(fig,savename,'Resolution',600)

%% Save reviewer rated optimal delta and hard code save to data struct

subject = 'M5'; % Define patient specifically each time
locfolder = strcat('P:\Personal\Irena\BPDN_methods_paper\realdata\',subject,'\deltaselection\');
filename = strcat(locfolder,'MSE_2D_DataDuration12_PercentDrop25_repeats10_N3000_numsamp360  1080  1800stabletraintestdata_2.mat');
load(filename);

% % Adding the original data to the file post-hoc
% % filename = strcat('/mnt/eplab/Personal/Irena/RCS/fromdashboard/',subject,'/periodicity/',subject,'_periodicities_spikeseries_left_nozscore.csv');
% filename2 = strcat('P:\Personal\Irena\RCS\fromdashboard\',subject,'\periodicity\',subject,'_periodicities_spikeseries_left_nozscore.csv');
% monthsdata = 12;
% data = readtable(filename2);  
% signal = data.spikes_interp';
% time = (data.start_uutc./10^6)'; 
% % cut time down to specified duration
% indmax = round((monthsdata * 30 *24 *60 *60) / (time(2) - time(1)));
% signal = signal(1:indmax);
% time = time(1:indmax);

% deltaselection.signal = signal;
% deltaselection.time = time;
deltaselection.delta = []; % Minimum from MSE curve% min(deltaselection.MSE{nn}); % with nn specific to the sample density

save(filename,'deltaselection');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data drops - generating drops in data, try to detect spectral peaks 

% BPDN_datadrops_TestScript

% Percent hits vs datadrops
% Loads delta and other parameters from 2D parameter sweeps

% NOTE TP, FN etc. applies to whether or not it got all 3 component osc.
% If all 3, then true positive

% FIRST DEFINE SUBJECTS
subjects = {'M1','M4','M5'};

bandranges_days = [0.8,1.2; 13,22; 25,35; 80,120];
percentilethreshold = 99; 
amplitudepercentilethreshold = 90; 
numreps = 50;
samplesperday = 5;
gaplengths = [1,12,30,60];
continlengths = [2,24,90,120];

for i = 1%: length(subjects)
    
    subject = subjects{i};
    
    % Load reference data from parameter sweep (N, dt, data, etc.)
    locfolder = strcat('/mnt/eplab/Personal/Irena/BPDN_methods_paper/realdata/',subject,'/deltaselection/');
    destfolder = strcat('/mnt/eplab/Personal/Irena/BPDN_methods_paper/realdata/',subject,'/percentdetectionsVdatadrops/');
    load(strcat(locfolder,'MSE_2D_DataDuration12_PercentDrop25_repeats10_N3000_numsamp360  1080  1800stabletraintestdata.mat'));
    
    numsamps4sighits = samplesperday*(deltaselection.monthsdata*30);
    dt_ogdata = deltaselection.time(2) - deltaselection.time(1);

    % Re-sample 100X, t_hats, etc. and add to struct.
    measureddatas_100block = [];
    phis_100block = [];
    ts_100block = [];
    measureddatas_100block_cwt = [];
    ts_100block_cwt = [];
    period_days_cwt = [];
    meanpowerreg_cwt = [];
    peakvalues_cwt = [];
    peaklocs_cwt = [];
    peakperiods_cwt = [];
    
    for nn = 1:length(gaplengths)
        mds = [];
        ps = [];
        ts = [];
        mds_cwt = [];
        ts_cwt = [];
        period_days_cwt_temp = [];
        meanpowerreg_cwt_temp = [];
        peakvalues_cwt_temp = [];
        peaklocs_cwt_temp = [];
        peakperiods_cwt_temp = [];
        
        for ii = 1:100
            numsamp = numsamps4sighits;

            % SETUP FOR BPDN AND INTERPOLATED CWT
            [md, t, ~,~, ~,~,md_cwt, t_cwt] = BPDN_create_datadrops(deltaselection.signal, deltaselection.time, gaplengths(nn), continlengths(nn), 'both', samplesperday, dt_ogdata);
            phi = floor((t-min(deltaselection.time))/deltaselection.dt) + 1; % indices for when each time data were sampled
            mds{ii} = md; % measureddatas
            ps{ii} = phi; % phi's 
            ts{ii} = t;
            mds_cwt{ii} = md_cwt;
            ts_cwt{ii} = t_cwt;
            
            % Calculate CWT outputs for each block that was created and
            % interpolated 
            [perioddays_cwt,meanpowerCWTreg, peakvalues, peaklocs, peakperiods] = BPDN_cwt2findcycles_matlab2021(md_cwt, t_cwt(2)-t(1)); % redundant, but just for peak finding
            period_days_cwt_temp{ii} = perioddays_cwt;
            meanpowerreg_cwt_temp{ii} = meanpowerCWTreg;
            peakvalues_cwt_temp{ii} = peakvalues;
            peaklocs_cwt_temp{ii} = peaklocs;
            peakperiods_cwt_temp{ii} = peakperiods;
            
        end 
        measureddatas_100block{nn} = mds;
        phis_100block{nn} = ps;
        ts_100block{nn} = ts;
        measureddatas_100block_cwt{nn} = mds_cwt;
        ts_100block_cwt{nn} = ts_cwt;
        period_days_cwt{nn} = period_days_cwt_temp;
        meanpowerreg_cwt{nn} = meanpowerreg_cwt_temp;
        peakvalues_cwt{nn} = peakvalues_cwt_temp;
        peaklocs_cwt{nn} = peaklocs_cwt_temp;
        peakperiods_cwt{nn} = peakperiods_cwt_temp;
    end
    deltaselection.measureddatas_100block = measureddatas_100block;
    deltaselection.phis_100block = phis_100block;
    deltaselection.ts_100block = ts_100block;

    % Find "true" peaks based on cwt. Add to structure too. 
    [perioddays_cwt,meanpowerCWTreg, peakvalues, peaklocs, peakperiods] = BPDN_cwt2findcycles_matlab2021(deltaselection.signal, deltaselection.time(2) - deltaselection.time(1)); % redundant, but just for peak finding
    deltaselection.cwt.period_days = perioddays_cwt;
    deltaselection.cwt.meanpowerreg = meanpowerCWTreg;
    deltaselection.cwt.peakvalues = peakvalues;
    deltaselection.cwt.peaklocs = peaklocs;
    deltaselection.cwt.peakperiods = peakperiods;
    
    % ALSO OUTPUT A CWT PER INTERPOLATION.
    deltaselection.cwtinterp.measureddatas_100block_cwt = measureddatas_100block_cwt;
    deltaselection.cwtinterp.ts_100block_cwt = ts_100block_cwt;
    deltaselection.cwtinterp.period_days = period_days_cwt;
    deltaselection.cwtinterp.meanpowerreg = meanpowerreg_cwt;
    deltaselection.cwtinterp.peakvalues = peakvalues_cwt;
    deltaselection.cwtinterp.peaklocs = peaklocs_cwt;
    deltaselection.cwtinterp.peakperiods = peakperiods_cwt;
    
    percenthits_persamplings = [];
    
    for nn = 1:length(gaplengths) 
        
        % LOOP through the number of y and t_hat resamplings in simsig
        solved = [];
        foundosc_offsets_repeats = [];
        foundosc_offsets_repeats_amp = [];
        x_sighits_mat = [];
        x_reshufflings = [];
        x_percentiles = [];
        z_sighits_mat = [];
        z_reshufflings = [];
        z_percentiles = [];
        TPs = [];
        FPs = [];
        allpeaksdetected = [];

        for rr = 1: numreps%size(deltaselection.measureddatas_100block{nn},2) 

            % DEFINE BASES, PHI, ETC.
            %dt = ((simsig.t_hat(end,rr) - simsig.t_hat(1,rr)) * sparsity)/ length(simsig.y(:,rr));
            [x, xreshuffled, xpercents, z, zreshuffled, zpercents, reshuff_optval] = BPDN_wReshuffling_short(deltaselection.delta, deltaselection.measureddatas_100block{nn}{rr}', deltaselection.phis_100block{nn}{rr}, deltaselection.dct_basis, deltaselection.poly_basis);
%             [x, z, A, B, cvx_status, cvx_optval] = BPDN(deltaselection.delta, deltaselection.measureddatas_100block{nn}(:,rr), deltaselection.phis_100block{nn}(:,rr), deltaselection.dct_basis, deltaselection.poly_basis);% quick debug
            x_reshufflings(:,:,rr) = xreshuffled; % will be 3d matrix
            z_reshufflings(:,:,rr) = zreshuffled;
            x_percentiles(:,rr) = xpercents;
            z_percentiles(:,rr) = zpercents;
            solved{rr} = reshuff_optval;

            % DETERMINE PERCENT SIGNIFICANCE BASED ON PERCENTILE
            % THRESHOLD-------------------------------------------------------
            [x_sighits] = percentilehits(xpercents,percentilethreshold);
            x_sighits_mat(:,rr) = x_sighits; 
            [z_sighits] = percentilehits(zpercents,percentilethreshold);
            z_sighits_mat(:,rr) = z_sighits; 

            % HOW FAR WERE SIG HITS (PEAKS) FROM KNOWN PERIODs (closest found peak)
            realosc = peakperiods;
            [foundosc_ind] = find(x_sighits ==1);
            foundosc = deltaselection.desiredperioddays(foundosc_ind);
            oscdifferences = []; % find distance of real osc from identified osc per real osc
            for p = 1: length(realosc)
                try
                    oscdifferences(p) = min(abs(foundosc - realosc(p)));
                catch 
                    oscdifferences(p) = NaN;
                end
            end
            foundosc_offsets_repeats(:,rr) = oscdifferences; 

            % "TRUE HITS" where offset found osc is less than 2 days
            temp = find(oscdifferences < 2);
            truepos = zeros(1,length(peakperiods));
            truepos(temp) = 1;
            TPs(:,rr) = truepos;
            FPs(:,rr) = length(foundosc) - sum(truepos);
            if sum(truepos) == length(realosc)
                allpeaksdetected(:,rr) = 1;
            elseif sum(truepos) ~= length(realosc)
                allpeaksdetected(:,rr) = 0;
            end 
              
            % WERE THERE SIG HITS IN THE MEANINGFUL BAND RANGES? % AND were
            % they high amplitude? ** should turn this into a function
            pow = x.^2;
            ampthresh = prctile(pow, amplitudepercentilethreshold); % set amplitude percentile threshold 
            foundoscpow = pow(foundosc_ind);
            foundoscperband = [];
            for b = 1:length(bandranges_days)
                peaksinband = [];
                for os = 1: length(foundosc)
                    %found osc value in band AND associated power above X percentile 
                    if (foundosc(os) <= bandranges_days(b,2)) && (foundosc(os) >= bandranges_days(b,1)) && (foundoscpow(os) >= ampthresh)
                        peaksinband(os) = 1;
                    else
                        peaksinband(os) = 0;
                    end 
                end
                foundoscperband(b) = sum(peaksinband);
            end 
            
        end 
    
        % little figure for debugging/visualizing
%         figure;plot(deltaselection.desiredperioddays,x.^2);hold on;plot(deltaselection.desiredperioddays(find(x_sighits ==1)),x(find(x_sighits ==1)).^2,'*')
        
        % Calculate percent significant hits per frequency -> with basic percentile
        xnumhits = sum(x_sighits_mat,2); 
        xpercenthits = (xnumhits./ size(x_sighits_mat,2)) * 100; 
        znumhits = sum(z_sighits_mat,2); 
        zpercenthits = (znumhits./ size(z_sighits_mat,1)) * 100; 

        % Update simsig with results 
        percenthits.samplesperday = samplesperday;
        percenthits.numsamps4sighits = numsamps4sighits;
        percenthits.gaplengths = gaplengths;
        percenthits.continlengths = continlengths;
        percenthits.numreps = numreps;
        percenthits.solverstatus = solved;
        percenthits.foundosc_offsets_reshufflings = foundosc_offsets_repeats; % closest found peaks to the real/known oscillations 
        percenthits.foundosc = foundosc;
        percenthits.foundoscpow = foundoscpow;
        percenthits.foundosc_ind = foundosc_ind;

        percenthits.x = x;
        percenthits.xpercenthits = xpercenthits;
        percenthits.x_sighits = x_sighits_mat;
        percenthits.x_reshufflings = x_reshufflings;
        percenthits.x_percentiles = x_percentiles;

        percenthits.z = z;
        percenthits.zpercenthits = zpercenthits;
        percenthits.z_sighits = z_sighits_mat;
        percenthits.z_reshufflings = z_reshufflings;
        percenthits.z_percentiles = z_percentiles;
        
        percenthits.TP = TPs;
        percenthits.FP = FPs;
        percenthits.allpeaksdetected = allpeaksdetected;
        
        percenthits.bandranges_days = bandranges_days;
        percenthits.foundoscperband = foundoscperband; % this had percentile threshold too
        percenthits.powerpercentilethreshold = amplitudepercentilethreshold;
        
        percenthits_persampling{nn} = percenthits;
       
    end
    
    deltaselection.percenthits = percenthits_persampling;
    
    loopfilename = strcat(destfolder,'MSE_2D_DataDuration12_PercentDrop25_repeats10_N3000_numsamp360  1080  1800stabletraintestdata','_withpercenthits_datadropstestrun','.mat');
    save(loopfilename,'deltaselection'); 
    display(strcat('Finished sig-hits run for subject',subject))
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIMULATIONS

%% Generate simulated data + run delta parameter sweep: One year spike rate timeseries 

% DEFINE BASIC PARAMS
subject = 'S1';
years = 1 ;
n = 3*24*365*years; % 3 samples per hour for one year 
totalsec = years*60*60*365*24;
dt = round( totalsec / n); % year in seconds/n - % confirm this is indeed 20 min***
% manually define desired cycles, their amplitude, and noise (cycles chosen here align with well conserved spike cycles in literature)
poi = [1, 7,15, 21, 30, 50, 100]; % periods (in days) of interest to include in simulated signal
maxsigamp = [60,20,10,10,40,25,20]; % amplitude associated with each cycle in poi
polycoef = [-.005 1*10^-50 10^-100  0 0]; % polynomial coefficients for first, second, third, fourth, fifth degree 
noiseparams = [16,10,10,10,6,4,4]; % for each cycle in poi

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

figure
plot(simsig,'Color','black')

simdata.spikes = simsig;
simdata.time = t;
simdata.periods_days = poi;
simdata.polycoef = polycoef;

% check CWT output to eval aligment with actual spike data spectra
[meanpowerCWTL2norm,perioddays_cwt,WT, fwt] = wavelet_decomp_L2norm(simsig, 1/dt);
meanpowerCWTreg = nanmean(abs(WT).^2,2);
reconsigCWT = icwt(WT);
figure;plot(perioddays_cwt, meanpowerCWTreg)

save('P:\Personal\Irena\BPDN_methods_paper\simulateddata\S1\rawdata\S1_simdata.mat', 'simdata')

% RUN PAREMETER SWEEP - DEFINE FIXED PARAMS PER SUBJECT
N = 3000; 
monthsdata = 12;
numrepeats = 10; % number of times to recalculate MSE
percent2drop = 25;
maxdegree = 3;
maxperiod = 120;
minperiod = 10; 
sampletype = 'densesubrange';
subjects = {'S1'};
sampperdays = [1,3,5];
deltas = 0:50000:1000000;

for p = 1:length(subjects)  
    
    sub = subjects{p};

    % LOAD DATA
    filename = strcat('/mnt/eplab/Personal/Irena/BPDN_methods_paper/simulateddata/',sub,'/rawdata/',sub,'_simdata.mat');

    load(filename);  

    signal = simdata.spikes';
    time = simdata.time'; 
    
    % cut time down to specified duration
    indmax = round((monthsdata * 30 *24 *60 *60) / (time(2) - time(1)));
    signal = signal(1:indmax);
    time = time(1:indmax);
    dt = round((max(time) - min(time))) / (N-1);
    t_full = 1:N;
    
    % Building the dct basis can be slow if its large, good to do outside loop
    [f_k, desiredperioddays] = frequency_sampling(N,maxperiod, minperiod, dt, sampletype);
    [dct_basis, scalefactor] = DCT2_basis(N, f_k); % making the basis takes a long time. want to pre-define it before running BPDN for loops and loops 
    [poly_basis] = polynomial_basis(N, maxdegree);

    MSEpersamp = [];
    Xspersamp = [];
    zspersamp = [];
    measureddatapersamp = [];
    tpersamp = [];
    phipersamp = [];
    numsamps = [];
    training_md_persamp = [];
    training_t_persamp = [];
    training_phi_persamp = []; 
    testing_md_persamp = [];
    testing_t_persamp = [];
    testing_phi_persamp = [];
    
    for s = 1:length(sampperdays)
        
        % Define numsamp needed 
        numsamp = monthsdata * 30 * sampperdays(s);
        numsamps(s) = numsamp;

        [measureddata_init, t_init, phi_init, training_md, training_t, training_phi, testing_md, testing_t, testing_phi] = BPDN_samplerealdata_traintest(signal,time,numsamp, dt, percent2drop, numrepeats);
        % debug fig
        % figure;plot(phi_init, measureddata_init);hold on;plot(training_phi(1,:), training_md(1,:),'*');hold on;plot(training_phi(5,:), training_md(5,:),'o')
        
        % Initialize variables 
        MSE = zeros(length(deltas), numrepeats);
        xs = zeros(length(deltas),numrepeats,N);
        zs = zeros(length(deltas),numrepeats, maxdegree+1);
        
        parfor dd = 1: length(deltas)

            delta = deltas(dd);

            for nn = 1: numrepeats

                % Training data
                measureddata = training_md(nn,:);
                t = training_t(nn,:);
                phi = training_phi(nn,:);

                % Testing data
                measured_missed = testing_md(nn,:);
                t_missed = testing_t(nn,:);
                phi_missed = testing_phi(nn,:);

                %figure;plot(phi_init,measureddata_init,'o');hold on;plot(phi,
                %measureddata,'*');hold on;plot(phi_missed,measured_missed,'x') %
                %to visualize for debugging

                % Run BPDN
                [x, z, ~, ~, ~, ~] = BPDN(delta, measureddata', phi, dct_basis, poly_basis);
                [reconsig] = BPDN_reconsig(f_k, x,scalefactor, maxdegree, z, t_full);
                xs(dd,nn,:) = x;
                zs(dd,nn,:) = z;

                % Get missed samples from reconstruction
                measured_missed_recon = reconsig(phi_missed); 

                % Calculate MSE
                MSE(dd,nn) = immse(measured_missed, measured_missed_recon); 
       
            end

        end 
        
        MSEpersamp{s} = MSE;
        Xspersamp{s} = xs;
        zspersamp{s} = zs;
        measureddatapersamp{s} = measureddata_init;
        tpersamp{s} = t_init;
        phipersamp{s} = phi_init;
        training_md_persamp{s} = training_md;
        training_t_persamp{s} = training_t;
        training_phi_persamp{s} = training_phi; 
        testing_md_persamp{s} = testing_md;
        testing_t_persamp{s} = testing_t;
        testing_phi_persamp{s} = testing_phi;
        
    end

    deltaselection.subject = sub;
    deltaselection.signal = signal;
    deltaselection.time = time;
    deltaselection.deltas = deltas;
    deltaselection.sampperdays = sampperdays;
    deltaselection.numsamps = numsamps;
    deltaselection.note = 'Repeated cross validation for each combo of delta and sample density (number of samples). Each cell corresponds to cross validation MSE vs delta for a particular density';
    deltaselection.MSE = MSEpersamp;
    deltaselection.N = N;
    deltaselection.dt = dt;
    deltaselection.measuredata_init = measureddatapersamp;
    deltaselection.t_init = tpersamp;
    deltaselection.phi_init = phipersamp;
    deltaselection.training_md_persamp = training_md_persamp;
    deltaselection.training_md_persamp= training_md_persamp;
    deltaselection.training_phi_persamp = training_phi_persamp; 
    deltaselection.testing_md_persamp = testing_md_persamp;
    deltaselection.testing_t_persamp = testing_t_persamp;
    deltaselection.testing_phi_persamp = testing_phi_persamp;
    deltaselection.dct_basis = dct_basis;
    deltaselection.maxdegree = maxdegree;
    deltaselection.minperiod = minperiod;
    deltaselection.maxperiod = maxperiod;
    deltaselection.sampletype = sampletype;
    deltaselection.poly_basis = poly_basis;
    deltaselection.f_k = f_k;
    deltaselection.desiredperioddays = desiredperioddays;
    deltaselection.x = Xspersamp;
    deltaselection.z = zspersamp;
    deltaselection.scalefactor = scalefactor;
    deltaselection.t_full = t_full;
    deltaselection.monthsdata = monthsdata;
    deltaselection.percent2drop = percent2drop;

    savename = strcat('/mnt/eplab/Personal/Irena/BPDN_methods_paper/simulateddata/',sub,'/deltaselection/MSE_2D_DataDuration',num2str(monthsdata),'_PercentDrop',num2str(percent2drop),'_repeats',num2str(numrepeats),'_N',num2str(N),'_numsamp',num2str(numsamps),'stabletraintestdata.mat');
    save(savename,'deltaselection');
    
end 

%% plot simulations delta param sweep

subjects = {'S1'};

colors = {'red','black','blue','cyan','green','magenta'};
fsize = 15;
fig = figure('Position',[10,10,700,1000]);
for i = 1:length(subjects)
    
    subplot(1,length(subjects),i)
    subject = subjects{i};
   
    filename =  strcat('P:\Personal\Irena\BPDN_methods_paper\simulateddata\',subject,'\deltaselection\MSE_2D_DataDuration12_PercentDrop25_repeats10_N3000_numsamp360  1080  1800_stabletraintestdata.mat');
    load(filename);
    
    for nn = 1: length(deltaselection.sampperdays)
    
        x = deltaselection.deltas;
        dat = deltaselection.MSE{nn}'; % Number of ‘Experiments’ In Data Set
        yMean = mean(dat,1);  % Mean Of All Experiments At Each Value Of ‘x’
        ySEM = std(dat,0,1)/sqrt(size(dat,1));   % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
        CI95 = tinv([0.025 0.975], (size(dat,1)-1)); % Calculate 95% Probability Intervals Of t-Distribution
        yCI95 = bsxfun(@times, ySEM, CI95(:)); % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’

        %plot(x, yCI95+yMean, 'Color','blue') % Plot 95% Confidence Intervals Of All Experiments
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
end 

savename = 'H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\deltaparametersweeps_simdata.tiff';
exportgraphics(fig,savename,'Resolution',600)

%% update file with preferred delta 

subject = 'S1';
locfolder = strcat('P:\Personal\Irena\BPDN_methods_paper\simulateddata\',subject,'\deltaselection\');
    filename =  strcat('P:\Personal\Irena\BPDN_methods_paper\simulateddata\',subject,'\deltaselection\MSE_2D_DataDuration12_PercentDrop25_repeats10_N3000_numsamp360  1080  1800_stabletraintestdata.mat');
load(filename);

deltaselection.delta = []; 

save(filename,'deltaselection');

%% Run significance testing/repeats AND update simulation with cwt struct - reshuffling 100x

subjects = {'S1'};
sampperdaydesired = 5;

for i = 1: length(subjects)
    
    sub = subjects{i};
    % Load reference data from parameter sweep (N, dt, data, etc.)
    destfolder = strcat('P:\Personal\Irena\BPDN_methods_paper\simulateddata\',subject,'\modeloutput\');
    filename =  strcat('P:\Personal\Irena\BPDN_methods_paper\simulateddata\',subject,'\deltaselection\MSE_2D_DataDuration12_PercentDrop25_repeats10_N3000_numsamp360  1080  1800stabletraintestdata.mat');
    load(filename);

    % Define basics
    time = deltaselection.time; 
    dt_ogdata = time(2)-time(1);
    % Select data based on desired num per day
    ind = find(deltaselection.sampperdays == sampperdaydesired);

    measureddata_wds = deltaselection.measuredata_init{ind}(:);
    time_wds =  deltaselection.t_init{ind}(:);
    phi = round((time_wds-min(time))/deltaselection.dt) + 1; % indices for when each time data were sampled
    numsamp = deltaselection.numsamps(ind);

    % First, BPDN
    % [phi_wds, N, t_full, f_k, desiredperioddays, dct_basis, scalefactor, poly_basis] = BPDN_fullSetup(time_wds,dt,maxperiod,minperiod,sampletype,maxdegree);
    % [x_wds, z_wds, A, B, cvx_status, cvx_optval] = BPDN(delta, measureddata_wds', phi_wds, dct_basis, poly_basis);
    [x_wds, xreshuffled, xpercentiles, z_wds, zreshuffled, zpercentiles, reshuff_optval] = BPDN_wReshuffling_short(deltaselection.delta, measureddata_wds, phi, deltaselection.dct_basis, deltaselection.poly_basis);
    sighits = find(xpercentiles > 99);
    [reconsig_wds] = BPDN_reconsig(deltaselection.f_k, x_wds, deltaselection.scalefactor, deltaselection.maxdegree, z_wds, deltaselection.t_full);

    % Rescale power by frequency (% Change from l1 to l2 norm for cwt.  l2 norm reduces amplitude of high freq data, while l1 does not. (scale is inversely related to frequency))
    x2 = [];
    for i=1:length(x_wds) 
        x2(i) = x_wds(i).*1/sqrt(deltaselection.f_k(i));  
    end
    x_power_rescale = abs(x2);

    % Nick CWT - with L2 adjustment
    [meanpowerCWTL2norm, perioddays_cwt,WT, fwt] = wavelet_decomp_L2norm(deltaselection.signal, 1/dt_ogdata);
    [perioddays_cwt,meanpowerCWTreg,pks,locs,reconsigCWT] = BPDN_cwt2findcycles(deltaselection.signal, dt_ogdata); % redundant, but just for peak finding

    % Modify time scale for plots
    time_full_days = (time - min(time))./(60*60*24);
    t_days = (time_wds - min(time_wds))./ (60*60*24);
    t_recon_days = (deltaselection.t_full*deltaselection.dt) ./ (60*60*24);

    % Saving data for fig
    deltaselection.fig.sampperday = sampperdaydesired;
    deltaselection.fig.x = x_wds;
    deltaselection.fig.xreshuffled = xreshuffled;
    deltaselection.fig.xpercentiles = xpercentiles;
    deltaselection.fig.z = z_wds;
    deltaselection.fig.zpercentiles = zpercentiles;
    deltaselection.fig.reshuff_optval = reshuff_optval;
    deltaselection.fig.measureddata = measureddata_wds;
    deltaselection.fig.time = time_wds;
    deltaselection.fig.sighits = sighits;
    deltaselection.fig.reconsig = reconsig_wds;
    deltaselection.fig.x_power_rescale = x_power_rescale;
    deltaselection.fig.meanpowerCWTreg = meanpowerCWTreg;
    deltaselection.fig.perioddays_cwt = perioddays_cwt;
    deltaselection.fig.WT = WT;
    deltaselection.fig.fwt = fwt;
    deltaselection.fig.reconsigCWT = reconsigCWT;
    deltaselection.fig.pks = pks;
    deltaselection.fig.locs = locs;
    deltaselection.fig.time_full_days = time_full_days;
    deltaselection.fig.t_days = t_days;
    deltaselection.fig.t_recon_days = t_recon_days;

    % Append and save
    savename = strcat(destfolder,'Outputfor',num2str(deltaselection.duration),'_N',num2str(deltaselection.N),'_numsamp',num2str(numsamp),'_traintest_forfig.mat');
    save(savename,'deltaselection');
end 

%% Run significance testing/repeats AND update simulation with cwt struct - reshuffling 1,000x - per reviewer

subjects = {'M1'};
sampperdaydesired = 5;

for i = 1: length(subjects)
    
    tic
    subject = subjects{i};
    % Load reference data from parameter sweep (N, dt, data, etc.)
%     destfolder = strcat('P:\Personal\Irena\BPDN_methods_paper\simulateddata\',subject,'\modeloutput\');
%     filename =  strcat('P:\Personal\Irena\BPDN_methods_paper\simulateddata\',subject,'\deltaselection\MSE_2D_DataDuration12_PercentDrop25_repeats10_N3000_numsamp360  1080  1800_stabletraintestdata.mat');
    %destfolder = strcat('/mnt/eplab/Personal/Irena/BPDN_methods_paper/simulateddata/',sub,'/modeloutput/');
    %filename =  strcat('/mnt/eplab/Personal/Irena/BPDN_methods_paper/simulateddata/',sub,'/deltaselection/MSE_2D_DataDuration12_PercentDrop25_repeats10_N3000_numsamp360  1080  1800_stabletraintestdata.mat');
    
    destfolder = strcat('/mnt/eplab/Personal/Irena/BPDN_methods_paper/realdata/',subject,'/modeloutput/');
    filename =  strcat('/mnt/eplab/Personal/Irena/BPDN_methods_paper/realdata/',subject,'/deltaselection/MSE_2D_DataDuration12_PercentDrop25_repeats10_N3000_numsamp360  1080  1800stabletraintestdata.mat');
%     
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
    
    [x_wds, xreshuffled, xpercentiles, z_wds, zreshuffled, zpercentiles, reshuff_optval] = BPDN_wReshuffling_long(delta, measureddata_wds, phi, dct_basis, poly_basis);

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
    
    elapsedtime = toc
    % Saving data for fig
    deltaselection.fig.sampperday = sampperdaydesired;
    deltaselection.fig.x = x_wds;
    deltaselection.fig.reshufflingruntime = elapsedtime;
    deltaselection.fig.xreshuffled = xreshuffled;
    deltaselection.fig.xpercentiles = xpercentiles;
    deltaselection.fig.z = z_wds;
    deltaselection.fig.zpercentiles = zpercentiles;
    deltaselection.fig.reshuff_optval = [];
    deltaselection.fig.measureddata = measureddata_wds;
    deltaselection.fig.time = time_wds;
    deltaselection.fig.sighits = sighits;
    deltaselection.fig.reconsig = reconsig_wds;
    deltaselection.fig.x_power_rescale = x_power_rescale;
    deltaselection.fig.meanpowerCWTreg = meanpowerCWTreg;
    deltaselection.fig.perioddays_cwt = perioddays_cwt;
    deltaselection.fig.WT = WT;
    deltaselection.fig.fwt = fwt;
    deltaselection.fig.reconsigCWT = reconsigCWT;
    deltaselection.fig.pks = pks;
    deltaselection.fig.locs = locs;
    deltaselection.fig.time_full_days = time_full_days;
    deltaselection.fig.t_days = t_days;
    deltaselection.fig.t_recon_days = t_recon_days;

    % Append and save
    savename = strcat(destfolder,'Outputfor',num2str(duration),'_N',num2str(N),'_numsamp',num2str(numsamp),'_traintest_forfig_1000x_test.mat');
    save(savename,'deltaselection');
    
end 

%% Run significance testing/repeats AND update simulation with cwt struct - reshuffling 10,000x 

% NOTE LINE 1008 starting index should be 1

subjects = {'S1'};
sampperdaydesired = 5;

for i = 1: length(subjects)
    
    tic
    subject = subjects{i};
    % Load reference data from parameter sweep (N, dt, data, etc.)
%     destfolder = strcat('P:\Personal\Irena\BPDN_methods_paper\simulateddata\',subject,'\modeloutput\');
%     filename =  strcat('P:\Personal\Irena\BPDN_methods_paper\simulateddata\',subject,'\deltaselection\MSE_2D_DataDuration12_PercentDrop25_repeats10_N3000_numsamp360  1080  1800_stabletraintestdata.mat');
    %destfolder = strcat('/mnt/eplab/Personal/Irena/BPDN_methods_paper/simulateddata/',sub,'/modeloutput/');
    %filename =  strcat('/mnt/eplab/Personal/Irena/BPDN_methods_paper/simulateddata/',sub,'/deltaselection/MSE_2D_DataDuration12_PercentDrop25_repeats10_N3000_numsamp360  1080  1800_stabletraintestdata.mat');
    
    % if running on lambda and saving to hydrogen
    destfolder = strcat('/mnt/Hydrogen/Irena/BPWP/realdata/',subject,'/modeloutput/');
    if strcmp(subject,'S1') == 1
        filename =  strcat('/mnt/Hydrogen/Irena/BPWP/realdata/',subject,'/deltaselection/MSE_2D_DataDuration12_PercentDrop25_repeats10_N3000_numsamp360  1080  1800_stabletraintestdata.mat');
    else 
        filename =  strcat('/mnt/Hydrogen/Irena/BPWP/realdata/',subject,'/deltaselection/MSE_2D_DataDuration12_PercentDrop25_repeats10_N3000_numsamp360  1080  1800stabletraintestdata.mat');
    end  
    
    % if running via eplab
    %destfolder = strcat('/mnt/eplab/Personal/Irena/BPDN_methods_paper/realdata/',subject,'/modeloutput/');
    %filename =  strcat('/mnt/eplab/Personal/Irena/BPDN_methods_paper/realdata/',subject,'/deltaselection/MSE_2D_DataDuration12_PercentDrop25_repeats10_N3000_numsamp360  1080  1800stabletraintestdata.mat');

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
    
    % Running 10,000x reshuffling in parallel with BPDN_wReshuffling_long
    % causes memory issues. Therefore running it BPDN_wReshuffling_long
    % (1,000x) for 10 times then rebuilding tables with 10,000x format
    numit = 10; % number of separately saved batches 
    for rr = 8:numit  % CHANGE BACK
        [x_wds_temp, xreshuffled_temp, ~, z_wds_temp, zreshuffled_temp, ~, ~] = BPDN_wReshuffling_long(delta, measureddata_wds, phi, dct_basis, poly_basis);
        % save each separate one separately
        temp.x_wds = x_wds_temp;
        temp.xreshuffled = xreshuffled_temp;
        temp.z_wds = z_wds_temp;
        temp.zreshuffled = zreshuffled_temp;
        temp.deltaselection = deltaselection;
        savename = strcat(destfolder,'Outputfor',num2str(duration),'_N',num2str(N),'_numsamp',num2str(numsamp),'_traintest_forfig_1000x_',num2str(rr),'of',num2str(numit),'_for10000.mat');
        save(savename,'temp');
        clear temp
    end 
    %load individuals and merge into one
    x_wds = [];
    xreshuffled = [];
    z_wds = [];
    zreshuffled = [];
    for rr = 1:numit 
        load(strcat(destfolder,'Outputfor',num2str(duration),'_N',num2str(N),'_numsamp',num2str(numsamp),'_traintest_forfig_1000x_',num2str(rr),'of',num2str(numit),'_for10000.mat'))
        x_wds = vertcat(x_wds, temp.x_wds);
        xreshuffled = vertcat(xreshuffled,temp.xreshuffled);
        z_wds = horzcat(z_wds,temp.z_wds);
        zreshuffled = vertcat(zreshuffled, temp.zreshuffled);
        clear temp
    end 
    % calculate percentiles based on the whole of the outputs (instead of
    % within function)
    % get percentile status for each output in x from original data
    [xpercentiles] = BPDN_getPercentile(xreshuffled,x_wds(:,1)); % just taking the first x
    [zpercentiles] = BPDN_getPercentile(zreshuffled,z_wds(:,1)); % just taking the first z
 
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
    
    elapsedtime = toc
    % Saving data for fig
    deltaselection.fig.sampperday = sampperdaydesired;
    deltaselection.fig.x = x_wds;
    deltaselection.fig.reshufflingruntime = elapsedtime;
    deltaselection.fig.xreshuffled = xreshuffled;
    deltaselection.fig.xpercentiles = xpercentiles;
    deltaselection.fig.z = z_wds;
    deltaselection.fig.zpercentiles = zpercentiles;
    deltaselection.fig.reshuff_optval = [];
    deltaselection.fig.measureddata = measureddata_wds;
    deltaselection.fig.time = time_wds;
    deltaselection.fig.sighits = sighits;
    deltaselection.fig.reconsig = reconsig_wds;
    deltaselection.fig.x_power_rescale = x_power_rescale;
    deltaselection.fig.meanpowerCWTreg = meanpowerCWTreg;
    deltaselection.fig.perioddays_cwt = perioddays_cwt;
    deltaselection.fig.WT = WT;
    deltaselection.fig.fwt = fwt;
    deltaselection.fig.reconsigCWT = reconsigCWT;
    deltaselection.fig.pks = pks;
    deltaselection.fig.locs = locs;
    deltaselection.fig.time_full_days = time_full_days;
    deltaselection.fig.t_days = t_days;
    deltaselection.fig.t_recon_days = t_recon_days;

    % Append and save
    savename = strcat(destfolder,'Outputfor',num2str(duration),'_N',num2str(N),'_numsamp',num2str(numsamp),'_traintest_forfig_10000x.mat');
    save(savename,'deltaselection');
    
end 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 1

% First part - Code N/A: Schematic diagram and some raw data

subject ='M5';

% Load reference data from parameter sweep (N, dt, data, etc.) + sighits
destfolder = strcat('P:\Personal\Irena\MoodCyclesSpikeCycles_Paper\realdata\',subject,'\modeloutput\');
load(strcat(destfolder,'Outputfor_MSE_percent2drop10_N1500_traintest_forfig.mat'));

% Load spike data 
spikefile = strcat('P:\Personal\Irena\MoodCyclesSpikeCycles_Paper\realdata\',subject,'\spikes\spikeinfo.mat');
load(spikefile)

% Centralize time
tbase = spikeinfo.time(1);
spiketimedays = (spikeinfo.time - tbase)./(60*60*24);
spikes = spikeinfo.signal;
dates = datetime(spikeinfo.time, 'ConvertFrom','posixtime','TimeZone','local');

moodtimedays = (deltaselection.time - tbase)./(60*60*24);
moodrecontimedays = (deltaselection.t_full.*deltaselection.dt) ./ (60*60*24);
sighits = deltaselection.fig.sighits;

% couple days of spikes
indrange = [2000:2300];
timeshort = spiketimedays(indrange) - spiketimedays(min(indrange));
datesshort = dates(indrange);
hoursshort = hour(datesshort);
% get start/ends of nighttime windows - will be using to make nighttime
% gray
starts = find(hoursshort == 22);
tocut = find(diff(starts)==1);
starts(tocut) = [];
ends = find(hoursshort == 8);
tocut = find(diff(ends)==1);
ends(tocut) = [];
spikeshort = spikes(indrange);
tocut = find(ends < starts(1));
ends(tocut) = [];

% couple months of spikes
indrange = [500:6400];
timemed = spiketimedays(indrange) - spiketimedays(min(indrange));
spikesmed = spikes(indrange);

fsize = 6;

% Figure
fig = figure('Position',[10 10 300 200])
% 3 days
subplot(2,1,1)
plot(timeshort,spikeshort, 'Color','black')
xlim([0,max(timeshort)])
xlabel('Time (days)')
ylabel('Spike rate')
title('Circadian cycles')
set(gca,'FontSize',fsize)
% 5 months
subplot(2,1,2)
plot(timemed,spikesmed,'Color','black')
xlim([0,max(timemed)])
xlabel('Time (days)')
ylabel('Spike rate')
title('Monthly+ cycles')
set(gca,'FontSize',fsize)

savename = 'H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\m5cyclesectiondemo.png';
exportgraphics(fig,savename,'Resolution',600)

% couple circadian days 
fig = figure('Position',[10,10,200,100]);
shadecolor = [0 0 0 0.2];
for jj = 1: length(starts)
    h1 = xline(timeshort(starts(jj)),'Color','white');
    hold on
    h2 = xline(timeshort(ends(jj)),'Color','white');
    hold on
    fill([timeshort(starts(jj)) timeshort(starts(jj)) timeshort(ends(jj)) timeshort(ends(jj))], [0 max(spikeshort) max(spikeshort) 0], 'blue','FaceAlpha',0.1,'EdgeAlpha',0.1);
end 
hold on
plot(timeshort,spikeshort, 'Color','black')
% plot(datesshort,spikeshort, 'Color','black')
xlim([0,max(timeshort)])
set(gca,'FontSize',fsize)
savename = 'H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\circadian.png';
exportgraphics(fig,savename,'Resolution',600)

% then just a long spike period figure too
fig = figure('Position',[10,10,1200,150]);
plot(spiketimedays(435:end) - spiketimedays(435),spikes(435:end),'Color','black')
xlim([0,max(spiketimedays)])
xlabel('Time (days)')
ylabel('Spikes per hour')
set(gca,'FontSize',fsize)

savename = 'H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\m5fullspikesdemo.png';
exportgraphics(fig,savename,'Resolution',600)

% example of a single spike
s = readmatrix('H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\spikeexamp.csv');
sig = s(2,2:end)./1000;
tvec = s(3,2:end);
fig = figure('Position',[10,10,150,100]);
% ylim([-500 800])
plot(tvec,sig,'Color','black')
set(gca,'FontSize',fsize)

savename = 'H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\examplespike.png';
exportgraphics(fig,savename,'Resolution',600)

%% Figure 2

% outputting sinusoids to illustrate what the assumed signal looks like

subject = 'simulation';
subsamplebasis = 'yes';
zeropadding = 'no';
pmin = 1;
pmax = 120;
maxsigamp = [15,5,10,10];
poi = [1, 7, 30, 70]; %  DEFINE PERIODS IN SIM SIGNALperiods (in days) of interest to include in simulated signal
poi = [3, 7, 30, 70]; %  DEFINE PERIODS IN SIM SIGNALperiods (in days) of interest to include in simulated signal

polycoef = [-.005 5*10^-6 0  0 0]; % for first, second, third, fourth, fifth, 
delta= 1;
dt = 1*60*60; %%% UDPATE THIS DT TO A COMPARABLE ONE FROM PATIENTS
noiseparams = 0.25:0.5:(5+0.25);

% DEFINE TIME ARRAY AND RANDOM SAMPLE SPACIN
years = 0.5'; % was 1.5
numsamp = 500;
% numsamp = 860;
timedurationseconds = years * 365 * 24* 60* 60;
fulltime = 1:timedurationseconds;
t = sort(datasample(fulltime, numsamp, 'Replace',false)); % randomly sample the full seconds timeseries
diffs = diff(t);
toremove = find(diffs < dt);% convert one hour to seconds
t(toremove) = [];
 % Adjust sample spacing (fill in zeros) based on temporal spacing
T_0 = min(t); % scale sample time array (based on distance in time from start-time)
t_hat = floor((t-T_0)/dt) + 1; % indices for when each time data were sampled
phi = t_hat;
t_full = min(t_hat): max(t_hat); % initialize full time array (temporal context, as if regularly sampled)

nn = 10;
noiseparam = noiseparams(nn);
    
% SIMULATE SIGNAL WITH NOISE
Fs = 1/dt; % sampling period and rate are inverse 
n = length(t_full);
N=n;
f = Fs*(0:n-1)/n;
f_adj = 1./(f*24*60*60); 
foi = [];
for i = 1: length(poi)
    [val, idx] = min(abs(f_adj-poi(i)));
    foi(i) = f(idx);
end 

tmod = linspace(min(t),max(t), length(t_full));
basesig = zeros(1,length(t_full));
components = [];
for i= 1:length(foi)
    component = maxsigamp(i).*sin(2*pi*foi(i)*tmod); 
    basesig = basesig + component;
    components(i,:) = component;
end
addednoise = noiseparam.*randn(1, length(t_full));
t4pol = [1:length(tmod)];
polynomialcomponent = polycoef(1)*t4pol + polycoef(2)*t4pol.^2 + polycoef(3)*t4pol.^3 +polycoef(4)*t4pol.^4 +polycoef(5)*t4pol.^5;
noisecomponent = noiseparam.*randn(1, length(t_full));
simsig = basesig + noisecomponent + polynomialcomponent; % adding noise
measureddata = simsig(phi)';

ind=4000;

fig=figure
subplot(4,1,1)
plot(basesig(1:ind), 'Color','black')
ylim([-40,40])

subplot(4,1,2)
plot(polynomialcomponent(1:ind), 'Color','black','LineWidth',1)

subplot(4,1,3)
plot(noisecomponent(1:ind), 'Color','black')
ylim([-40,40])

subplot(4,1,4)
plot(simsig(1:ind), 'Color','black')
ylim([-55,90])

% saveas(fig,'H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\ExampleSignalComponents.png');
print(gcf,'H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\ExampleSignalComponents.png','-dpng','-r600');   

F=fft(eye(n));

% basic example of single polynomial trend
fig = figure('Position',[10 10 100 100])
plot(polynomialcomponent(1:ind), 'Color','black','LineWidth',1)
ylim([-20, 65])
savename = 'H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\polytrend.tiff';
exportgraphics(fig,savename,'Resolution',600)

% example of sparsely sampled OG signal
fig=figure('Position',[10 10 500 80])
sig = simsig;
tveg = 1:length(sig);
tosamp = 1:length(sig);
tosamp(200:000) = 0;% add drops to tosamp
tosamp(500:550) = 0;
tosamp(850:1000) = 0;
tosamp(1500:1600) = 0;
tosamp(2600:2650) = 0;
tosamp(3000:3500) = 0;
tosamp(tosamp == 0) = [];
ind = datasample(tosamp,500);
plot(tvec(ind),sig(ind),'.', 'Color','black')
savename = 'H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\sparsesamplesim.tiff';
exportgraphics(fig,savename,'Resolution',600)

% series of oscillations (high, med, low freq)
ind = 1500;
fig=figure('Position',[10 10 80 150]);
subplot(3,1,1)
plot(components(1,1:ind),'Color','black')
subplot(3,1,2)
plot(components(2,1:ind),'Color','black')
ylim([-11 11])
subplot(3,1,3)
plot(components(3,1:ind),'Color','black')
ylim([-20 20])
savename = 'H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\coupleosc.tiff';
exportgraphics(fig,savename,'Resolution',600)

% example of just noise 
ind=500;
fig = figure('Position',[10 10 100 100])
plot(noisecomponent(1:ind), 'Color','black')
savename = 'H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\whitenoise.tiff';
exportgraphics(fig,savename,'Resolution',600)

% example of signal combined
fig = figure('Position',[10 10 600 80])
plot(simsig, 'Color','black')
ylim([-55,90])
savename = 'H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\combosignal.tiff';
exportgraphics(fig,savename,'Resolution',600)

% example of samples from signal combined 
sig = simsig;
tvec = 1:length(sig);
tosamp = 1:length(sig);
tosamp(200:000) = 0;% add drops to tosamp
tosamp(500:550) = 0;
tosamp(850:1000) = 0;
tosamp(1500:1600) = 0;
tosamp(2600:2650) = 0;
tosamp(3000:3500) = 0;
tosamp(tosamp == 0) = [];
ind = datasample(tosamp,100);
fig = figure('Position',[10 10 600 80])
plot(tvec, simsig, 'Color',[0 0 0 0.2])
hold on
plot(tvec(ind),simsig(ind),'.','Color','black') 
ylim([-55,90])
savename = 'H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\combosignal_withsamp.tiff';
exportgraphics(fig,savename,'Resolution',600)


subject ='M5';

% Load reference data from parameter sweep (N, dt, data, etc.) + sighits
destfolder = strcat('P:\Personal\Irena\MoodCyclesSpikeCycles_Paper\realdata\',subject,'\modeloutput\');
load(strcat(destfolder,'Outputfor_MSE_percent2drop10_N1500_traintest_forfig.mat'));

fsize = 6;
    
% basic example of BPDN spectrum
desiredperioddays = deltaselection.desiredperioddays;
x = deltaselection.fig.x;
fig = figure('Position',[10 10 50 50])
plot(desiredperioddays,x.^2,'Color','black','LineWidth',1)
xlim([-1 40])
savename = 'H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\minispectrum.tiff';
exportgraphics(fig,savename,'Resolution',600)

% example of DCT basis 
fig = figure('Position',[10 10 100 100])
simplebasis = dct(eye(100,100));
imagesc(simplebasis)
colormap('gray')
savename = 'H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\dctbasis.tiff';
exportgraphics(fig,savename,'Resolution',600)


%% Figure 3

% making a figure of outputs from performance in simulation - for paper 
% INCLUDED IN PLOT = LIKE REAL DATA ALL IN ONE FIGURE (minus seizure phase because NA)
% OG data
% spectrogram of OG data
% wavelet spec of OG data
% subsampling of og data samples
% BPDN-based reconstruction
% cwt over BPDN spectral comparison

% Jumbo - all in one per subject figure
% for publication
% has raw data, cwt spec, bpdn, recons, spectral overlay, slow seizure
% cycles, and phase plots 

subject = 'S1';
includepoly = 0;

fig=figure('Position',[10 10 1500 2000]);
fsize = 10;

% PART 1: BASIC RAW DATA AND CWT SPECTROGRAM 
duration = 12;
N = 3000;
numsamp = 1800;
destfolder = strcat('P:\Personal\Irena\BPDN_methods_paper\simulateddata\',subject,'\modeloutput\');
% filename = strcat(destfolder,'Outputfor',num2str(duration),'_N',num2str(N),'_numsamp',num2str(numsamp),'_forfig.mat');
%filename = strcat(destfolder,'Outputfor',num2str(duration),'_N',num2str(N),'_numsamp',num2str(numsamp),'_traintest_forfig.mat'); % original
filename = strcat(destfolder,'Outputfor',num2str(duration),'_N',num2str(N),'_numsamp',num2str(numsamp),'_traintest_forfig_1000x_test.mat'); % 1000x comparison

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
power = zscore(abs(deltaselection.fig.WT)); % MIGHT GO BACK AND RESCALE
wave_fam = 1./deltaselection.fig.fwt;
fs = 1 / deltaselection.dt;

% original''; signal 
subplot(4,6,[1,2,3,4])
plot(time_full_days,signal, 'Color','black','DisplayName','raw data')
ylabel('Spike rate')
%     title(strcat(subject,{' '},'raw spike rate'))
xlim([0 max(time_full_days)])
colorbar()
% ylim([0 max(signal) + 2000])
% legend()
ax=gca; ax.YAxis.Exponent = 3;
set(gca,'FontSize',fsize)

% spectrogram for cwt
subplot(4,6,[7,8,9,10])
contourf(time_full_days, log(perioddays_cwt), power,'edgecolor','none')
% xlabel('Time (days)')
ylabel('Period (days)')
%     title({'';'Spectrogram, CWT of raw spike rate'})
xlim([0 max(time_full_days)])
logged = log(perioddays_cwt);
inds = [49,98];
yticks([logged(inds)])
yticklabels({num2str(round(perioddays_cwt(inds)))})
ylim([logged(42) max(logged)])
colormap(parula);
colorbar()
set(gca,'FontSize',fsize)

% Delta selection figure 
subplot(4,6,[5,6,11,12])
colors = {'red','black','blue','cyan','green','magenta'};
for nn = 1: length(deltaselection.sampperdays)

    xx = deltaselection.deltas;
    dat = deltaselection.MSE{nn}'; % Number of ‘Experiments’ In Data Set
    yMean = mean(dat,1);  % Mean Of All Experiments At Each Value Of ‘x’
    ySEM = std(dat,0,1)/sqrt(size(dat,1));   % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
    CI95 = tinv([0.025 0.975], (size(dat,1)-1)); % Calculate 95% Probability Intervals Of t-Distribution
    yCI95 = bsxfun(@times, ySEM, CI95(:)); % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’

    %plot(x, yCI95+yMean, 'Color','blue') % Plot 95% Confidence Intervals Of All Experiments
    CIs = yCI95+yMean;
    curve1 = CIs(1,:);
    curve2 = CIs(2,:);
    xs = [xx, fliplr(xx)];
    inbetween = [curve1,fliplr(curve2)];
    fill(xs,inbetween, colors{nn},'FaceAlpha',0.5,'DisplayName','','EdgeColor','none','HandleVisibility','off')
    hold on
    labelname = strcat(num2str(deltaselection.sampperdays(nn)),{' '},'samples per day');
    plot(xx, yMean, 'Color',colors{nn}, 'LineWidth',2,'DisplayName',labelname{:}) % Plot Mean Of All Experiments
    hold on
    grid
end 
xlabel('Delta')
ylabel('MSE')
title(strcat(subject,{' '},'75/25 cross validation'))
legend
set(gca,'FontSize',fsize)

% PART 2: BPDN SAMPLES AND MODEL OUTPUT

% original''; signal with samples
subplot(4,6,[13,14,15,16])
% plot(time_full_days,signal, 'Color','black', 'DisplayName','raw data')
% hold on
plot(t_days, measureddata_wds,'.','Color','black', 'DisplayName','sparse samples','MarkerSize',2.5)
%     if s == length(subjects)
%         xlabel('Time (days)')
%     end
ylabel('Spike rate')
%     title({'';''})
xlim([0 max(time_full_days)])
ax=gca; ax.YAxis.Exponent = 3;
colorbar()
% ylim([0 max(signal) + 2000])
% legend()
set(gca,'FontSize',fsize)

% recon overlayed on original signal (or recon alone)
subplot(4,6,[19,20,21,22])
plot(time_full_days,signal, 'Color',[0 0 0 0.05],'DisplayName','raw data')
hold on
plot(t_recon_days, reconsig_wds,'Color',[0.8500 0.3250 0.0980],'DisplayName','BPDN recon')
% xlabel('Time (days)')
ylabel('Spike rate')
xlim([0 max(time_full_days)])
%     title({'';'Reconstructed signal, BPDN'})
ax=gca; ax.YAxis.Exponent = 3;
xlim([0 max(time_full_days)])
% ylim([0 max(reconsig_wds) + 1000])
% legend()
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
% xlabel('Period (days)')
ylabel('Power, BPDN')
set(gca,'FontSize',fsize)
% spectrum for cwt
yyaxis left
plot(perioddays_cwt, meanpowerCWTreg, 'Color',[0 0 0 0.3], 'LineWidth',1.5); % blue color with 0.6 alpha
% hold on
% plot(perioddays_cwt(locs),pks,'Color','cyan','o', 'MarkerSize',2)
xlabel('Period (days)')
ylabel({'';'';'Power,CWT'})
%     title('BPDN and CWT spectra')
ax=gca; ax.YAxis(1).Exponent = 5; ax.YAxis(1).Color = 'black';
xlim([-10 max(perioddays_cwt)])
set(gca,'FontSize',fsize)

savename = strcat('H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\',subject,'_allinone.tiff');
exportgraphics(fig,savename,'Resolution',600)

% plot the noise floor ------ Noise floor insert 
fig = figure('Position',[10,10,400, 200]);
for iii = 1:size(deltaselection.fig.xreshuffled,1)
    plot(desiredperioddays,deltaselection.fig.xreshuffled(iii,:).^2,'-','Color',[0.4940 0.1840 0.5560 0.2])
    hold on
end 
plot(desiredperioddays, x.^2,'-', 'Color',[0.8500 0.3250 0.0980], 'LineWidth',1);
hold on
plot(desiredperioddays(sighits),x(sighits).^2,'*','Color','black','MarkerSize',4)
set(gca,'YAxisLocation','right','ycolor',[0.8500 0.3250 0.0980])
xlim([0 80])
ylim([0 10*10^6])
savename = strcat('H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\',subject,'_allinone_noisefloorinsert.tiff');
exportgraphics(fig,savename,'Resolution',1200)

%% Figure 4

% for a real subject, plot...raw data, spectrogram, CWT, BPWP output, BPWP
% recon, slow filter on recon and seizure periodicities/polar plots

subject = 'M5';
includepoly = 0;

fig=figure('Position',[10 10 1500 2000]);
fsize = 10;

% PART 1: BASIC RAW DATA AND CWT SPECTROGRAM 
duration = 12;
N = 3000;
numsamp = 1800;
destfolder = strcat('P:\Personal\Irena\BPDN_methods_paper\realdata\',subject,'\modeloutput\');
% filename = strcat(destfolder,'Outputfor',num2str(duration),'_N',num2str(N),'_numsamp',num2str(numsamp),'_forfig.mat');
filename = strcat(destfolder,'Outputfor',num2str(duration),'_N',num2str(N),'_numsamp',num2str(numsamp),'_traintest_forfig.mat');
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
power = zscore(abs(deltaselection.fig.WT)); % MIGHT GO BACK AND RESCALE
wave_fam = 1./deltaselection.fig.fwt;
fs = 1 / deltaselection.dt;

% original''; signal 
subplot(6,6,[1,2,3,4])
plot(time_full_days,signal, 'Color','black','DisplayName','raw data')
ylabel('Spike rate')
%     title(strcat(subject,{' '},'raw spike rate'))
xlim([0 max(time_full_days)])
colorbar()
% ylim([0 max(signal) + 2000])
% legend()
ax=gca; ax.YAxis.Exponent = 3;
set(gca,'FontSize',fsize)

% spectrogram for cwt
subplot(6,6,[7,8,9,10])
contourf(time_full_days, log(perioddays_cwt), power,'edgecolor','none')
% xlabel('Time (days)')
ylabel('Period (days)')
%     title({'';'Spectrogram, CWT of raw spike rate'})
xlim([0 max(time_full_days)])
logged = log(perioddays_cwt);
inds = [49,98];
yticks([logged(inds)])
yticklabels({num2str(round(perioddays_cwt(inds)))})
ylim([logged(42) max(logged)])
colormap(parula);
colorbar()
set(gca,'FontSize',fsize)

% spectrum for cwt
subplot(6,6,[5,6,11,12])
plot(perioddays_cwt, meanpowerCWTreg, 'Color','black', 'LineWidth',1.5);
xlabel('Period (days)')
ylabel({'';'';''; 'Average power'})
%     title('Average CWT spectrum')
title({'';''})
xlim([-10 max(perioddays_cwt)])
%     ylim([0, 2.5*10^5])
ax=gca; ax.YAxis.Exponent = 5;
set(gca,'FontSize',fsize)

% PART 2: BPDN SAMPLES AND MODEL OUTPUT

% original''; signal with samples
subplot(6,6,[13,14,15,16])
% plot(time_full_days,signal, 'Color','black', 'DisplayName','raw data')
% hold on
plot(t_days, measureddata_wds,'.','Color','black', 'DisplayName','sparse samples','MarkerSize',2.5)
%     if s == length(subjects)
%         xlabel('Time (days)')
%     end
ylabel('Spike rate')
%     title({'';''})
xlim([0 max(time_full_days)])
ax=gca; ax.YAxis.Exponent = 3;
colorbar()
% ylim([0 max(signal) + 2000])
% legend()
set(gca,'FontSize',fsize)

% recon overlayed on original signal (or recon alone)
subplot(6,6,[19,20,21,22])
plot(time_full_days,signal, 'Color',[0 0 0 0.05],'DisplayName','raw data')
hold on
plot(t_recon_days, reconsig_wds,'Color',[0.8500 0.3250 0.0980],'DisplayName','BPDN recon')
% xlabel('Time (days)')
ylabel('Spike rate')
xlim([0 max(time_full_days)])
%     title({'';'Reconstructed signal, BPDN'})
ax=gca; ax.YAxis.Exponent = 3;
xlim([0 max(time_full_days)])
% ylim([0 max(reconsig_wds) + 1000])
% legend()
colorbar()
set(gca,'FontSize',fsize)
% 
% % CALCULATE COEFFICIENT OF DETERMINATION (ordinary, unadjusted)
% % (a measure of explained variance, based on linear correlation)
% % https://www.mathworks.com/help/stats/coefficient-of-determination-r-squared.html
% downsamplefactor = floor(length(deltaselection.signal) / deltaselection.N);
% downsampledsig = downsample(deltaselection.signal, downsamplefactor);
% % remove any excess samples randomly
% excess = length(downsampledsig) - length(reconsig_wds);
% r = randperm(deltaselection.N,excess);
% downsampledsig(r) = [];
% [rsquared] = BPDN_explainedVarianceMetrics(reconsig_wds, downsampledsig)


% overlay BPDN and CWT spectra
subplot(6,6,[17,18,23,24])
% spectrum for BPDN
yyaxis right
% no rescaling
plot(desiredperioddays, x.^2,'-', 'Color',[0.8500 0.3250 0.0980], 'LineWidth',1);
hold on
plot(desiredperioddays(sighits),x(sighits).^2,'*','Color','black','MarkerSize',4)

% xlabel('Period (days)')
ylabel('Power, BPDN')
set(gca,'FontSize',fsize)

% spectrum for cwt
yyaxis left
plot(perioddays_cwt, meanpowerCWTreg, 'Color',[0 0 0 0.3], 'LineWidth',1.5); % blue color with 0.6 alpha
% hold on
% plot(perioddays_cwt(locs),pks,'Color','cyan','o', 'MarkerSize',2)
xlabel('Period (days)')
ylabel({'';'';'Power,CWT'})
%     title('BPDN and CWT spectra')
ax=gca; ax.YAxis(1).Exponent = 5; ax.YAxis(1).Color = 'black';
xlim([-10 max(perioddays_cwt)])
set(gca,'FontSize',fsize)

% PART 3: Seizures, slow cycles, and phase plots 

 % Load seizures, convert to days - some seizure data processing
if strcmp(subject,'M4')
    szfile = strcat('P:\Personal\Irena\RCS\fromdashboard\',subject,'\periodicity\',subject,'_leftsz_ToryProbabilities.mat');
else
    szfile = strcat('P:\Personal\Irena\RCS\fromdashboard\',subject,'\periodicity\',subject,'_periodicities_sz.mat');
end 
load(szfile)
sz_days = (double(sz)./10^6 - deltaselection.time(1)) ./ (24*60*60);
% find indices in time of seizures
szind = [];
for i = 1: length(sz_days)
    [val,loc] = min(abs(t_recon_days - sz_days(i)));
    szind(i) = loc;
end 
% cut off ends 
tocut = find(szind >= length(t_recon_days));
szind(tocut) = [];
sz_days(tocut) = [];

% Recon just the sighits (zero out everything except the x of interest)
singlecyclerecons = [];
for i = 1:length(sighits)
    xmini = zeros(length(x),1);
    xmini(sighits(i)) = x(sighits(i));
    if includepoly == 0 % if not including polytrend
        ztemp = zeros(deltaselection.maxdegree + 1,1);
    elseif includepoly == 1 % if including poly trend
        ztemp = deltaselection.fig.z;
    end
    reconmini = BPDN_reconsig(deltaselection.f_k, xmini, deltaselection.scalefactor, deltaselection.maxdegree, ztemp, deltaselection.t_full);  
    singlecyclerecons(i,:) = reconmini;
end 

% add neighboring significant cycles (high frequency block, low frequency
% block)
periods = deltaselection.desiredperioddays(sighits);
lowfreqind = find(periods > 2);
toplot = sum(singlecyclerecons(lowfreqind,:),1);
% plot seizures over raw spike rates
% subplot(6,6,[25,26,27,28])
% plot(time_full_days,signal, 'Color','black')
% hold on
% plot(sz_days,toplot(szind),'*','MarkerSize',4,'DisplayName','Seizure','Color','green','LineWidth',1.5)
% ylabel('Spike rate')
% %     title(strcat(subject,{' '},'raw spike rate'))
% xlim([0 max(time_full_days)])
% colorbar()
% ax=gca; ax.YAxis.Exponent = 3;
% set(gca,'FontSize',fsize)

% now plot slow cycles and seizures 
% subplot(6,6,[31,32,33,34])
subplot(6,6,[25,26,27,28,31,32,33,34])
plot(t_recon_days,toplot,'DisplayName','cycles > 2days','LineWidth',1.5,'Color',[0.8500 0.3250 0.0980])
hold on
plot(sz_days,toplot(szind),'o','MarkerSize',3,'DisplayName','sz','Color',[0 0.4470 0.7410].*0.25,'LineWidth',1)
%     hold on
%     plot(sz_days,toplot(szind),'o','MarkerSize',4,'Color','black','LineWidth',1)
% title(strcat(subject,', Seizures with sum of "slow cycles" (>2 days)'))
% xlabel('Time (days)')
ylabel('Spike rate')
colorbar()
ax=gca; ax.YAxis.Exponent = 3;
set(gca,'FontSize',fsize)
xlim([0,365])
% ylim([-650 max(toplot(szind))+ 600])
% legend()
xlabel('Time (days)')

savename = strcat('H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\',subject,'_allinone.png');
exportgraphics(fig,savename,'Resolution',1200)

% ---------- plot noise floor separately as 2 little cutouts ------
% ---- one shows high freq end, other shows overall noisefloor

% plot the noise floor overall
fig = figure('Position',[10,10,400, 200]);
for iii = 1:size(deltaselection.fig.xreshuffled,1)
    plot(desiredperioddays,deltaselection.fig.xreshuffled(iii,:).^2,'-','Color',[0.4940 0.1840 0.5560 0.2])
    hold on
end 
plot(desiredperioddays, x.^2,'-', 'Color',[0.8500 0.3250 0.0980], 'LineWidth',1);
hold on
plot(desiredperioddays(sighits),x(sighits).^2,'*','Color','black','MarkerSize',4)
set(gca,'YAxisLocation','right','ycolor',[0.8500 0.3250 0.0980])
xlim([0 110])
ylim([0 max(x(sighits).^2)/3])
savename = strcat('H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\',subject,'_allinone_noisefloorinsert.png');
exportgraphics(fig,savename,'Resolution',1200)

% plot the noise floor left end (to show high amp spike isn't just some
% artifact of high freq end)
fig = figure('Position',[10,10,400, 200]);
for iii = 1:size(deltaselection.fig.xreshuffled,1)
    plot(desiredperioddays,deltaselection.fig.xreshuffled(iii,:).^2,'-','Color',[0.4940 0.1840 0.5560 0.2])
    hold on
end 
plot(desiredperioddays, x.^2,'-', 'Color',[0.8500 0.3250 0.0980], 'LineWidth',1);
hold on
plot(desiredperioddays(sighits),x(sighits).^2,'*','Color','black','MarkerSize',4)
set(gca,'YAxisLocation','right','ycolor',[0.8500 0.3250 0.0980])
xlim([0 5])
%ylim([0 5*10^8])
ylim([0 max(x(sighits).^2)/3])
savename = strcat('H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\',subject,'_allinone_noisefloorinsert_highend.png');
exportgraphics(fig,savename,'Resolution',1200)

% -------- Separate figure for seizure phase plots -------------
% individual seizure phase plots
polarcolor = 'black';%[0 0.4470 0.7410].*0.5;
fig = figure('Position',[10 10 500 500])
nbins = 24;
% intro processing
cycleranges = {[0.995,1.005];[5,10];[15,25];[25,35]}; % in days
% Use recon-sig without polynomial component
reconsig_justosc =  idct_custom(deltaselection.f_k, deltaselection.fig.x, deltaselection.scalefactor);
% Get seizure phase phase at the different ranges per Leguia 
t = t_recon_days .* (60*60*24);

[filtereds, phase] = firls_nick(reconsig_justosc, fs, cycleranges);
% Get phase per seizure
szphases = phase(:,szind);

% for the one day phase plot - actually just circadian hours
subplot(2,2,1)
% HOURS circadian version
% how to notate phase with time of day for circadian cycle?
% get seizure hours in subject local time
if subject == 'M1'
    tz = 'America/Denver';
elseif subject == 'M4'
    tz = 'America/New_York';
elseif subject == 'M5'
    tz = 'America/Chicago';
end 
dates = datetime(sz./10^6, 'ConvertFrom','posixtime','TimeZone','UTC');
dates.TimeZone = tz;
hours = hour(dates);
% convert hours to phase?
hoursdeg = hours.*15; % convert to degrees 
polarhistogram(deg2rad(hoursdeg), deg2rad(0:15:360),'FaceColor',polarcolor,'FaceAlpha',0.5);
pax = gca;
pax.ThetaAxisUnits = 'radians';
pax.ThetaDir = 'counterclockwise';
pax.ThetaZeroLocation = 'right'; % from playing around, very top is midnight, bottom is noon. *** THIS ORIENTATION IS ONLY TRUE FOR 1D
pax.ThetaTickLabel = [];
alpha = deg2rad(hoursdeg); % convert to radians
[pval, m] = circ_otest(alpha,[],ones(size(alpha))); % from berens circstat toolbox
%   Computes Omnibus or Hodges-Ajne test for non-uniformity of circular data.
if pval < 0.001
%             titlestring = [strcat(num2str(cycleranges{f}(2)),'-',num2str(cycleranges{f}(1)),{' '},'day cycle');strcat('Omnibus test, p = ***')];
    titlestring = strcat('1 day',{' '}, '***');
elseif pval < 0.01
    titlestring = [strcat('1 day cycle');strcat('Omnibus test, p = **')];
elseif pval < 0.05
    titlestring = [strcat('1 day cycle');strcat('Omnibus test, p = *')];
elseif pval >= 0.05
    titlestring = [strcat('1 day cycle');strcat('Omnibus test, p = ns')];
end 
title(titlestring)
% thetaticks(0:pi/2:2*pi)
set(gca,'FontSize',fsize-2)

% one week phase plot
subplot(2,2,2)
f=2;
polarhistogram(szphases(f,:),nbins,'FaceColor',polarcolor,'FaceAlpha',0.5)
pax = gca;
pax.ThetaAxisUnits = 'radians';
pax.ThetaZeroLocation = 'left';
%     pax.ThetaTick = [0 0.5*pi pi 1.5*pi];
pax.ThetaDir = 'clockwise';
alpha = (szphases(f,:).*pi) ./180; % convert to radians
[pval, m] = circ_otest(alpha,[],ones(size(alpha))); % from berens circstat toolbox
%   Computes Omnibus or Hodges-Ajne test for non-uniformity of circular data.
if pval < 0.001
%             titlestring = [strcat(num2str(cycleranges{f}(2)),'-',num2str(cycleranges{f}(1)),{' '},'day cycle');strcat('Omnibus test, p = ***')];
    titlestring = strcat(num2str(cycleranges{f}(1)),'-',num2str(cycleranges{f}(2)),{' '},'d',{' '}, '***');
elseif pval < 0.01
    titlestring = [strcat(num2str(cycleranges{f}(1)),'-',num2str(cycleranges{f}(2)),{' '},'day cycle');strcat('Omnibus test, p = **')];
elseif pval < 0.05
    titlestring = [strcat(num2str(cycleranges{f}(1)),'-',num2str(cycleranges{f}(2)),{' '},'day cycle');strcat('Omnibus test, p = *')];
elseif pval >= 0.05
    titlestring = [strcat(num2str(cycleranges{f}(1)),'-',num2str(cycleranges{f}(2)),{' '},'day cycle');strcat('Omnibus test, p = ns')];
end 
title(titlestring)
thetaticks(0:pi/2:2*pi)
set(gca,'FontSize',fsize-2)

% two-three week phase plot
subplot(2,2,3)
f=3;
polarhistogram(szphases(f,:),nbins,'FaceColor',polarcolor,'FaceAlpha',0.5)
pax = gca;
% pax.ThetaAxisUnits = 'radians';
% %         pax.ThetaDir = 'clockwise';
% pax.ThetaZeroLocation = 'left';
pax.ThetaAxisUnits = 'radians';
pax.ThetaZeroLocation = 'left';
% pax.ThetaTickLabel = [];
%     pax.ThetaTick = [0 0.5*pi pi 1.5*pi];
pax.ThetaDir = 'clockwise';
alpha = (szphases(f,:).*pi) ./180; % convert to radians
[pval, m] = circ_otest(alpha,[],ones(size(alpha))); % from berens circstat toolbox
%   Computes Omnibus or Hodges-Ajne test for non-uniformity of circular data.
if pval < 0.001
%             titlestring = [strcat(num2str(cycleranges{f}(2)),'-',num2str(cycleranges{f}(1)),{' '},'day cycle');strcat('Omnibus test, p = ***')];
    titlestring = strcat(num2str(cycleranges{f}(1)),'-',num2str(cycleranges{f}(2)),{' '},'d',{' '}, '***');
elseif pval < 0.01
    titlestring = [strcat(num2str(cycleranges{f}(1)),'-',num2str(cycleranges{f}(2)),{' '},'day cycle');strcat('Omnibus test, p = **')];
elseif pval < 0.05
    titlestring = [strcat(num2str(cycleranges{f}(1)),'-',num2str(cycleranges{f}(2)),{' '},'day cycle');strcat('Omnibus test, p = *')];
elseif pval >= 0.05
    titlestring = [strcat(num2str(cycleranges{f}(1)),'-',num2str(cycleranges{f}(2)),{' '},'day cycle');strcat('Omnibus test, p = ns')];
end 
title(titlestring)
thetaticks(0:pi/2:2*pi)
set(gca,'FontSize',fsize-2)

% four week phase plot
subplot(2,2,4)
f=4;
polarhistogram(szphases(f,:),nbins,'FaceColor',polarcolor,'FaceAlpha',0.5)
pax = gca;
pax.ThetaAxisUnits = 'radians';
pax.ThetaZeroLocation = 'left';
% pax.ThetaTickLabel = [];
%     pax.ThetaTick = [0 0.5*pi pi 1.5*pi];
pax.ThetaDir = 'clockwise';
alpha = (szphases(f,:).*pi) ./180; % convert to radians
[pval, m] = circ_otest(alpha,[],ones(size(alpha))); % from berens circstat toolbox
%   Computes Omnibus or Hodges-Ajne test for non-uniformity of circular data.
if pval < 0.001
%             titlestring = [strcat(num2str(cycleranges{f}(2)),'-',num2str(cycleranges{f}(1)),{' '},'day cycle');strcat('Omnibus test, p = ***')];
    titlestring = strcat(num2str(cycleranges{f}(1)),'-',num2str(cycleranges{f}(2)),{' '},'d',{' '}, '***');
elseif pval < 0.01
    titlestring = [strcat(num2str(cycleranges{f}(1)),'-',num2str(cycleranges{f}(2)),{' '},'d');strcat('Omnibus test, p = **')];
elseif pval < 0.05
    titlestring = [strcat(num2str(cycleranges{f}(1)),'-',num2str(cycleranges{f}(2)),{' '},'d');strcat('Omnibus test, p = *')];
elseif pval >= 0.05
    titlestring = [strcat(num2str(cycleranges{f}(1)),'-',num2str(cycleranges{f}(2)),{' '},'d');strcat('Omnibus test, p = ns')];
end 
title(titlestring)
thetaticks(0:pi/2:2*pi)
set(gca,'FontSize',fsize-2)

% FINALIZE AND SAVE
savename = strcat('H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\',subject,'_polarplots4allinone.png');
exportgraphics(fig,savename,'Resolution',1200)


%% Figure 5 - Varying sample density and impact on model performance

% initialize
subject = 'M5';
droplengthsdesired = [1,3,5];
droplengthinds = [2,4,6]; % DEFINE INDICES BASED ON VECTOR - double check that correct

% load data   
locfolder = strcat('P:\Personal\Irena\BPDN_methods_paper\realdata\',subject,'\percentdetectionsVsampledensity\');
%     filename = strcat(locfolder,'MSE_2D_DataDuration12_PercentDrop25_repeats10_N3000_numsamp180   360   720  1080  1440  1800_withpercenthits.mat');
filename = strcat(locfolder,'MSE_2D_DataDuration12_PercentDrop25_repeats10_N3000_numsamp360  1080  1800stabletraintestdata_withpercenthits_sampledensity_outputs.mat');
load(filename);

% colorspec for offset subplot
c1 = linspace(0,0.85,12); % dark blue
c2 = linspace(0,0.325,12);
c3 = linspace(1,0.098,12); % yellow

colors = horzcat(c1',c2',c3');
bpdnoutputcolor = [0.8500 0.3250 0.0980];

% generate plot for each drop length (loop)
fsize=8;
rows = 2*length(droplengthsdesired) +2;

fig=figure('Position',[10 10 1500 2000]);

for ii = 1: length(droplengthsdesired)
    
    i = droplengthinds(ii);
    droplength = droplengthsdesired(ii);
    suptitle(strcat(subject,{' '},'output at varying sample density'))

    % spectrum "truth" from continuous signal without drops
    truecwtspec = deltaselection.cwt.meanpowerreg;
    truecwtspecperiods =  deltaselection.cwt.period_days;
    truesig = deltaselection.signal;
    truesigtime = (deltaselection.time - deltaselection.time(1))./ (60*60*24); % converted to days
    
    cwtinterpspec = deltaselection.cwt.meanpowerreg;
    cwtinterpspectperiods = deltaselection.cwt.period_days;
    
    % BPDN measureddata, recon, and spectrum
    bpdnsamples = deltaselection.measureddatas_100block{1,i}(:,1);
    bpdnsamplestime = (deltaselection.ts_100block{1,i}(:,1) - deltaselection.time(1)) ./ (60*60*24); % converted to days
    % reconstruct signal
    [bpdnrecon] = BPDN_reconsig(deltaselection.f_k, deltaselection.percenthits{1,i}.x(:,1), deltaselection.scalefactor, deltaselection.maxdegree, deltaselection.percenthits{1,i}.z(:,1), deltaselection.t_full);
    bpdnrecontime = deltaselection.t_full *deltaselection.dt ./(60*60*24);
    bpdnspec = deltaselection.percenthits{1,i}.x(:,1).^2; % CONFIRM THIS IS RIGHT
    sighits = deltaselection.percenthits{1,i}.x_sighits(:,1);% CONFIRM THIS IS RIGHT 
    sighits_ind = find(sighits == 1);
    bpdnspecperiods = deltaselection.desiredperioddays;
    
    % plot original signal, unmarked
    subplot(rows,6,[1,2,3,4,7,8,9,10])
    plot(truesigtime, truesig, 'Color','black')
    set(gca,'FontSize',fsize)
    xlim([0, 365])
    
    % plot performance fig for offset/component period
    subplot(rows,6,[5,6,11,12])
    % DEFINE DATA - average offset of nearest sig peak vs sampling density
    samplesperdays = deltaselection.percenthits{1,1}.numsampsperdays4sighits; 
    component_periods = deltaselection.cwt.peakperiods;
    numreps = []; 
    toplot_perc = zeros(length(samplesperdays), length(component_periods));
    toplot_lb = zeros(length(samplesperdays), length(component_periods)); % lower and upper bounds for 95% CI
    toplot_ub = zeros(length(samplesperdays), length(component_periods));
    for i4 = 1: length(samplesperdays)
        samplesperday = samplesperdays(i4);
        for i5 = 1:length(component_periods)
            component_period = component_periods(i5);
            scaledmean = log((mean(deltaselection.percenthits{1,i4}.foundosc_offsets_reshufflings(i5,:))/component_period) / deltaselection.percenthits{1,i4}.numreps); % SCALED TP by
            stde = std(log(deltaselection.percenthits{1,i4}.foundosc_offsets_reshufflings(i5,:)./component_period));
            n = deltaselection.percenthits{1,i4}.numreps;
            toplot_lb(i4,i5) = scaledmean - 1.96*(stde / sqrt(n));
            toplot_perc(i4,i5) = scaledmean;
            toplot_ub(i4,i5) = scaledmean + 1.96*(stde / sqrt(n));
        end 
    end 
    for iii = 1:length(component_periods)
        xconf = [samplesperdays samplesperdays(end:-1:1)] ;         
        yconf = [toplot_ub(:,iii)' fliplr(toplot_lb(:,iii)')];
        p = fill(xconf,yconf,colors(iii,:),'HandleVisibility','off');
        p.FaceColor = colors(iii,:);  
        p.FaceAlpha = 0.2;
        p.EdgeColor = 'none';   
        hold on
        plot(samplesperdays, toplot_perc(:,iii),'o-','DisplayName',strcat(num2str(round(component_periods(iii),2)),'d cycle CWT'),'LineWidth',2,'Color', colors(iii,:));
%         xlabel('Random samples per day')
%         ylabel('Average offset (days) / period')
        set(gca,'FontSize',fsize)
        hold on
    end 
    legend('Location','East')
     set(gca, 'YScale', 'log')
    ylim([-10, 0])
    xlim([0,9])
    
    % plot bpdn source signal
    subplot(rows,6,[13,14,15,16] + (ii-1)*12)
    plot(truesigtime, truesig, 'Color',[0 0 0 0.05])
    hold on
    plot(bpdnsamplestime, bpdnsamples,'*','Color',bpdnoutputcolor,'MarkerSize',1.5)
%     xlabel('Time (days)')
%     ylabel('Spikes per hour')
%     title({'';'Data input to BPDN'})
    set(gca,'FontSize',fsize)
    xlim([0, 365])
%     ylim([0,2000])
    
    % plot bpdn recon
    subplot(rows,6,[19,20,21,22]+ (ii-1)*12)
    plot(truesigtime, truesig, 'Color',[0 0 0 0.05])
    hold on
    plot(bpdnrecontime, bpdnrecon,'Color',bpdnoutputcolor)
%     xlabel('Time (days)')
%     ylabel('Spikes per hour')
%     title({'';'BPDN reconstructed signal'})
    set(gca,'FontSize',fsize)
    xlim([0, 365])
%     ylim([0,700])

    % plot bpdn spectrum (with gaps) and cwt overlay
    subplot(rows,6,[17,18,23,24] + (ii-1)*12)
    % spectrum for BPDN
    yyaxis right
    % no rescaling
    plot(bpdnspecperiods, bpdnspec, 'Color',bpdnoutputcolor,'DisplayName','BPDN', 'LineWidth',1)
    hold on
    plot(bpdnspecperiods(sighits_ind),bpdnspec(sighits_ind),'*','Color','black','MarkerSize',4)
    % xlabel('Period (days)')
%     ylabel('Power, BPDN')
%     ylim([-500 max(bpdnspec)+50])
    set(gca,'FontSize',fsize)
    % spectrum for cwt
    yyaxis left
    plot(truecwtspecperiods, truecwtspec, 'Color',[0 0 0 0.3],'DisplayName','CWT original signal','LineWidth',1)
    ax=gca;ax.YAxis(1).Exponent = 3;ax.YAxis(1).Color = 'black'; % ax.YAxis(1).Exponent = 5; 
%     xlabel('Period (days)')
%     ylabel({'';'';'Power,CWT'})
%     title('BPDN spectrum')
%     ylim([-500 max(truecwtspec)+50])
    xlim([-10 max(truecwtspecperiods)])
    set(gca,'FontSize',fsize)    
    
end 
savename = strcat('H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\',subject,'_varyingSamplesAllinone.tiff');
exportgraphics(fig,savename,'Resolution',600)


% Plot performance figure out separately to break the y axis

% the main block of lines
fig = figure('Position',[10, 10, 600,400]);
% DEFINE DATA - average offset of nearest sig peak vs sampling density
samplesperdays = deltaselection.percenthits{1,1}.numsampsperdays4sighits; 
component_periods = deltaselection.cwt.peakperiods;
numreps = []; 
toplot_perc = zeros(length(samplesperdays), length(component_periods));
toplot_lb = zeros(length(samplesperdays), length(component_periods)); % lower and upper bounds for 95% CI
toplot_ub = zeros(length(samplesperdays), length(component_periods));
for i4 = 1: length(samplesperdays)
    samplesperday = samplesperdays(i4);
    for i5 = 1:length(component_periods)
        component_period = component_periods(i5);
        scaledmean = log((mean(deltaselection.percenthits{1,i4}.foundosc_offsets_reshufflings(i5,:))/component_period) / deltaselection.percenthits{1,i4}.numreps); % SCALED TP by
        stde = std(log(deltaselection.percenthits{1,i4}.foundosc_offsets_reshufflings(i5,:)./component_period));
        n = deltaselection.percenthits{1,i4}.numreps;
        toplot_lb(i4,i5) = scaledmean - 1.96*(stde / sqrt(n));
        toplot_perc(i4,i5) = scaledmean;
        toplot_ub(i4,i5) = scaledmean + 1.96*(stde / sqrt(n));
    end 
end 
for iii = 2:length(component_periods)
    xconf = [samplesperdays samplesperdays(end:-1:1)] ;         
    yconf = [toplot_ub(:,iii)' fliplr(toplot_lb(:,iii)')];
    p = fill(xconf,yconf,colors(iii,:),'HandleVisibility','off');
    p.FaceColor = colors(iii,:);  
    p.FaceAlpha = 0.2;
    p.EdgeColor = 'none';   
    hold on
    plot(samplesperdays, toplot_perc(:,iii),'o-','DisplayName',strcat(num2str(round(component_periods(iii),2)),'d cycle CWT'),'LineWidth',2,'Color', colors(iii,:));
%         xlabel('Random samples per day')
%         ylabel('Average offset (days) / period')
    set(gca,'FontSize',fsize)
    hold on
end 
legend('Location','East')
 set(gca, 'YScale', 'log')
ylim([-9, -2.5])
xlim([0,8])
savename = strcat('H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\',subject,'_offsetinsert_main.tiff');
exportgraphics(fig,savename,'Resolution',600)

% the separate floater line 
fig = figure('Position',[10, 10, 600,400]);

% DEFINE DATA - average offset of nearest sig peak vs sampling density
samplesperdays = deltaselection.percenthits{1,1}.numsampsperdays4sighits; 
component_periods = deltaselection.cwt.peakperiods;
numreps = []; 
toplot_perc = zeros(length(samplesperdays), length(component_periods));
toplot_lb = zeros(length(samplesperdays), length(component_periods)); % lower and upper bounds for 95% CI
toplot_ub = zeros(length(samplesperdays), length(component_periods));
for i4 = 1: length(samplesperdays)
    samplesperday = samplesperdays(i4);
    for i5 = 1:length(component_periods)
        component_period = component_periods(i5);
        scaledmean = log((mean(deltaselection.percenthits{1,i4}.foundosc_offsets_reshufflings(i5,:))/component_period) / deltaselection.percenthits{1,i4}.numreps); % SCALED TP by
        stde = std(log(deltaselection.percenthits{1,i4}.foundosc_offsets_reshufflings(i5,:)./component_period));
        n = deltaselection.percenthits{1,i4}.numreps;
        toplot_lb(i4,i5) = scaledmean - 1.96*(stde / sqrt(n));
        toplot_perc(i4,i5) = scaledmean;
        toplot_ub(i4,i5) = scaledmean + 1.96*(stde / sqrt(n));
    end 
end 
for iii = 1
    xconf = [samplesperdays samplesperdays(end:-1:1)] ;         
    yconf = [toplot_ub(:,iii)' fliplr(toplot_lb(:,iii)')];
    p = fill(xconf,yconf,colors(iii,:),'HandleVisibility','off');
    p.FaceColor = colors(iii,:);  
    p.FaceAlpha = 0.2;
    p.EdgeColor = 'none';   
    hold on
    plot(samplesperdays, toplot_perc(:,iii),'o-','DisplayName',strcat(num2str(round(component_periods(iii),2)),'d cycle CWT'),'LineWidth',2,'Color', colors(iii,:));
%         xlabel('Random samples per day')
%         ylabel('Average offset (days) / period')
    set(gca,'FontSize',fsize)
    hold on
end 
legend('Location','NorthEast')
 set(gca, 'YScale', 'log')
ylim([-10, 0])
xlim([0,8])
savename = strcat('H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\',subject,'_offsetinsert_singleline.tiff');
exportgraphics(fig,savename,'Resolution',600)

%% Figure 6 - Varying gaps/drops in data and associated model performance

% initialize
subject = 'M5';
droplengthsdesired = [12,30,60];
droplengthinds = [2,4,5]; % DEFINE INDICES BASED ON VECTOR - double check that correctand matches the droplengths desired.

% load data 
locfolder = strcat('P:\Personal\Irena\BPDN_methods_paper\realdata\',subject,'\percentdetectionsVdatadrops\');
filename = strcat(locfolder,'MSE_2D_DataDuration12_PercentDrop25_repeats10_N3000_numsamp360  1080  1800stabletraintestdata_withpercenthits_datadrops_outputs_finalfix.mat');
load(filename);

% color spec for performance figure
c1 = linspace(0,0.85,12); % dark blue
c2 = linspace(0,0.325,12);
c3 = linspace(1,0.098,12); % yellow
 
colors = horzcat(c1',c2',c3');
bpdnoutputcolor = [0.8500 0.3250 0.0980];

% generate plot for each drop length (loop)
fsize=8;
rows = 2*length(droplengthsdesired) +2;

fig=figure('Position',[10 10 1500 2000]);
for ii = 1: length(droplengthsdesired)
    
    i = droplengthinds(ii);
    droplength = deltaselection.percenthits{1,1}.gaplengths(i);
    suptitle(strcat(subject,{' '},'output at varying drop length'))

    % spectrum "truth" from continuous signal without drops
    truecwtspec = deltaselection.cwt.meanpowerreg;
    truecwtspecperiods =  deltaselection.cwt.period_days;
    truesig = deltaselection.signal;
    truesigtime = (deltaselection.time - deltaselection.time(1))./ (60*60*24); % converted to days
    
    % continuous (interpolated) signal used for cwt, recon, and spectrum (picking just a single example from the 10 total)
    cwtinterpsig = deltaselection.cwtinterp.measureddatas_100block_cwt{1,i}{1,1};
    cwtinterpsigtime = (deltaselection.cwtinterp.ts_100block_cwt{1,i}{1,1} -deltaselection.time(1))./ (60*60*24); % converted to days
    cwtinterprecon = [];
    cwtinterprecontime = [];
    cwtinterpspec = deltaselection.cwtinterp.meanpowerreg{1,i}{1,1};
    cwtinterpspectperiods = deltaselection.cwtinterp.period_days{1,i}{1,1};
    
    % BPDN measureddata, recon, and spectrum
    bpdnsamples = deltaselection.measureddatas_100block{1,i}{1,1};
    bpdnsamplestime = (deltaselection.ts_100block{1,i}{1,1} - deltaselection.time(1)) ./ (60*60*24); % converted to days
    % reconstruct signal
    [bpdnrecon] = BPDN_reconsig(deltaselection.f_k, deltaselection.percenthits{1,i}.x(:,2), deltaselection.scalefactor, deltaselection.maxdegree, deltaselection.percenthits{1,i}.z(:,1), deltaselection.t_full);
    bpdnrecontime = deltaselection.t_full *deltaselection.dt ./(60*60*24);
    bpdnspec = deltaselection.percenthits{1,i}.x(:,2).^2;
    sighits = deltaselection.percenthits{1,i}.x_sighits(:,2); % CONFIRM ACCURACRY
    sighits_ind = find(sighits == 1);
    bpdnspecperiods = deltaselection.desiredperioddays;
    
    % plot original signal, unmarked
    subplot(rows,6,[1,2,3,4,7,8,9,10])
    plot(truesigtime, truesig, 'Color','black')
%     xlabel('Time (days)')
%     ylabel('Spikes per hour')
%     title({'';'Raw data'})
    set(gca,'FontSize',fsize)
    xlim([0, 365])

    % plot bpdn source signal
    subplot(rows,6,[13,14,15,16] + (ii-1)*12)
    plot(truesigtime, truesig, 'Color',[0 0 0 0.05])
    hold on
    plot(bpdnsamplestime, bpdnsamples,'*','Color',bpdnoutputcolor, 'MarkerSize',1.5)
%     xlabel('Time (days)')
%     ylabel('Spikes per hour')
%     title({'';'Data input to BPDN'})
    set(gca,'FontSize',fsize)
    xlim([0, 365])
    
    % plot bpdn recon
    subplot(rows,6,[19,20,21,22]+ (ii-1)*12)
    plot(truesigtime, truesig, 'Color',[0 0 0 0.05])
    hold on
    plot(bpdnrecontime, bpdnrecon,'Color',bpdnoutputcolor)
%     xlabel('Time (days)')
%     ylabel('Spikes per hour')
%     title({'';'BPDN reconstructed signal'})
    set(gca,'FontSize',fsize)
    xlim([0, 365])
    
    % plot bpdn spectrum (with gaps) and cwt overlay
    subplot(rows,6,[17,18,23,24] + (ii-1)*12)
    % spectrum for BPDN
    yyaxis right
    % no rescaling
    plot(bpdnspecperiods, bpdnspec, 'Color',bpdnoutputcolor,'DisplayName','BPDN signal with data drops','LineWidth',1)
    hold on
    plot(bpdnspecperiods(sighits_ind),bpdnspec(sighits_ind),'*','Color','black','MarkerSize',4)
    % xlabel('Period (days)')
%     ylabel('Power, BPDN')
    set(gca,'FontSize',fsize)
    % spectrum for cwt
    yyaxis left
    plot(truecwtspecperiods, truecwtspec, 'Color','black','DisplayName','CWT original signal')
    ax=gca;ax.YAxis(1).Exponent = 3;ax.YAxis(1).Color = 'black'; 
%     xlabel('Period (days)')
%     ylabel({'';'';'Power,CWT'})
%     title('BPDN spectrum (with drops)')
    xlim([-10 max(truecwtspecperiods)])
    set(gca,'FontSize',fsize)    
    
end 
savename = strcat('H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\',subject,'_varyingDropsAllinone.tiff');
exportgraphics(fig,savename,'Resolution',600)

%% Figure 7 (+ supplemental 10,11) - Method performance (frequency-specific correlation estimate vs original)

%Error/correlation between recon and original signal (filtered by band)

% Load data
subject = 'M4';
duration = 12;
N = 3000;
numsamp = 1800;
destfolder = strcat('P:\Personal\Irena\BPDN_methods_paper\realdata\',subject,'\modeloutput\');
filename = strcat(destfolder,'Outputfor',num2str(duration),'_N',num2str(N),'_numsamp',num2str(numsamp),'_methodperformance.mat');
load(filename);

fsize = 10;

% TO BUILD SCHEMATIC ILLUSTRATION -----------------

bandindex = 10; % band indices roughly correlated with central frequency
fig = figure('Position',[10,10,600, 200]);

% Plot of raw signal
subplot(2,2,1)
plot(performance.signal_time, performance.signal,'Color','black')
xlim([0 360])
ylim([0 6000])
set(gca,'FontSize',fsize)
ylabel('IED rate - original signal')

% Plot of reconstructed signal
subplot(2,2,3)
plot(performance.recon_time, performance.recon,'Color',[0.8500 0.3250 0.0980]) % this is the orange color used throughout
xlim([0 360])
ylim([0 6000])
set(gca,'FontSize',fsize)
ylabel('IED rate - reconstruction')

% Plot of filtered raw signal
subplot(2,2,2)
plot(performance.recon_time, performance.signalfiltered(bandindex,:),'-','Color','black','LineWidth',1);
xlim([0 360])
ylim([-1500 2500])
set(gca,'FontSize',fsize)
ylabel(strcat('Filt',{' '},num2str(performance.bands{bandindex})))

% Plot of filtered reconstructed signal
subplot(2,2,4)
plot(performance.recon_time,performance.reconfiltered(bandindex,:),'-','Color',[0.8500 0.3250 0.0980],'LineWidth',1);
xlim([0 360])
ylim([-1500 2500])
set(gca,'FontSize',fsize)
ylabel(strcat('Filt',{' '},num2str(performance.bands{bandindex})))

savename = strcat('C:\Users\m164085\OneDrive - Mayo Clinic\H_Mig\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\','method_performance_methodschematic_',subject,'.tiff');
exportgraphics(fig,savename,'Resolution',1200)

% TO PLOT OVERALL PERFORMANCE -------------------------------

fig = figure('Position',[10,10,300, 200]);
plot(performance.centralfreqs,performance.corrcoefs,'-','Color','black')
hold on
plot(performance.centralfreqs,performance.corrcoefs,'.','MarkerSize',10,'Color','black')
ylabel('Correlation coeff')
ylim([0 1.05])
xlabel('Central Period (days)')
set(gca,'FontSize',fsize)
savename = strcat('H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\method_performance_corrcoefvscentralfreq_',subject,'.tiff');
exportgraphics(fig,savename,'Resolution',1200)

% TO PLOT EXAMPLES -------------------------------------------

fig = figure('Position',[10,10,500, 300]);

% high freq
bi = 1;
subplot(3,1,1)
plot(performance.recon_time, performance.signalfiltered(bi,:),'Color','black','LineWidth',1);
hold on
plot(performance.recon_time,performance.reconfiltered(bi,:),'Color',[0.8500 0.3250 0.0980],'LineWidth',1);
ylabel(strcat('bp',{' '},num2str(performance.bands{bi})))
xlim([0 360])
ylim([-1200 2000])
set(gca,'FontSize',fsize)

% med freq
bi = 10;
subplot(3,1,2)
plot(performance.recon_time, performance.signalfiltered(bi,:),'Color','black','LineWidth',1);
hold on
plot(performance.recon_time,performance.reconfiltered(bi,:),'Color',[0.8500 0.3250 0.0980],'LineWidth',1);
ylabel(strcat('bp',{' '},num2str(performance.bands{bi})))
xlim([0 360])
set(gca,'FontSize',fsize)

% low freq
bi = 50;
subplot(3,1,3)
plot(performance.recon_time, performance.signalfiltered(bi,:),'Color','black','LineWidth',1);
hold on
plot(performance.recon_time,performance.reconfiltered(bi,:),'Color',[0.8500 0.3250 0.0980],'LineWidth',1);
ylabel(strcat('bp',{' '},num2str(performance.bands{bi})))
xlim([0 360])
set(gca,'FontSize',fsize)

savename = strcat('H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\method_performance_examplesperband_',subject,'.tiff');
exportgraphics(fig,savename,'Resolution',1200)

% Zoomed-in version of 1 day cycle (fewer days represented)
fig = figure('Position',[10,10,400, 100]);
bi = 1;
plot(performance.recon_time, performance.signalfiltered(bi,:),'Color','black','LineWidth',1);
hold on
plot(performance.recon_time,performance.reconfiltered(bi,:),'Color',[0.8500 0.3250 0.0980],'LineWidth',1);
ylabel(strcat('bp',{' '},num2str(performance.bands{bi})))
xlim([0 25])
set(gca,'FontSize',fsize)
savename = strcat('H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\method_performance_examplesperband_zoomedin1day_',subject,'.tiff');
exportgraphics(fig,savename,'Resolution',1200)

%% Supplemental figure 1 - 75/25 cross validation output

subjects = {'S1'}; % etc

colors = {'red','black','blue','cyan','green','magenta'};
fsize = 15;
fig = figure('Position',[10,10,700,1000]);

for i = 1:length(subjects)
    
    subplot(1,length(subjects),i)
    subject = subjects{i};
   
    filename =  strcat('P:\Personal\Irena\BPDN_methods_paper\simulateddata\',subject,'\deltaselection\MSE_2D_DataDuration12_PercentDrop25_repeats10_N3000_numsamp360  1080  1800stabletraintestdata.mat');
    load(filename);
    
    for nn = 1: length(deltaselection.sampperdays)
    
        x = deltaselection.deltas;
        dat = deltaselection.MSE{nn}'; % Number of ‘Experiments’ In Data Set
        yMean = mean(dat,1);  % Mean Of All Experiments At Each Value Of ‘x’
        ySEM = std(dat,0,1)/sqrt(size(dat,1));   % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
        CI95 = tinv([0.025 0.975], (size(dat,1)-1)); % Calculate 95% Probability Intervals Of t-Distribution
        yCI95 = bsxfun(@times, ySEM, CI95(:)); % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’

        %plot(x, yCI95+yMean, 'Color','blue') % Plot 95% Confidence Intervals Of All Experiments
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
end 

savename = 'H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\deltaparametersweeps_simdata.tiff';
exportgraphics(fig,savename,'Resolution',600)

%% Supplemental figure 2

% CREATE SIMULATED SIGNALS ---------------------------------------------
% Note: can customize looped signal params to have variance, SNR, etc. 

% Fixed signal params
polycoef = [-0.0025 -10^-7 0 0 0]; 
samplewindowduration = 16;% hours (window during the day when samples can be collected)

% Destination folder
destfolder = 'P:\Personal\Irena\BPDN_methods_paper\simulated_signals_focusOnSparseness\';

% Looped signal params
component_periods = [5,50,100]; 
numresamp = 100; % number of y and t_hats to generate
durations = [12];%[3, 12, 36]; % months
phaseshifts = [0];%[0,1]; % yes/no, is there a phase shift in the middle of the timeseries?
samplingtypes = {'irregular'}; %{'irregular','regular'}; % Irregular is still just sampled during the day
samplesperdays = [1/7, 4/7];%[2, 1, 3/7, 1/7]; % = sparsity essentially 
SNR_models = [1, 5]; % based on amplitude of signal vs amplitude of added noise 
variances = [1, 10]; % this technically signal amplitude - ties directly to variance 
dtinput = []; % this will default signal building to initial 15 min dt
numsamp = []; % this will default to samples per day as definition
signaltype = 'oscillation';
percentcontaining = 100;

% loop through loops
for i1 = 1:length(durations)
    duration = durations(i1);
    for i2 = 1:length(phaseshifts)
        phaseshift = phaseshifts(i2);
        for i3 = 1:length(samplingtypes)
            samplingtype = samplingtypes{i3};
            for i4 = 1:length(samplesperdays)
                samplesperday = samplesperdays(i4);
                for i5 = 1:length(SNR_models)
                    SNR_model = SNR_models(i5);
                    for i6 = 1:length(variances)
                        variance = variances(i6);
                        for i7 = 1: length(component_periods)
                            poi = component_periods(i7);

                            % build signal function 
                            [signal, t_full] = build_simsig_BPDN(signaltype, poi, polycoef, duration, phaseshift, SNR_model, variance, dtinput, percentcontaining);

                            % sample signal function
                            y = [];
                            t_hat = [];
                            for rr = 1: numresamp
                                [y_new, t_hat_new] = sample_simsig_BPDN(samplingtype, samplewindowduration,samplesperday, numsamp, dtinput,  t_full, signal);
                                y(:,rr) = y_new;
                                t_hat(:,rr) = t_hat_new;
                            end 
                            simsig.component_osc_periods_days = poi;
                            simsig.polynomial_coefficients = polycoef;
                            simsig.months = duration;
                            simsig.samplewindowduration = samplewindowduration;
                            simsig.phaseshift = phaseshift; 
                            simsig.samplingtype = samplingtype;
                            simsig.samplesperday = samplesperday; 
                            simsig.SNR_model = SNR_model;
                            simsig.variance = variance; 
                            simsig.signal = signal;
                            simsig.t_full = t_full; % array in seconds
                            simsig.t_hat = t_hat; % indices of samples, array
                            simsig.y = y; % samples associated with t_hat, array

                            filename = strcat(destfolder,signaltype,'Duration',num2str(duration),'months_phaseshift',num2str(phaseshift),'_samplingtype',samplingtype,'_sampperday',num2str(samplesperday),'_SNRmodel',num2str(SNR_model), '_variance',num2str(variance),'_osc',num2str(poi),'days.mat'); % FINISH FILE NAME 
                            save(filename, 'simsig'); 
                    
                        end
                    end 
                end 
            end
        end
    end
end
% 

% ANALYZE SIMULATED SIGNALS ---------------------------------------------
% Simulation analysis script for percent hits vs frequency: APPLIED TO EVERY SIMULATED SIGNAL

% MADE SOME CHANGES TO GET IT TO RUN FASTER
% changed bpdn reshuffling temporalily to only do 100 repeats instead of
% 1000
% CHANGE BACK WHEN DONE.... 
% chose to fix dt, commented out sparsity, but don't have to. 

locfolder = 'P:\Personal\Irena\BPDN_methods_paper\simulated_signals_focusOnSparseness\'; % CHANGE BACK
percentilethreshold = 99; 
amplitudepercentilethreshold = []; 
sigpeakmax = 10;
delta = 7500; % selected based on parameter sweeps
dt = 5*60*60; % fixed at one hour. check once and see how long that makes things *******************
maxperiod = 120; 
minperiod = 10; 
maxdegree = 2;
sampletype = 'densesubrange';

% JUST RUN THROUGH ALL THE SIMULATED SIGNALS
files = dir(locfolder);
files(1:2) = []; 

for i = 1: length(files) 
    
    tic
    loopfilename = strcat(files(i).folder,'\',files(i).name); 
    load(loopfilename); 

    % LOOP through the number of y and t_hat resamplings in simsig
    foundosc_offsets_repeats = [];
    foundosc_offsets_repeats_amp = [];
    x_sighits_mat = [];
    x_reshufflings = [];
    x_percentiles = [];
    z_sighits_mat = [];
    z_reshufflings = [];
    z_percentiles = [];
    TPs = [];
    FPs = [];
    TNs = [];
    FNs = [];
    TPRs = [];
    FPRs = [];
    precisions = [];
    recalls = [];
    
    TPs_a = [];
    FPs_a = [];
    TNs_a = [];
    FNs_a = [];
    TPRs_a = [];
    FPRs_a = [];
    precisions_a = [];
    recalls_a = [];
     
    for rr = 1: size(simsig.y,2) 

        % DEFINE BASES, PHI, ETC.
        %dt = ((simsig.t_hat(end,rr) - simsig.t_hat(1,rr)) * sparsity)/ length(simsig.y(:,rr));
        [phi, N, t_full, f_k, desiredperioddays, dct_basis, scalefactor, poly_basis] = BPDN_fullSetup(simsig.t_hat(:,rr),dt,maxperiod,minperiod,sampletype,maxdegree);
        [x, xreshuffled, xpercents, z, zreshuffled, zpercents,reshuff_optval] = BPDN_wReshuffling_short(delta, simsig.y(:,rr), phi, dct_basis, poly_basis);
        x_reshufflings(:,:,rr) = xreshuffled; % will be 3d matrix
        z_reshufflings(:,:,rr) = zreshuffled;
        x_percentiles(:,rr) = xpercents;
        z_percentiles(:,rr) = zpercents;
        
        % DETERMINE PERCENT SIGNIFICANCE BASED ON PERCENTILE
        % THRESHOLD-------------------------------------------------------
        [x_sighits] = percentilehits(xpercents,percentilethreshold);
        x_sighits_mat(:,rr) = x_sighits; % TODO CONFIRM DIMENSIONS
        [z_sighits] = percentilehits(zpercents,percentilethreshold);
        z_sighits_mat(:,rr) = z_sighits; % TODO CONFIRM DIMENSIONS

        % HOW FAR WERE SIG HITS (PEAKS) FROM KNOWN PERIODs (closest found peak)
        realosc = simsig.component_osc_periods_days;
        [foundosc] = find(x_sighits ==1);
        if realosc < 1 
            foundosc = round(desiredperioddays(foundosc),1); % rounding here is OK for whole freq, not 0.5 periods.CONFIRM rounding is OK******
        elseif realosc >=1
            foundosc = round(desiredperioddays(foundosc)); % rounding here is OK for whole freq, not 0.5 periods.CONFIRM rounding is OK******
        end
        foundosc = unique(foundosc);
        oscdifferences = []; % find distance of real osc from identified osc per real osc
        for p = 1: length(realosc)
            try
                oscdifferences(p) = min(abs(foundosc - realosc(p)));
            catch 
                oscdifferences(p) = NaN;
            end
        end
        foundosc_offsets_repeats(:,rr) = oscdifferences; 

        % GET TRUE POS + FALSE POS etc.
        [TP,FN,TN,FP, TPR, FPR, precision, recall] = get_OnOffTargetPeaks(x_sighits, desiredperioddays, simsig.component_osc_periods_days);
        
        TPs(rr) = TP; 
        FPs(rr) = FP;
        FNs(rr) = FN;
        TNs(rr) = TN;
        TPRs(rr) = TPR;
        FPRs(rr) = FPR;
        precisions(rr) = precision;
        recalls(rr) = recall;
        
        % DETERMINE PERCENT SIGNIFICANCE BASED ON PERCENTILE THRESHOLD AND
        % ADDED AMPLITUDE THRESHOLD -------------------------------------
        [x_sighits_amp] = percentilehits_withamplitudethreshold(xpercents, percentilethreshold,x, amplitudepercentilethreshold, sigpeakmax);
        x_sighits_amp_mat(:,rr) = x_sighits_amp; % TODO CONFIRM DIMENSIONS
        
        % HOW FAR WERE SIG HITS (PEAKS) FROM KNOWN PERIODs (closest found peak)
        realosc = simsig.component_osc_periods_days;
        [foundosc] = find(x_sighits_amp ==1);
        if realosc < 1 
            foundosc = round(desiredperioddays(foundosc),1); % rounding here is OK for whole freq, not 0.5 periods.CONFIRM rounding is OK******
        elseif realosc >=1
            foundosc = round(desiredperioddays(foundosc)); % rounding here is OK for whole freq, not 0.5 periods.CONFIRM rounding is OK******
        end
        foundosc = unique(foundosc);
        oscdifferences = []; % find distance of real osc from identified osc per real osc
        for p = 1: length(realosc)
            try
                oscdifferences(p) = min(abs(foundosc - realosc(p)));
            catch 
                oscdifferences(p) = NaN;
            end
        end
        foundosc_offsets_repeats_amp(:,rr) = oscdifferences; 
        
        % GET TRUE POS + FALSE POS etc. 
        [TP,FN,TN,FP, TPR, FPR, precision, recall] = get_OnOffTargetPeaks(x_sighits_amp, desiredperioddays, simsig.component_osc_periods_days);
        
        TPs_a(rr) = TP; % check this matrix filling, if it does so correctly
        FPs_a(rr) = FP;
        FNs_a(rr) = FN;
        TNs_a(rr) = TN;
        TPRs_a(rr) = TPR;
        FPRs_a(rr) = FPR;
        precisions_a(rr) = precision;
        recalls_a(rr) = recall;
        
    end 
    
    % Calculate percent significant hits per frequency -> with basic percentile
    xnumhits = sum(x_sighits_mat,2); 
    xpercenthits = (xnumhits./ size(x_sighits_mat,2)) * 100; 
    znumhits = sum(z_sighits_mat,2); 
    zpercenthits = (znumhits./ size(z_sighits_mat,1)) * 100; 

    % Calculate percent significant hits per frequency -> with additional amplitude threshold 
    xnumhits_a = sum(x_sighits_amp_mat,2); 
    xpercenthits_a = (xnumhits./ size(x_sighits_amp_mat,2)) * 100; % CONFIRM DIMENSIONS, DENOMINATOR SHOULD BE length 1000

    % Update simsig with results 
    simsig.foundosc_offsets_reshufflings = foundosc_offsets_repeats; % closest found peaks to the real/known oscillations 
    simsig.foundosc_offsets_reshufflings_amp = foundosc_offsets_repeats_amp; % closest found peaks to the real/known oscillations 

    simsig.x = x;
    simsig.xpercenthits = xpercenthits;
    simsig.xpercenthits_a = xpercenthits_a;
    simsig.x_sighits = x_sighits_mat;
    simsig.x_sighits_amp = x_sighits_amp_mat;
    simsig.x_reshufflings = x_reshufflings;
    simsig.x_percentiles = x_percentiles;

    simsig.z = z;
    simsig.zpercenthits = zpercenthits;
    simsig.z_sighits = z_sighits_mat;
    simsig.z_reshufflings = z_reshufflings;
    simsig.z_percentiles = z_percentiles;

    simsig.desiredperioddays = desiredperioddays;
    simsig.f_k = f_k;
    simsig.dct_basis = dct_basis;
    simsig.poly_basis = poly_basis;
    simsig.BPDN_params.delta = delta; 
    simsig.BPDN_params.pmin = minperiod;
    simsig.BPDN_params.pmax = maxperiod;
    simsig.BPDN_params.maxdegree = maxdegree;
    simsig.BPDN_params.n = N; 
    simsig.BPDN_params.dt = dt;

    simsig.TP = TPs;
    simsig.FP = FPs;
    simsig.TN = TNs;
    simsig.FN = FNs;
    simsig.TPR = TPRs;
    simsig.FPR = FPRs;
    simsig.precision = precisions;
    simsig.recalls = recalls;
    
    simsig.TP_amp = TPs_a;
    simsig.FP_amp = FPs_a;
    simsig.TN_amp = TNs_a;
    simsig.FN_amp = FNs_a;
    simsig.TPR_amp = TPRs_a;
    simsig.FPR_amp = FPRs_a;
    simsig.precision_amp = precisions_a;
    simsig.recalls_amp = recalls_a;
    
    % TO FINISH - SAVE UPDATED SIMSIG FILE
    save(loopfilename,'simsig'); 
    toc
    display(strcat('Finished sig-hits run for signal',{' '},num2str(i),{' '},'of',{' '},num2str(length(files))))
end

% CREATE FIGURES SHOWING MODEL PERFORMANCE --------------------------------

% Figures based on the significant hits simulations - frequency dependence
% each figure is y axis percent hits, x axis poi in figure. 
% with lots of params, most need to be fixed for a particular
% representation/question
% be careful with this code. there is LOTS of hardcoding in terms of
% indexing the different datasets. 

destfolder = 'P:\Personal\Irena\BPDN_methods_paper\simulated_signals_focusOnSparseness\'; 

% Vector of frequencies
phaseshift = 0;
fsize=15;
lw=3;
samplesperdays = [1/7, 4/7]; %[1/28,5/7]; % 1, 2, 7, 10
SNR_models = [1,5];
variances = [1,10];
component_periods = [1, 5, 30, 50, 80, 100, 120]; 
colors = {'blue','red'};

% FOCUS = DURATION --------------------------------------------------

% % Focus timeseries
% focusfeatstring = 'durations';
% durations = [3, 12, 36]; % will have a line for each.
% focusfeat = durations;
% signaltype = 'oscillation'; % alternative is white noise.... 
% % Fixed params:
% phaseshift = 1;
% samplingtype = 'irregular';
% samplesperday = 3/7; % middle of the road
% SNR_model = 5; % middle of the road
% variance = 10; % middle of the road
% % Vector of frequencies
% 
% toplot = zeros(length(durations),length(component_periods));
% 
% % loop through durations
% for d = 1: length(durations)
%     duration = durations(d);
%     % loop through period
%     for p = 1: length(component_periods)
%         poi = component_periods(p);
%         load(strcat(destfolder,signaltype,'Duration',num2str(duration),'months_phaseshift',num2str(phaseshift),'_samplingtype',samplingtype,'_sampperday',num2str(samplesperday),'_SNRmodel',num2str(SNR_model), '_variance',num2str(variance),'_osc',num2str(poi),'days.mat')); 
%         toplot(d,p) = simsig.xpercenthits; 
%     end 
% end 
% 
% figure
% for i = 1:length(size(toplot,1))
%     plot(component_periods, toplot(i,:),'DisplayName',num2str(focusfeat(i)))
%     xlabel('Period (days)')
%     ylabel('Percent detected')
%     title(focusfeatstring)
% end 

% FOCUS = Variance ------------------------------------------------
% This version has two subplots, one for high and one for low SNR
% Focus timeseries
focusfeatstring = 'variances';
focusfeat = variances;
% Fixed params:
signaltype = 'oscillation'; % alternative is white noise.... 
phaseshift = 0;
samplingtype = 'irregular';
duration = 12;% middle of the road

figure
for ss = 1:length(samplesperdays)
    subplot(length(samplesperdays),2,ss*2-1)
    sampstring = strcat(num2str(samplesperdays(ss)*7),{' '},'s per week');
    
    % LOW SNR SUBPLOT
    SNR_model = 1; %LOW SNR
    toplot = zeros(length(variances),length(component_periods));
    for d = 1: length(variances)
        variance = variances(d);
        % loop through period
        for p = 1: length(component_periods)
            poi = component_periods(p);
            load(strcat(destfolder,signaltype,'Duration',num2str(duration),'months_phaseshift',num2str(phaseshift),'_samplingtype',samplingtype,'_sampperday',num2str(samplesperdays(ss)),'_SNRmodel',num2str(SNR_model), '_variance',num2str(variance),'_osc',num2str(poi),'days.mat')); 
            toplot(d,p) =  (sum(simsig.TP)/length(simsig.TP))*100;%simsig.xpercenthits; % percent hits is weird variable. use True positives per iteration to get 
        end 
    end 
    for i = 1:size(toplot,1)
        plot(component_periods, toplot(i,:),'o-','DisplayName',strcat('var=',num2str(focusfeat(i))),'LineWidth',lw)
        xlabel('Period (days)')
        ylabel({sampstring{:}; 'Percent detected'})
        ylim([0, 110])
        xlim([-5, 130])
        set(gca,'FontSize',fsize)
        title(strcat('Low SNR (model based):',{' '},num2str(SNR_model)))
        hold on
    end 
    legend('Location','southwest')

    subplot(length(samplesperdays),2,ss*2)
    % HIGH SNR SUBPLOT
    SNR_model = 5;
    toplot = zeros(length(variances),length(component_periods));
    % loop through durations
    for d = 1: length(variances)
        variance = variances(d);
        % loop through period
        for p = 1: length(component_periods)
            poi = component_periods(p);
            load(strcat(destfolder,signaltype,'Duration',num2str(duration),'months_phaseshift',num2str(phaseshift),'_samplingtype',samplingtype,'_sampperday',num2str(samplesperdays(ss)),'_SNRmodel',num2str(SNR_model), '_variance',num2str(variance),'_osc',num2str(poi),'days.mat')); 
            toplot(d,p) =  (sum(simsig.TP)/length(simsig.TP))*100;%simsig.xpercenthits; % percent hits is weird variable. use True positives per iteration to get 
        end 
    end 
    for i = 1:size(toplot,1)
        plot(component_periods, toplot(i,:),'o-','DisplayName',strcat('var=',num2str(focusfeat(i))), 'LineWidth',lw)
        xlabel('Period (days)')
        ylabel({sampstring{:}; 'Percent detected'})
        ylim([0, 110])
        xlim([-5, 130])
        title(strcat('High SNR (model based):',{' '},num2str(SNR_model)))
        set(gca,'FontSize',fsize)
        hold on

    end 
    legend('Location','southwest')
end 
    

% FOCUS = Variance ------------------------------------------------
% This version has two subplots, one for high and one for low SNR
% Focus timeseries
focusfeatstring = 'variances';
focusfeat = variances;
% Fixed params:
signaltype = 'oscillation'; % alternative is white noise.... 
phaseshift = 0;
samplingtype = 'irregular';
duration = 12;% middle of the road

figure
for ss = 1:length(samplesperdays)
    subplot(length(samplesperdays),2,ss*2-1)
    sampstring = strcat(num2str(samplesperdays(ss)*7),{' '},'s per week');
    
    % LOW SNR SUBPLOT
    SNR_model = 1; %LOW SNR
    toplot = zeros(length(variances),length(component_periods));
    for d = 1: length(variances)
        variance = variances(d);
        % loop through period
        for p = 1: length(component_periods)
            poi = component_periods(p);
            load(strcat(destfolder,signaltype,'Duration',num2str(duration),'months_phaseshift',num2str(phaseshift),'_samplingtype',samplingtype,'_sampperday',num2str(samplesperdays(ss)),'_SNRmodel',num2str(SNR_model), '_variance',num2str(variance),'_osc',num2str(poi),'days.mat')); 
            toplot(d,p) =  (sum(simsig.TP)/length(simsig.TP))*100;%simsig.xpercenthits; % percent hits is weird variable. use True positives per iteration to get 
        end 
    end 
    for i = 1:size(toplot,1)
        plot(component_periods, toplot(i,:),'o-','DisplayName',strcat('var=',num2str(focusfeat(i))),'LineWidth',lw)
        xlabel('Period (days)')
        ylabel({sampstring{:}; 'Percent detected'})
        ylim([0, 110])
        xlim([-5, 130])
        set(gca,'FontSize',fsize)
        title(strcat('Low SNR (model based):',{' '},num2str(SNR_model)))
        hold on
    end 
    legend('Location','southwest')

    subplot(length(samplesperdays),2,ss*2)
    % HIGH SNR SUBPLOT
    SNR_model = 5;
    toplot = zeros(length(variances),length(component_periods));
    % loop through durations
    for d = 1: length(variances)
        variance = variances(d);
        % loop through period
        for p = 1: length(component_periods)
            poi = component_periods(p);
            load(strcat(destfolder,signaltype,'Duration',num2str(duration),'months_phaseshift',num2str(phaseshift),'_samplingtype',samplingtype,'_sampperday',num2str(samplesperdays(ss)),'_SNRmodel',num2str(SNR_model), '_variance',num2str(variance),'_osc',num2str(poi),'days.mat')); 
            toplot(d,p) =  (sum(simsig.TP)/length(simsig.TP))*100;%simsig.xpercenthits; % percent hits is weird variable. use True positives per iteration to get 
        end 
    end 
    for i = 1:size(toplot,1)
        plot(component_periods, toplot(i,:),'o-','DisplayName',strcat('var=',num2str(focusfeat(i))), 'LineWidth',lw)
        xlabel('Period (days)')
        ylabel({sampstring{:}; 'Percent detected'})
        ylim([0, 110])
        xlim([-5, 130])
        title(strcat('High SNR (model based):',{' '},num2str(SNR_model)))
        set(gca,'FontSize',fsize)
        hold on

    end 
    legend('Location','southwest')
end 

% FOCUS - SNR_model--WITH SIGNAL EXAMPLES TOO ----------------------------
% This version has two subplots, one for high and one for low variance
% Focus timeseries
samplesperdaysmini =[1/7, 4/7]; %[1/28,5/7]; % 1, 2, 7, 10
focusfeatstring = 'SNR_models';
focusfeat = SNR_models;
% Fixed params:
signaltype = 'oscillation'; % alternative is white noise.... 
samplingtype = 'irregular';
duration = 12;% middle of the road

% FOR THE LOW VAR CONDITION ------------------
fig = figure('Position',[10,10,900,900]);
suptitle('Low variance condition')
% GET LOW VAR DATA (PERCENT HITS) AND SIGNALS
ss=1;
variance = 1; 
for d = 1: length(SNR_models)
    SNR_model = SNR_models(d);
    % loop through period
    sigs = [];
    for p = 1: length(component_periods)
        poi = component_periods(p);
        load(strcat(destfolder,signaltype,'Duration',num2str(duration),'months_phaseshift',num2str(phaseshift),'_samplingtype',samplingtype,'_sampperday',num2str(samplesperdaysmini(ss)),'_SNRmodel',num2str(SNR_model), '_variance',num2str(variance),'_osc',num2str(poi),'days.mat')); 
        toplot(d,p) =  (sum(simsig.TP)/length(simsig.TP))*100;%simsig.xpercenthits; % percent hits is weird variable. use True positives per iteration to get 
        sigs(p,:) = simsig.signal;
    end 
    toplotsig{d} = sigs;
end 

% PLOT THE SOURCE SIGNALS
for i = 1: size(toplot,2)
    
    % SNR column 1 (low)
    subplot(7,8,[1 2 3 4] + ((i-1)*8))
    plot(simsig.t_full ./(60*60*14),toplotsig{1,1}(i,:),'Color','blue')
    xlim([0 600])
    ylim([-40 40])
    ylabel(strcat(num2str(component_periods(i)),{' '}, 'd. c.'))
    if i == 1
        title('Low SNR')
    end
    if i == size(toplot,2)
        xlabel('Time (days)')
    end 

    % SNR column 2 (high)
    subplot(7,8,[5,6,7,8] + ((i-1)*8))
    plot(simsig.t_full ./(60*60*14),toplotsig{1,2}(i,:),'Color','red')
    xlim([0 600])
    ylim([-40 40])
    if i == 1
        title('High SNR')
    end
    if i == size(toplot,2)
        xlabel('Time (days)')
    end 
end 
savename = strcat('H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\sim_exampsigs_lowvar.tiff');
exportgraphics(fig,savename,'Resolution',1200)

% PLOT THE PERFORMANCES
fig = figure('Position',[10,10,900,900]);
lw = 0.25;
% FOR LOW SAMPLING
ss=1;
variance=1;
for d = 1: length(SNR_models)
    SNR_model = SNR_models(d);
    % loop through period
    for p = 1: length(component_periods)
        poi = component_periods(p);
        load(strcat(destfolder,signaltype,'Duration',num2str(duration),'months_phaseshift',num2str(phaseshift),'_samplingtype',samplingtype,'_sampperday',num2str(samplesperdaysmini(ss)),'_SNRmodel',num2str(SNR_model), '_variance',num2str(variance),'_osc',num2str(poi),'days.mat')); 
        scaledmean = (sum(simsig.TP)/length(simsig.TP))*100;
        stde = std(simsig.TP); % TBD IF THIS INPUT IS CORRECT
        n = length(simsig.TP);
        toplot_lb(d,p) = scaledmean - 1.96*(stde / sqrt(n));
        toplot(d,p) = scaledmean;
        toplot_ub(d,p) = scaledmean + 1.96*(stde / sqrt(n));
    end 
end 
s = [9,10,11,12];
newvec = zeros(size(toplot,2),4);
for i = 1:size(toplot,2)
    newvec(i,:) = s + ((i-1)*16);
end 
newvec = reshape(newvec,[1,numel(newvec)]);
% create plot indices
% subplot(7,16,newvec)
subplot(1,2,1)

for i = 1:size(toplot,1)
    xconf = [component_periods component_periods(end:-1:1)] ;         
    yconf = [toplot_ub(i,:) fliplr(toplot_lb(i,:))];
    p = fill(xconf,yconf,colors{i},'HandleVisibility','off');
    p.FaceColor = colors{i};  
    p.FaceAlpha = 0.2;
    p.EdgeColor = 'none';   
    hold on
    if focusfeat(i) == min(focusfeat)
        labels = 'low';
    elseif focusfeat(i) == max(focusfeat)
        labels = 'high';
    end 
    labelstring = strcat('SNR =',{' '},labels);
    plot(component_periods, toplot(i,:),'o-','Color',colors{i},'DisplayName',labelstring{:},'LineWidth',lw)
    xlabel('Period (days)')
    ylabel({'';'';'Percent detected'})
    ylim([0, 110])
    xlim([-5, 130])
%     set(gca,'FontSize',fsize)
    title(strcat('Low variance:',{' '},num2str(variance)))
    hold on
end 
title(strcat('Samples per week =',{' '},num2str(samplesperdaysmini(ss)*7)));
legend('Location','southwest')

% FOR DENSER SAMPLING
ss=2;
variance = 1; 
for d = 1: length(SNR_models)
    SNR_model = SNR_models(d);
    % loop through period
    for p = 1: length(component_periods)
        poi = component_periods(p);
        load(strcat(destfolder,signaltype,'Duration',num2str(duration),'months_phaseshift',num2str(phaseshift),'_samplingtype',samplingtype,'_sampperday',num2str(samplesperdaysmini(ss)),'_SNRmodel',num2str(SNR_model), '_variance',num2str(variance),'_osc',num2str(poi),'days.mat')); 
        scaledmean = (sum(simsig.TP)/length(simsig.TP))*100;
        stde = std(simsig.TP); % TBD IF THIS INPUT IS CORRECT
        n = length(simsig.TP);
        toplot_lb(d,p) = scaledmean - 1.96*(stde / sqrt(n));
        toplot(d,p) = scaledmean;
        toplot_ub(d,p) = scaledmean + 1.96*(stde / sqrt(n));
    end 
end 

s = [13,14,15,16];
newvec = zeros(size(toplot,2),4);
for i = 1:size(toplot,2)
    newvec(i,:) = s + ((i-1)*16);
end 
newvec = reshape(newvec,[1,numel(newvec)]);

% create plot indices
% subplot(7,16,newvec)
subplot(1,2,2)
for i = 1:size(toplot,1)
    xconf = [component_periods component_periods(end:-1:1)] ;         
    yconf = [toplot_ub(i,:) fliplr(toplot_lb(i,:))];
    p = fill(xconf,yconf,colors{i},'HandleVisibility','off');
    p.FaceColor = colors{i};  
    p.FaceAlpha = 0.2;
    p.EdgeColor = 'none';   
    hold on
    if focusfeat(i) == min(focusfeat)
        labels = 'low';
    elseif focusfeat(i) == max(focusfeat)
        labels = 'high';
    end 
    labelstring = strcat('SNR =',{' '},labels);
    plot(component_periods, toplot(i,:),'o-','Color',colors{i},'DisplayName',labelstring{:},'LineWidth',lw)
    xlabel('Period (days)')
    ylabel({'';'';'Percent detected'})
    ylim([0, 110])
    xlim([-5, 130])
%     set(gca,'FontSize',fsize)
    title(strcat('Low variance:',{' '},num2str(variance)))
    hold on
end 
title(strcat('Samples per week =',{' '},num2str(samplesperdaysmini(ss)*7)));
legend('Location','southwest')
savename = strcat('H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\sim_performance_lowvar.tiff');
exportgraphics(fig,savename,'Resolution',1200)

% FOR THE HIGH VAR CONDITION _--------------------------------
fig = figure('Position',[10,10,1000,900]);
suptitle('High variance condition')
% GET LOW VAR DATA (PERCENT HITS) AND SIGNALS
ss=1;
variance = 10; 
for d = 1: length(SNR_models)
    SNR_model = SNR_models(d);
    % loop through period
    sigs = [];
    for p = 1: length(component_periods)
        poi = component_periods(p);
        load(strcat(destfolder,signaltype,'Duration',num2str(duration),'months_phaseshift',num2str(phaseshift),'_samplingtype',samplingtype,'_sampperday',num2str(samplesperdaysmini(ss)),'_SNRmodel',num2str(SNR_model), '_variance',num2str(variance),'_osc',num2str(poi),'days.mat')); 
        scaledmean = (sum(simsig.TP)/length(simsig.TP))*100;
        stde = std(simsig.TP); % TBD IF THIS INPUT IS CORRECT
        n = length(simsig.TP);
        toplot_lb(d,p) = scaledmean - 1.96*(stde / sqrt(n));
        toplot(d,p) = scaledmean;
        toplot_ub(d,p) = scaledmean + 1.96*(stde / sqrt(n));        sigs(p,:) = simsig.signal;
    end 
    toplotsig{d} = sigs;
end 

% PLOT THE SOURCE SIGNALS
for i = 1: size(toplot,2)
    
    % SNR column 1 (low)
    subplot(7,8,[1 2 3 4] + ((i-1)*8))
    plot(simsig.t_full ./(60*60*14),toplotsig{1,1}(i,:),'Color','blue')
    xlim([0 600])
    ylim([-40 40])
    ylabel(strcat(num2str(component_periods(i)),{' '}, 'd. c.'))
    if i == 1
        title('Low SNR')
    end
    if i == size(toplot,2)
        xlabel('Time (days)')
    end 

    % SNR column 2 (high)
    subplot(7,8,[5,6,7,8] + ((i-1)*8))
    plot(simsig.t_full ./(60*60*14),toplotsig{1,2}(i,:),'Color','red')
    xlim([0 600])
    ylim([-40 40])
    if i == 1
        title('High SNR')
    end
    if i == size(toplot,2)
        xlabel('Time (days)')
    end 
end 
savename = strcat('H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\sim_exampsigs_highvar.tiff');
exportgraphics(fig,savename,'Resolution',1200)


% PLOT THE PERFORMANCES
% FOR LOW SAMPLING
% create plot indices
fig = figure('Position',[10,10,900,900]);
lw = 0.25;
s = [9,10,11,12];
newvec = zeros(size(toplot,2),4);
for i = 1:size(toplot,2)
    newvec(i,:) = s + ((i-1)*16);
end 
newvec = reshape(newvec,[1,numel(newvec)]);

% subplot(7,16,newvec)
subplot(1,2,1)
for i = 1:size(toplot,1)
    xconf = [component_periods component_periods(end:-1:1)] ;         
    yconf = [toplot_ub(i,:) fliplr(toplot_lb(i,:))];
    p = fill(xconf,yconf,colors{i},'HandleVisibility','off');
    p.FaceColor = colors{i};  
    p.FaceAlpha = 0.2;
    p.EdgeColor = 'none';   
    hold on
%     plot(component_periods, toplot(i,:),'o-','Color',colors{i},'DisplayName',strcat('SNR=',num2str(focusfeat(i))),'LineWidth',lw)
%         xlabel('Period (days)')
    if focusfeat(i) == min(focusfeat)
        labels = 'low';
    elseif focusfeat(i) == max(focusfeat)
        labels = 'high';
    end 
    labelstring = strcat('SNR =',{' '},labels);
    plot(component_periods, toplot(i,:),'o-','Color',colors{i},'DisplayName',labelstring{:},'LineWidth',lw)
    xlabel('Period (days)')
    ylabel({'';'';'Percent detected'})
    ylim([0, 110])
    xlim([-5, 130])
%     set(gca,'FontSize',fsize)
    title(strcat('Low variance:',{' '},num2str(variance)))
    hold on
end 
title(strcat('Samples per week =',{' '},num2str(samplesperdaysmini(ss)*7)));
legend('Location','southwest')

% FOR DENSER SAMPLING
ss=2;
variance = 10; 
for d = 1: length(SNR_models)
    SNR_model = SNR_models(d);
    % loop through period
    for p = 1: length(component_periods)
        poi = component_periods(p);
        load(strcat(destfolder,signaltype,'Duration',num2str(duration),'months_phaseshift',num2str(phaseshift),'_samplingtype',samplingtype,'_sampperday',num2str(samplesperdaysmini(ss)),'_SNRmodel',num2str(SNR_model), '_variance',num2str(variance),'_osc',num2str(poi),'days.mat')); 
        scaledmean = (sum(simsig.TP)/length(simsig.TP))*100;
        stde = std(simsig.TP); % TBD IF THIS INPUT IS CORRECT
        n = length(simsig.TP);
        toplot_lb(d,p) = scaledmean - 1.96*(stde / sqrt(n));
        toplot(d,p) = scaledmean;
        toplot_ub(d,p) = scaledmean + 1.96*(stde / sqrt(n));
    end 
end 
% create plot indices 
s = [13,14,15,16];
newvec = zeros(size(toplot,2),4);
for i = 1:size(toplot,2)
    newvec(i,:) = s + ((i-1)*16);
end 
newvec = reshape(newvec,[1,numel(newvec)]);
% 
% subplot(7,16,newvec)
subplot(1,2,2)
for i = 1:size(toplot,1)
    xconf = [component_periods component_periods(end:-1:1)] ;         
    yconf = [toplot_ub(i,:) fliplr(toplot_lb(i,:))];
    p = fill(xconf,yconf,colors{i},'HandleVisibility','off');
    p.FaceColor = colors{i};  
    p.FaceAlpha = 0.2;
    p.EdgeColor = 'none';   
    hold on
%     plot(component_periods, toplot(i,:),'o-','Color',colors{i},'DisplayName',strcat('SNR=',num2str(focusfeat(i))),'LineWidth',lw)
    if focusfeat(i) == min(focusfeat)
        labels = 'low';
    elseif focusfeat(i) == max(focusfeat)
        labels = 'high';
    end
    labelstring = strcat('SNR =',{' '},labels);
    plot(component_periods, toplot(i,:),'o-','Color',colors{i},'DisplayName',labelstring{:},'LineWidth',lw)
    xlabel('Period (days)')
    ylabel({'';'';'Percent detected'})
    ylim([0, 110])
    xlim([-5, 130])
%     set(gca,'FontSize',fsize)
    title(strcat('Low variance:',{' '},num2str(variance)))
    hold on
end 
title(strcat('Samples per week =',{' '},num2str(samplesperdaysmini(ss)*7)));
legend('Location','southwest')

savename = strcat('H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\sim_performance_highvar.tiff');
exportgraphics(fig,savename,'Resolution',1200)

%% Supplemental figure 3
% See Figure 4

%% Supplemental figure 4
% See Figure 4

%% Supplemental figure 5
% See Figure 5

%% Supplemental figure 6
% See Figure 5

%% Supplemental figure 7
% See Figure 6

%% Supplemental figure 8
% See Figure 6

%% Supplemental figure 9
% See Supp. Figure 2 (re-ran parameter sweep with 90/10 cross val)
% ie percent2drop = 10; 

%% Supplemental figure 10
% See Figure 7

%% Supplemental figure 11 
% See Figure 7

%% Supplemental figure 12
% spectral outputs for 100x repeats vs 1000x repeats
% Spectra with sig stars + noise floor inserts for S1 and all real subjects
% just change subject name and run to get one for each

subject = 'M1';
includepoly = 0;
fsize = 10;

% PART 1: BASIC RAW DATA AND CWT SPECTROGRAM 
duration = 12;
N = 3000;
numsamp = 1800;
if strcmp(subject,'S1') == 1
    destfolder = strcat('P:\Personal\Irena\BPDN_methods_paper\simulateddata\',subject,'\modeloutput\');
else
    destfolder = strcat('P:\Personal\Irena\BPDN_methods_paper\realdata\',subject,'\modeloutput\');
end
filename = strcat(destfolder,'Outputfor',num2str(duration),'_N',num2str(N),'_numsamp',num2str(numsamp),'_traintest_forfig.mat'); % original

load(filename);

% Initialize variables
desiredperioddays = deltaselection.desiredperioddays;
x = deltaselection.fig.x;
sighits = deltaselection.fig.sighits;
perioddays_cwt = deltaselection.fig.perioddays_cwt;
meanpowerCWTreg = deltaselection.fig.meanpowerCWTreg;

fig = figure('Position',[10 10 1000 300]);
% --- FOR ORIGINAL 100 X ------------------------------
% overlay BPDN and CWT spectra
subplot(1,2,1)
% spectrum for BPDN
yyaxis right
% no rescaling
plot(desiredperioddays, x.^2, 'Color',[0.8500 0.3250 0.0980], 'LineWidth',1);
hold on
plot(desiredperioddays(sighits),x(sighits).^2,'*','Color','black','MarkerSize',4)
% xlabel('Period (days)')
ylabel('Power, BPDN')
set(gca,'FontSize',fsize)
% spectrum for cwt
yyaxis left
plot(perioddays_cwt, meanpowerCWTreg, 'Color',[0 0 0 0.3], 'LineWidth',1.5); % blue color with 0.6 alpha
% hold on
% plot(perioddays_cwt(locs),pks,'Color','cyan','o', 'MarkerSize',2)
xlabel('Period (days)')
ylabel({'';'';'Power,CWT'})
title('100x')
ax=gca; ax.YAxis(1).Exponent = 5; ax.YAxis(1).Color = 'black';
xlim([-10 max(perioddays_cwt)])
set(gca,'FontSize',fsize)

% % plot the noise floor ------ Noise floor insert 
% figmini = figure('Position',[10,10,400, 200]);
% for iii = 1:size(deltaselection.fig.xreshuffled,1)
%     plot(desiredperioddays,deltaselection.fig.xreshuffled(iii,:).^2,'-','Color',[0.4940 0.1840 0.5560 0.2])
%     hold on
% end 
% plot(desiredperioddays, x.^2,'-', 'Color',[0.8500 0.3250 0.0980], 'LineWidth',1);
% hold on
% plot(desiredperioddays(sighits),x(sighits).^2,'*','Color','black','MarkerSize',4)
% set(gca,'YAxisLocation','right','ycolor',[0.8500 0.3250 0.0980])
% xlim([0 80])
% ylim([0 10*10^6])
% savename = strcat('H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\',subject,'_noisefloorinsert_100x.tiff');
% exportgraphics(figmini,savename,'Resolution',1200)

% ------------ FOR 1000x -------------------------------------
filename = strcat(destfolder,'Outputfor',num2str(duration),'_N',num2str(N),'_numsamp',num2str(numsamp),'_traintest_forfig_1000x_test.mat'); % 1000x comparison
load(filename);

% Initialize variables
desiredperioddays = deltaselection.desiredperioddays;
x = deltaselection.fig.x;
sighits = deltaselection.fig.sighits;
perioddays_cwt = deltaselection.fig.perioddays_cwt;
meanpowerCWTreg = deltaselection.fig.meanpowerCWTreg;

% overlay BPDN and CWT spectra
subplot(1,2,2)
% spectrum for BPDN
yyaxis right
% no rescaling
plot(desiredperioddays, x.^2, 'Color',[0.8500 0.3250 0.0980], 'LineWidth',1);
hold on
plot(desiredperioddays(sighits),x(sighits).^2,'*','Color','black','MarkerSize',4)
% xlabel('Period (days)')
ylabel('Power, BPDN')
set(gca,'FontSize',fsize)
% spectrum for cwt
yyaxis left
plot(perioddays_cwt, meanpowerCWTreg, 'Color',[0 0 0 0.3], 'LineWidth',1.5); % blue color with 0.6 alpha
% hold on
% plot(perioddays_cwt(locs),pks,'Color','cyan','o', 'MarkerSize',2)
xlabel('Period (days)')
ylabel({'';'';'Power,CWT'})
title('1000x')
ax=gca; ax.YAxis(1).Exponent = 5; ax.YAxis(1).Color = 'black';
xlim([-10 max(perioddays_cwt)])
set(gca,'FontSize',fsize)

% % plot the noise floor ------ Noise floor insert 
% figmini = figure('Position',[10,10,400, 200]);
% for iii = 1:size(deltaselection.fig.xreshuffled,1)
%     plot(desiredperioddays,deltaselection.fig.xreshuffled(iii,:).^2,'-','Color',[0.4940 0.1840 0.5560 0.2])
%     hold on
% end 
% plot(desiredperioddays, x.^2,'-', 'Color',[0.8500 0.3250 0.0980], 'LineWidth',1);
% hold on
% plot(desiredperioddays(sighits),x(sighits).^2,'*','Color','black','MarkerSize',4)
% set(gca,'YAxisLocation','right','ycolor',[0.8500 0.3250 0.0980])
% xlim([0 80])
% ylim([0 10*10^6])
% savename = strcat('H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\',subject,'_noisefloorinsert_1000x.tiff');
% exportgraphics(figmini,savename,'Resolution',1200)

savename = strcat('H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\',subject,'_100xVs1000xSpectra.tiff');
exportgraphics(fig,savename,'Resolution',600)


%% Supplemental figure 13
% spectral outputs for 100x repeats vs 10,000x repeats
% Spectra with sig stars + noise floor inserts for S1 and all real subjects
% just change subject name and run to get one for each

subject = 'S1';
includepoly = 0;
fsize = 10;

% PART 1: BASIC RAW DATA AND CWT SPECTROGRAM 
duration = 12;
N = 3000;
numsamp = 1800;
if strcmp(subject,'S1') == 1
    destfolder = strcat('P:\Personal\Irena\BPDN_methods_paper\simulateddata\',subject,'\modeloutput\');
else
    destfolder = strcat('P:\Personal\Irena\BPDN_methods_paper\realdata\',subject,'\modeloutput\');
end
filename = strcat(destfolder,'Outputfor',num2str(duration),'_N',num2str(N),'_numsamp',num2str(numsamp),'_traintest_forfig.mat'); % original

load(filename);

% Initialize variables
desiredperioddays = deltaselection.desiredperioddays;
x = deltaselection.fig.x;
sighits = deltaselection.fig.sighits;
perioddays_cwt = deltaselection.fig.perioddays_cwt;
periodcap = max(perioddays_cwt);
meanpowerCWTreg = deltaselection.fig.meanpowerCWTreg;

fig = figure('Position',[10 10 1000 300]);
% --- FOR ORIGINAL 100 X ------------------------------
% overlay BPDN and CWT spectra
subplot(1,2,1)
% spectrum for BPDN
yyaxis right
% no rescaling
plot(desiredperioddays, x.^2, 'Color',[0.8500 0.3250 0.0980], 'LineWidth',1);
hold on
plot(desiredperioddays(sighits),x(sighits).^2,'*','Color','black','MarkerSize',4)
% xlabel('Period (days)')
ylabel('Power, BPDN')
set(gca,'FontSize',fsize)
% spectrum for cwt
yyaxis left
plot(perioddays_cwt, meanpowerCWTreg, 'Color',[0 0 0 0.3], 'LineWidth',1.5); % blue color with 0.6 alpha
% hold on
% plot(perioddays_cwt(locs),pks,'Color','cyan','o', 'MarkerSize',2)
xlabel('Period (days)')
ylabel({'';'';'Power,CWT'})
title('100x')
ax=gca; ax.YAxis(1).Exponent = 5; ax.YAxis(1).Color = 'black';
xlim([-10 periodcap])
set(gca,'FontSize',fsize)

% % plot the noise floor ------ Noise floor insert 
% figmini = figure('Position',[10,10,400, 200]);
% for iii = 1:size(deltaselection.fig.xreshuffled,1)
%     plot(desiredperioddays,deltaselection.fig.xreshuffled(iii,:).^2,'-','Color',[0.4940 0.1840 0.5560 0.2])
%     hold on
% end 
% plot(desiredperioddays, x.^2,'-', 'Color',[0.8500 0.3250 0.0980], 'LineWidth',1);
% hold on
% plot(desiredperioddays(sighits),x(sighits).^2,'*','Color','black','MarkerSize',4)
% set(gca,'YAxisLocation','right','ycolor',[0.8500 0.3250 0.0980])
% xlim([0 80])
% ylim([0 10*10^6])
% savename = strcat('H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\',subject,'_noisefloorinsert_100x.tiff');
% exportgraphics(figmini,savename,'Resolution',1200)

% ------------ FOR 10000x -------------------------------------
destfolder = strcat('Y:\Irena\BPWP\realdata\',subject,'\modeloutput\'); 
filename = strcat(destfolder,'Outputfor',num2str(duration),'_N',num2str(N),'_numsamp',num2str(numsamp),'_traintest_forfig_10000x.mat'); % 1000x comparison
load(filename);

% Initialize variables
desiredperioddays = deltaselection.desiredperioddays;
x = deltaselection.fig.x(:,1);
sighits = deltaselection.fig.sighits;
perioddays_cwt = deltaselection.fig.perioddays_cwt;
meanpowerCWTreg = deltaselection.fig.meanpowerCWTreg;

% overlay BPDN and CWT spectra
subplot(1,2,2)
% spectrum for BPDN
yyaxis right
% no rescaling
plot(desiredperioddays, x.^2, 'Color',[0.8500 0.3250 0.0980], 'LineWidth',1);
hold on
plot(desiredperioddays(sighits),x(sighits).^2,'*','Color','black','MarkerSize',4)
% xlabel('Period (days)')
ylabel('Power, BPDN')
set(gca,'FontSize',fsize)
% spectrum for cwt
yyaxis left
plot(perioddays_cwt, meanpowerCWTreg, 'Color',[0 0 0 0.3], 'LineWidth',1.5); % blue color with 0.6 alpha
% hold on
% plot(perioddays_cwt(locs),pks,'Color','cyan','o', 'MarkerSize',2)
xlabel('Period (days)')
ylabel({'';'';'Power,CWT'})
title('10000x')
ax=gca; ax.YAxis(1).Exponent = 5; ax.YAxis(1).Color = 'black';
xlim([-10 periodcap])
set(gca,'FontSize',fsize)

% % plot the noise floor ------ Noise floor insert 
% figmini = figure('Position',[10,10,400, 200]);
% for iii = 1:size(deltaselection.fig.xreshuffled,1)
%     plot(desiredperioddays,deltaselection.fig.xreshuffled(iii,:).^2,'-','Color',[0.4940 0.1840 0.5560 0.2])
%     hold on
% end 
% plot(desiredperioddays, x.^2,'-', 'Color',[0.8500 0.3250 0.0980], 'LineWidth',1);
% hold on
% plot(desiredperioddays(sighits),x(sighits).^2,'*','Color','black','MarkerSize',4)
% set(gca,'YAxisLocation','right','ycolor',[0.8500 0.3250 0.0980])
% xlim([0 80])
% ylim([0 10*10^6])
% savename = strcat('H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\',subject,'_noisefloorinsert_1000x.tiff');
% exportgraphics(figmini,savename,'Resolution',1200)

savename = strcat('H:\Research\Worrell\Manuscripts\BPDN_Methods_Paper\figures\',subject,'_100xVs10000xSpectra.tiff');
exportgraphics(fig,savename,'Resolution',600)


