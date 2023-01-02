%% DESCRIPTION: Run delta parameter sweep on simulated data 
% Also plot outputs, update data with delta
% Note: Runs CVX

% DEFINE FIXED PARAMS PER SUBJECT
loc = [];
dest = [];
N = 300; 
numrepeats = 1; % number of times to recalculate MSE
percent2drop = 25; % percent of data allocated to testing dataset ie 75/25 crossvalidation
maxdegree = 3;
maxperiod = 120;
minperiod = 10; 
sampletype = 'densesubrange';
subjects = {'S1'};
sampperdays = [1,3,5];
deltas = 0:50000:1000000;

for p = 1:length(subjects)  
    
    subject = subjects{p};

    % LOAD + INITIALIZE DATA
    filename = strcat(loc,subject,'_simdata.mat');
    load(filename);  
    monthsdata = simdata.years*12;
    signal = simdata.spikes';
    time = simdata.time'; 
    
    % Adjust time representation to suit desired N for basis dimensions
    dt = round((max(time) - min(time))) / (N-1); % this is different from the original dt 
    t_full = 1:N;
    
    % Building the dct basis can be slow if N is large, good to do outside loop
    [f_k, desiredperioddays] = frequency_sampling(N,maxperiod, minperiod, dt, sampletype);
    [dct_basis, scalefactor] = DCT2_basis(N, f_k); 
    % Making the basis takes a long time so want to pre-define it before running BPDN for loops and loops 
    [poly_basis] = polynomial_basis(N, maxdegree);

    MSEpersamp = zeros(length(sampperdays)); 
    Xspersamp = zeros(length(sampperdays)); 
    zspersamp = zeros(length(sampperdays)); 
    measureddatapersamp = zeros(length(sampperdays)); 
    tpersamp = zeros(length(sampperdays));
    phipersamp = zeros(length(sampperdays)); 
    numsamps = zeros(length(sampperdays));
    training_md_persamp = zeros(length(sampperdays)); 
    training_t_persamp = zeros(length(sampperdays));
    training_phi_persamp = zeros(length(sampperdays)); 
    testing_md_persamp = zeros(length(sampperdays)); 
    testing_t_persamp = zeros(length(sampperdays)); 
    testing_phi_persamp = zeros(length(sampperdays));
    
    for s = 1:length(sampperdays)
        
        % Define numsamp needed 
        numsamp = monthsdata * 30 * sampperdays(s);
        numsamps(s) = numsamp;

        [measureddata_init, t_init, phi_init, training_md, training_t, training_phi, testing_md, testing_t, testing_phi] = BPDN_samplerealdata_traintest(signal,time,numsamp, dt, percent2drop, numrepeats);
       
        % Initialize variables 
        MSE = zeros(length(deltas), numrepeats);
        xs = zeros(length(deltas),numrepeats,N);
        zs = zeros(length(deltas),numrepeats, maxdegree+1);
        
        for dd = 1: length(deltas)

            delta = deltas(dd);

            parfor nn = 1: numrepeats

                % Training data
                measureddata = training_md(nn,:);
                t = training_t(nn,:);
                phi = training_phi(nn,:);

                % Testing data
                measured_missed = testing_md(nn,:);
                t_missed = testing_t(nn,:);
                phi_missed = testing_phi(nn,:);

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
    deltaselection.time_days = simdata.time_days;
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

    savename = strcat(dest,subject,'/deltaselection/MSE_2D_PercentDrop',num2str(percent2drop),'_repeats',num2str(numrepeats),'_N',num2str(N),'_numsamp',num2str(numsamps),'.mat');
    save(savename,'deltaselection');
    
end 

%% Plot outputs of parameter sweeps

subjects = {'S1'};

colors = {'red','black','blue','cyan','green','magenta'};
fsize = 15;

figure
for i = 1:length(subjects)
    
    subplot(1,length(subjects),i)
    subject = subjects{i};
    filename =  strcat(dest,subject,'/deltaselection/MSE_2D_PercentDrop',num2str(percent2drop),'_repeats',num2str(numrepeats),'_N',num2str(N),'_numsamp',num2str(numsamps),'.mat');
    load(filename);
    
    for nn = 1: length(deltaselection.sampperdays)
    
        x = deltaselection.deltas;
        dat = deltaselection.MSE{nn}'; % Number of ‘Experiments’ In Data Set
        yMean = mean(dat,1);  % Mean Of All Experiments At Each Value Of ‘x’
        ySEM = std(dat,0,1)/sqrt(size(dat,1));   % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
        CI95 = tinv([0.025 0.975], (size(dat,1)-1)); % Calculate 95% Probability Intervals Of t-Distribution
        yCI95 = bsxfun(@times, ySEM, CI95(:)); % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’
    
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
    title(strcat(subject,{' '},' cross validation'))
    legend
    set(gca,'FontSize',fsize)
    
end 

%% Update data with preferred value for delta

subject = [];
filename = strcat(dest,subject,'/deltaselection/MSE_2D_PercentDrop',num2str(percent2drop),'_repeats',num2str(numrepeats),'_N',num2str(N),'_numsamp',num2str(numsamps),'.mat');
load(filename);

deltaselection.delta = []; 

save(filename,'deltaselection');
