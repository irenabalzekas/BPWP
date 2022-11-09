
% drop - percent to drop

function [measureddata_init, t_init, phi_init, training_md, training_t, training_phi, testing_md, testing_t, testing_phi] = BPDN_samplerealdata_traintest(signal,time,numsamp, dt, drop, numrepeats)
    
    % CREATE REFERENCE/BASE SIG
    numsamp_full = round(numsamp * (1+drop/100));
    
    totind = sort(datasample(1:length(signal), numsamp_full, 'Replace',false)); % randomly sample the full seconds timeseries
    measureddata_init = signal(totind);
    t_init = time(totind);
    phi_init = round((t_init-min(time))/dt) + 1; % indices for when each time data were sampled
    numsampnew = numsamp_full;
    while length(unique(phi_init)) < numsamp_full
        length(unique(phi_init))
        numsampnew = numsampnew + round(numsamp_full / 20);
        if numsampnew > numsamp_full % if its more than all the samples, just sample all of them 
            numsampnew = numsamp_full;
        end
        totind = sort(datasample(1:length(signal), numsampnew, 'Replace',false)); % randomly sample the full seconds timeseries
%         end 
        measureddata_init = signal(totind);
        t_init = time(totind);
        phi_init = round((t_init-min(time))/dt) + 1;  
    end 
    % once long enough, pull unqiue and cut excess.
    [~,indlocs] = ismember(unique(phi_init), phi_init);
    % remove excess samples randomly
    tocut = sort(datasample(1:length(indlocs), length(indlocs) - numsamp_full, 'Replace',false)); % randomly sample the full seconds timeseries
    indlocs(tocut) = [];
    
    % BASE SIG DEFINED
    measureddata_init = measureddata_init(indlocs);
    t_init = t_init(indlocs);
    phi_init = phi_init(indlocs);
    
    % NOW RANDOMLY SUBSAMPLE THE BASE SIG. 
    training_md = [];
    training_t = [];
    training_phi = [];
    testing_md = [];
    testing_t = [];
    testing_phi = [];
    
    for i = 1:numrepeats
        
        % Randomly sample subset of data 
        trainind = sort(datasample(1:length(measureddata_init), numsamp, 'Replace',false)); % randomly sample the full seconds timeseries
        testind =  setdiff(find(1:length(measureddata_init)),trainind);

        % Define training and testing blocks 
        training_md(i,:) = measureddata_init(trainind);
        training_t(i,:) = t_init(trainind);
        training_phi(i,:) = phi_init(trainind);
        
        testing_md(i,:) = measureddata_init(testind);
        testing_t(i,:) = t_init(testind);
        testing_phi(i,:) = phi_init(testind);

    end 

    
end