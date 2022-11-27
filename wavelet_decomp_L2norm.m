% DESCRIPTION: Function that runs wavelet transform and calculate L2 norm

function [meanpowerCWTL2norm, perioddays_cwt, WT, fwt] = wavelet_decomp_L2norm(spikeseries, fsspikes)

    [WT, fwt, coi] = cwt(spikeseries, 'amor', fsspikes);
    wave_fam = 1./fwt; % wavelet family defined by period duration (units of days)
    
    % % % % L1 may be more typical; matter of preference L1 vs. L2 norm
    for i = 1:size(WT,1) % Change from L1 to L2 norm for cwt. L2 norm reduces amplitude of high freq data, while L1 does not. 
                         % Scale is inversely related to frequency. 
        WT2(i,:) = WT(i,:).*1/sqrt(fwt(i));  
    end
    
    power = abs(WT2);
    % get indices in cwt for boundaries of COI
    plotcone = 1 ./ coi; % cone on plot
    coneid = zeros(1,length(plotcone));
    
    for i = 1:length(plotcone)
        [d, ix] = min(abs(wave_fam - plotcone(i)));
        coneid(1,i) = ix;
    end 
    
    % Replace everything beyond coi boundary with nans
    for i = 1:length(coneid)
        power(coneid(i):end,i) = nan;
    end
    
    meanpowerCWTL2norm = nanmean(power,2);
    perioddays_cwt = 1./(fwt*24*60*60);
    
end 
