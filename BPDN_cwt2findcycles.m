%% DESCRIPTION: Function to run CWT (on continuous data) and identify cycle peaks

function [perioddays_cwt, meanpowerCWTreg, peakvalues, peaklocs, peakperiods, reconsigCWT] = cwt2findcycles(signal, dt_ogdata)
    
    [meanpowerCWTL2norm, perioddays_cwt, WT, fwt] = wavelet_decomp_L2norm(signal, 1/dt_ogdata);
    meanpowerCWTreg = nanmean(abs(WT).^2,2);
    reconsigCWT = icwt(WT);
    
    [peakvalues, peaklocs] = findpeaks(meanpowerCWTreg);
    peakperiods = perioddays_cwt(peaklocs);
    
end 
