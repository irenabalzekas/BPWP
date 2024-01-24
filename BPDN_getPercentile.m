

% INPUT DATA
% distribs = distribution for each value of x (xreshuffled for example).
    % Each column is a unique instance (ie the a particular period - with each
    % row being the x calculated for that particular period with each
    % reshuffling)
% ref = values for which you're getting the percentile. (x for example)

% NOTE - double check that the dimensions along which percentiles are
% calculated suit your input data 

function [percentiles] = BPDN_getPercentile(distribs,ref)

   percentiles = [];
    for i = 1: size(distribs,2)
        historicalData = abs(distribs(:,i)); % NOTE changed to positive values 
        exogenousVariable = abs(ref(i));
        % Compute centile
        nless = sum(historicalData < exogenousVariable);
        nequal = sum(historicalData == exogenousVariable);
        percentiles(i) = 100 * (nless + 0.5*nequal) / length(historicalData);
    end 

end 