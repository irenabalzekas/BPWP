% IB 7/22/22

% INPUT
% - regular BPDN function inputs

% OUTPUT
% - BPDN output for the original signal
% - BPDN outputs (distribution) from re-shuffled original signal
% - associated percentiles of original vs reshuffled

% NOTE
% - short because only reshuffles 100 times instead of 1000

function [x, xreshuffled, xpercentiles, z, zreshuffled, zpercentiles, reshuff_optval] = BPDN_wReshuffling_short(delta, measureddata, phi, dct_basis, poly_basis)

    % Calculate BPDN for real data 
    [x, z, ~, ~, ~, cvx_optval] = BPDN(delta, measureddata, phi, dct_basis, poly_basis);
    
    % Calculate BPDN for reshuffled data 
    xreshuffled = [];
    reshuff_optval = [];
    zreshuffled = [];
    zreshuff_optval = [];
    parfor i = 1:100 % CHANGE BACK
        ynew = measureddata(randperm(length(measureddata))); % shuffled order of measureddata (though sampling matrix (ie temporal spacing) is unchanged)
        % SOLVE
        [xres, zres, ~, ~, ~, cvx_optval] = BPDN(delta, ynew, phi, dct_basis, poly_basis);
        xreshuffled(i,:) = xres; % fill array where each row is another X output/solution
        zreshuffled(i,:) = zres;
%         reshuff_optval{i} = cvx_optval;
    end 

    % get percentile status for each output in x from original data
    xpercentiles = [];
    for i = 1: size(xreshuffled,2)
        historicalData = abs(xreshuffled(:,i)); % NOTE changed to positive values 
        exogenousVariable = abs(x(i));
        % Compute centile
        nless = sum(historicalData < exogenousVariable);
        nequal = sum(historicalData == exogenousVariable);
        xpercentiles(i) = 100 * (nless + 0.5*nequal) / length(historicalData);
    end 
    
   % get percentile status for each output in z from original data
    zpercentiles = [];
    for i = 1: size(zreshuffled,2)
        historicalData = abs(zreshuffled(:,i)); % NOTE changed to positive values 
        exogenousVariable = abs(z(i));
        % Compute centile
        nless = sum(historicalData < exogenousVariable);
        nequal = sum(historicalData == exogenousVariable);
        zpercentiles(i) = 100 * (nless + 0.5*nequal) / length(historicalData);
    end 
    %[~,c] = find(percentiles >=99);

end 

%% OLD VERSION

% function [x, xreshuffled,xpercentiles, z, zreshuffled, zpercentiles, f_adj, dt, f_final] = BPDN_wReshuffling(delta, measureddata, time, pmin, pmax, maxdegree)
% 
%     % Calculate BPDN for real data 
%     [x, z, dt, ~, ~, ~, f_adj, f_final, ~, ~, ~] = BPDN(delta, measureddata, time, pmin, pmax, maxdegree);
% 
%     % Calculate BPDN for reshuffled data 
%     xreshuffled = [];
%     zreshuffled = [];
%     for i = 1:1000 % for 1000 iterations
%         ynew = measureddata(randperm(length(measureddata))); % shuffled order of measureddata (though sampling matrix (ie temporal spacing) is unchanged)
%         % SOLVE
%         [xres, zres, ~, ~, ~, ~, ~, ~, ~, ~, ~] = BPDN(delta, ynew, time, pmin, pmax, maxdegree);
%         xreshuffled(i,:) = xres; % fill array where each row is another X output/solution
%         zreshuffled(i,:) = zres;
%     end 
% 
%     % get percentile status for each output in x from original data
%     xpercentiles = [];
%     for i = 1: size(xreshuffled,2)
%         historicalData = abs(xreshuffled(:,i)); % NOTE changed to positive values 
%         exogenousVariable = abs(x(i));
%         % Compute centile
%         nless = sum(historicalData < exogenousVariable);
%         nequal = sum(historicalData == exogenousVariable);
%         xpercentiles(i) = 100 * (nless + 0.5*nequal) / length(historicalData);
%     end 
%     
%    % get percentile status for each output in z from original data
%     zpercentiles = [];
%     for i = 1: size(zreshuffled,2)
%         historicalData = abs(zreshuffled(:,i)); % NOTE changed to positive values 
%         exogenousVariable = abs(z(i));
%         % Compute centile
%         nless = sum(historicalData < exogenousVariable);
%         nequal = sum(historicalData == exogenousVariable);
%         zpercentiles(i) = 100 * (nless + 0.5*nequal) / length(historicalData);
%     end 
%     %[~,c] = find(percentiles >=99);
% 
% end 


