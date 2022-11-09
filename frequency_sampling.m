% IB 7/26/22

% INPUT
% N = number of points
% maxperiod = longest period length of interest (days)
% minperiod = shortest period length of interest (days)
% dt = difference in seconds between samples (inverse of fs)
% sampletype = 'fullrange', 'subrange' string indicating approach to
% building the frequency vector

% OUTPUT
% f_k = modified k vector for variable density frequency sampling. Typical
    % DCT uses k = 0:N-1. Here values of K are selectedt o 

% NOTE
% f_k is defined based on scaling adjustments specific to this
% implementation of DCT-II and conversion from period length in days to
% frequencies considered in Hz. For correct implementation, inputs must
% be compatible with the current framework. f_k could be modified to be a
% different function of frequency spacing. Here we have chosen linspace at
% varying sub-ranges of the default basis. 

% APPROACH 1: fullrange. minperiod:maxperiod linspaced
% -> BPDN with this frequency vector is trash. Don't use. Orthogonality of
% basis is lost.

% APPROACH 2: subrange. take default f, then linspace a subrange of it to
% improve resolution around that frequency. 
% -> Much better, though requires some tuning in simulation depending on
% frequencies of interest. 

% APPROACH 3: densesubrange.  take default f, then linspace a subrange of it to
% improve resolution around that frequency but more densely - increase
% resolution in subrange by factor of 3 at expense of random cutting some
% representation from the highest frequency end.
% -> Comparable to approach 2. Preferred method. 

function [f_k, desiredperioddays] = frequency_sampling(N,maxperiod, minperiod, dt, sampletype)
    
    fs = 1/dt;
    % APPROACH 1: Linspace entire range of interest -> 
    if strcmp(sampletype,'fullrange') == 1
        
        desiredperioddays = linspace(maxperiod, minperiod, N);
        f_k = 2*N ./ (desiredperioddays * fs *24 * 60 * 60);

    % APPROACH 2: Linspace sub-range of interest, insert into default f dct
    % Default freq vector for regular DCT
    elseif strcmp(sampletype,'subrange') == 1
        
        f_default = fs * (1:N) ./ (2*N);
        p = 1./f_default; % convert frequency to period
        p_default = p ./ (24*60*60); % convert seconds to days
        
        if minperiod < min(p_default)
            error('min period is below default period range')
        end 

        % Find where your periods of interest fit in the default vector
        i1 = max(find(p_default >= maxperiod));
        i2 = min(find(p_default <= minperiod));

        p_insert = linspace(maxperiod, minperiod, i2-i1+1);
        desiredperioddays = p_default;
        desiredperioddays(i1:i2) = p_insert;
        f_k =  2*N ./ (desiredperioddays * fs *24 * 60 * 60);
        

    % APPROACH 3:
    elseif strcmp(sampletype, 'densesubrange') == 1
        
        f_default = fs * (1:N) ./ (2*N);
        p = 1./f_default; % convert frequency to period
        p_default = p ./ (24*60*60); % convert seconds to days
        
        if minperiod < min(p_default)
            error('min period is below default period range')
        end 
        
        % Find where your periods of interest fit in the default vector
        i1 = max(find(p_default >= maxperiod));
        if length(i1) == 0
            error('maxperiod not present in default reference frequency vector')
        end 
        i2 = min(find(p_default <= minperiod));
        if length(i2) == 0
            error('minperiod not present in default reference frequency vector')
        end 
        p_insert = linspace(maxperiod, minperiod, (i2-i1+1) *3);        
        desiredperioddays = horzcat(p_default(1:i1), p_insert, p_default(i2:end)); % slice in the p_insert
        
%         % Cut out bottom end VERSION 1
%         desiredperioddays = desiredperioddays(1:length(p_default));
        
        % Cut random bits from bottom two thirds VERSION 2
        range2cut = [(length(p_default) - (floor(length(p_default) * (2/3)))): length(p_default)];
        toremove = randsample(range2cut, length(desiredperioddays) - length(p_default));
        desiredperioddays(toremove) = [];
        
  
        f_k =  2*N ./ (desiredperioddays * fs *24 * 60 * 60);
        
     
    end 
    
    fprintf(strcat("Oscillation representation may not be ideal for periods greater than", {' '}, num2str(maxperiod),{' '},"days. Check that this aligns with your application."))
    
    
end 