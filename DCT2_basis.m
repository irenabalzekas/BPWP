% DESCRIPTION: Function to create DCT-II basis based on N and defined frequencies

% INPUT
% N = number of points
% f_k = k vector for non-uniform (variable density) frequency sampling.
% should be pre-converted. Normally, k = 0:N-1 for example....

        % desiredperioddays = linspace(maxperiod, minperiod, N);
        % f_k = 2 * N ./ (desiredperioddays * fs * 24 * 60 * 60);

% OUTPUT
% basis = N x N dct-II basis based on frequency vector. Used same expression 
        % as Matlab: https://www.mathworks.com/help/signal/ref/dct.html

% NOTE
% Mutual coherence: One way to evaluate output basis. Essentially gives PSF
% of operator. Want to be as close to a delta function as possible.
% Measures relation from peak to next side lobe. If playing with different
% frequency representations and their impact on basis, play around with
% mutual coherence. [mc] = mut_coh(X)
        

function [basis, scalefactor] = DCT2_basis(N, f_k)
    
    A = eye(N,N); % input data 
    basis = zeros(N,N);
    for r = 1: size(A,1)

        tempdat = A(r,:);

        ks = zeros(1,size(A,2));
        for k = 1:length(f_k)

            if k == 1 
                    filler = 1;
            elseif k ~= 1
                    filler = 0;
            end

            n = 1:N;
            per_k= tempdat.*(1/sqrt(1+filler)).*cos((pi/(2*N)) * (2.*(n)-1) * (f_k(k)-1)); 
            ks(k) = sqrt(2/N).*sum(per_k);  

        end 
        basis(r,:) = ks;
    end 
    % rescale basis (non-uniform sampling makes dct_manual non-orthogonal, 
            % dividing by 2-norm rescales it and matches coefficients to orthogonal,
            % uniformly sampled approach)
    scalefactor = normest(basis);
    basis = basis./scalefactor; 

end
