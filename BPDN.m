% Basis pursuit denoising (dct and polynomial components)
% IB 7.1.22

% INPUTS
% delta = integer; % delta parameter (noise term that is aggregate
% measurement and model error)
% measureddata = column vector of behavioral sample values
% time = column vector of timestamps (in seconds) that behavioral samples
% occured 
% pmin = minimum length period of interest in days (for subsampling basis)
% pmax = maximum length period of interest in days for subsampling basis)
% maxdegree = integer, maximum degree of polynomial to account for in
% polynomial basis 

% OUTPUTS
% x = dct coefficients
% z = polynomial coefficients


function [x, z, A, B, cvx_status, cvx_optval] = BPDN(delta, measureddata, phi, dct_basis, poly_basis)


% DCT basis 
A = dct_basis(phi,:); % dct_basis = psi

% Polynomial basis
B = poly_basis(phi,:); % poly_basis = tau
Bplus = pinv(B); % pinv is moore-penrose inverse 

I = eye(size(A,1));

% CONVEX OPTIMIZATION
cvx_begin quiet
cvx_solver sdpt3
    variable x(size(A,2))
    minimize power(2,norm((I-B*Bplus)*(measureddata - A*x),2)) 
    subject to
        norm(x,1) <= delta

cvx_end
cvx_status;
cvx_optval;

% solve for z (from variable projection)
z = pinv(B)*(measureddata - A*x);

end 



%% OLD VERSION

% 
% function [x, z, dt, A, B, phi, f_adj, f_final, t_full, cvx_status, cvx_optval] = BPDN(delta, measureddata, time, pmin, pmax, maxdegree)
% 
% % ORGANIZE INPUTS
% t = time; % timestamps for each sample in measured % ADD MODD TO CHANGE SCALE BASED ON NUMBER OF DIGITS IN POSIX TIME
% % % remove any instances where spacing between samples was less than 1 hour
% % diffs = diff(t);
% % toremove = find(diffs < (1 * 60* 60));% convert one hour to seconds
% % measureddata(toremove) = [];
% % t(toremove) = [];
% y = measureddata; 
% 
% % BUILD OPTIMAL (SUBSAMPLED BASIS) * note - could probably do more elegantly
% yearsrepresented = round((max(t)-min(t))/(60*60*24) / 365,2);
% if yearsrepresented < .99
%     fprintf('BPDN: Less than 1 year represented in data, use caution in evaluating multi-month cycles \n%s', 'subsampled basis will have limited representation of low frequency components');
% end
% 
% n_init = yearsrepresented * 30000; % arbitrarily chosen 19999
% dt_init =(max(t)-min(t))/n_init; % arbitrarily chosen, was originally 3
% Fs = 1/dt_init; % sampling period and rate are inverse 
% f = Fs*(0:n_init-1)/(2*n_init); % ORIG. unclear if should be 0:n-1 or 1:n.
% f_adj_init = 1./(f*24*60*60); 
% psi_init = dct(eye(n_init,n_init)); % transform(identity matrix), this dct is orthogonal
% fmax_ind = min(find(f_adj_init< pmin));
% fmin_ind = max(find(f_adj_init > pmax));
% f_adj = f_adj_init(fmin_ind:fmax_ind);
% psi_final = psi_init(fmin_ind:fmax_ind, fmin_ind:fmax_ind); % subsample basis
% f_final = f(fmin_ind:fmax_ind);
% 
% % USE NEW N TO DEFINE NEW DT
% n = fmax_ind-fmin_ind+2; 
% dt = (max(t) - min(t))/(n-1);
% dt = (max(t) - min(t))/(n-2);
% 
% % REDEFINE SAMPLING BASED ON DTs
% % Adjust sample spacing (fill in zeros) based on temporal spacing
% T_0 = min(t); % scale sample time array (based on distance in time from start-time)
% t_hat = floor((t-T_0)/dt)+1; % indices for when each time data were sampled
% %t_full = min(t_hat): max(t_hat); % initialize full time array (temporal context, as if regularly sampled)
% t_full = min(t_hat): n; % initialize full time array (temporal context, as if regularly sampled)
% phi = t_hat;% phi = sampling matrix. points in time where y samples occured. 
% 
% % DCT basis 
% A = psi_final(phi,:);
% 
% % Polynomial basis
% tau = fliplr(vander([1:n]));
% tau = tau(:,1:(maxdegree+1)); % degrees are represented in vandermonde as ith column having i-1 degree 
% B = tau(phi,:);
% 
% % CONVEX OPTIMIZATION
% cvx_begin
%     variable x(size(A,2))
%     variable z(size(B,2))
%     minimize norm(x,1)
%     subject to
%         power(2,norm(y - A*x - B*z,2)) <= delta
% 
% cvx_end
% cvx_status;
% cvx_optval;
% 
% end 


% %% Another old version
% 
% % CONVEX OPTIMIZATION
% cvx_begin
% cvx_solver sdpt3
%     variable x(size(A,2))
%     variable z(size(B,2))
%     minimize norm(x,1)
%     subject to
%         power(2,norm(measureddata - A*x - B*z,2)) <= delta
