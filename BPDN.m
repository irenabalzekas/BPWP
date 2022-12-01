% DESCRIPTION: Basis pursuit denoising (dct and polynomial components) 
% Function to calculate basis pursuit denoising with polynomial trend

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
A = dct_basis(phi,:); 

% Polynomial basis
B = poly_basis(phi,:); 
Bplus = pinv(B); % pinv is Moore-Penrose inverse 

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
