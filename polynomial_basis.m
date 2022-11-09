% IB 7/26/22

function [poly_basis] = polynomial_basis(N, maxdegree)


    poly_basis = fliplr(vander([1:N]));
    poly_basis = poly_basis(:,1:(maxdegree+1)); % degrees are represented in vandermonde as ith column having i-1 degree 

end 


