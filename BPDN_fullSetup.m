% DESCRIPTION: Wrapper for other BPDN functions, setup, calculating bases, etc. 

function [phi, N, t_full, f_k, desiredperioddays, dct_basis, scalefactor, poly_basis] = BPDN_fullSetup(t, dt, maxperiod, minperiod, sampletype, maxdegree)

    [phi, N, t_full] = BPDN_setupData(t, dt);
    [f_k, desiredperioddays] = frequency_sampling(N, maxperiod, minperiod, dt, sampletype);
    [dct_basis, scalefactor] = DCT2_basis(N, f_k); % making the basis takes a long time. want to pre-define it before running BPDN for loops and loops 
    [poly_basis] = polynomial_basis(N, maxdegree);

end 
