% DESCRIPTION: Function to reconstruct timedomain signal based on BDPN outputs
% Custom idct plus polynomial representation

function [reconsig] = BPDN_reconsig(f_k, x, scalefactor, maxdegree, z, t_full)

    recon_osc = idct_custom(f_k, x, scalefactor); % custom function to do idct based on f_k, matlab idct would give different output. 
    polycomp = [];
    for i = 1:(maxdegree+1)
        polycomp(i,:) = t_full.^(i-1);
    end 
    reconsig = recon_osc + sum(z.*polycomp);

end
