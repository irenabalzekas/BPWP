% DESCRIPTION: Function to format samples and timestamps for method; i.e. define phi, N, etc.

function [phi, N, t_full] = BPDN_setupData(t, dt)

    N = floor((max(t) - min(t))/dt)+1;
    T_0 = min(t); % scale sample time array (based on distance in time from start-time)
    phi = floor((t-T_0)/dt) + 1; % indices for when each time data were sampled
    t_full = min(phi):N; % initialize full time array (temporal context, as if regularly sampled)

end 
