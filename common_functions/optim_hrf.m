function [hrf] = optim_hrf(p,p2)
% Haemodynamic response function to be optimized for the time series of
% each voxel
% FORMAT [hrf] = spm_hrf(p)
% p    - parameters of the response function (two Gamma functions)
%                                                           
%        p(1) - delay of response (relative to onset)          
%        p2 - delay of undershoot (relative to onset)                          
%        p(2) - ratio of response to undershoot        

% Eval = false; Only used for finding optimum parameters for gamma fit


eval = false;
p_clip = 14;

if nargin<2
    p2 = 16;
end

if eval
    RT  = 1.4/16; %TR. May be different from the actual TR of the experiment. Depends upon the fir_order, length in the previous model.
else
    RT = 1;
end

fMRI_T = 16; % Microtime resolution.
p_start = 0;
p_end = 32;

dt  = RT/fMRI_T;
u   = [0:ceil(p_end/dt)] - p_start/dt;
hrf = spm_Gpdf(u,p(1),dt) - spm_Gpdf(u,p2,dt)/p(2);

boxwind = tukeywin(length(hrf),0.2);
hrf = hrf.*boxwind';

hrf = hrf([0:floor(p_end/RT)]*fMRI_T + 1);
hrf = hrf'/sum(hrf);

if ~eval
hrf = hrf(1:p_clip);
end

end
