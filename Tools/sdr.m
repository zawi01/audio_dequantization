function  [sdr_val] = sdr (u, v)
% function that computes signal-to-distortion ratio (SDR), where u is 
% considered the original signal and v is the processed signal.
%
% Pavel Z�vi�ka, Brno University of Technology, 2020

sdr_val = 20*log10(norm(u)/norm(u-v));

end
