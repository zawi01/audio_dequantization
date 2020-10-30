function [proj] = proj_parse_frame(c, F, data_quantized, delta)
% PROJ_PARSE_FRAME perfoms projection of coefficients c
% onto a multidimensional interval for the use of audio dequantization.
% 
% Input parameters
%       c               vector of input coefficients
%       F               used frame for synthesis and analysis
%       data_quantized  original clipped signal
%
% This projection
%       proj(z) = argmin_{u} ||z - u||_2 s.t. Au \in [b_1, b_2]
% can be evaluated as 
%       proj(z) = z-A^+(Az - proj(Az));  %here A^+ denotes pseudoinverse
% 
% The projection proj(Az) is computed by function proj_time.
% 
% Please note that this particular function works only for Parseval tight frame 
% (only in this case the pseudoinverse is identical to the analysis of the signal)
%
% Pavel Z�vi�ka, Brno University of Technology, 2020


% Synthesis of the signal 
syn = postpad(frsyn(F, c), length(data_quantized));

% Compute proj(Az)
proj_temp = proj_time(syn, data_quantized, delta);

% Final projection
proj = c - frana(F, syn-proj_temp);

end
