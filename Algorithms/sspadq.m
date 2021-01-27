function [data_rec, sdr_iter, obj_val] = sspadq(data_quantized, param, paramsolver, data_orig, gana)
% SSPADQ computes the synthesis version of dequantization algorithm SPADQ.
% 
% Input parameters
%       data_quantized vector of quantized signal
%       param          structure of parameters containing:
%                            Ls         length of the original signal
%                            theta      clipping threshold
%                            w          window length (in samples)
%                            a          window shift (in samples)
%                            wtype      string defining the type of the window, see http://ltfat.github.io/doc/sigproc/firwin.html
%                            F          frame (usually DFT frame)
%                            algorithm  string used to select A-SPADE or S-SPADE
%                            masks      structure containing clipping masks
%       paramsolver    structure of parameters for SPADE containing:
%                            s          relaxation stepsize
%                            r          relaxation steprate
%                            epsilon    stopping threshold of termination function
%                            maxit      maximal possible number of iterations with particular settings
%                            comp_sdr  switch to enable computing SDR in each iteration
%                            comp_obj  switch to enable storing the value of termination function in each iteration
%       data_orig      vector of original clean signal used to compute SDR during iterations
%       gana           analysis window used in the projection
%
% Output parameters
%       data_rec       vector of restored (declipped) block of signal
%       sdr_iter       SDR values during iterations the blocks of signal
%       obj_iter       Termination function values during iterations for the block of signal
%
%
%
%   z_hat (reconstructed signal coefficients)
%   z_bar (signal coefficients after hard thresholding)
%   u (dual variable, vector of signal coefficients)
%   k (required sparsity)
% 
% Pavel Záviška, Brno University of Technology, 2020


% initialization of variables
y = data_quantized;
z_hat = frana(param.F, y);
u = 0;
k = paramsolver.s;
cnt = 1;
bestObj = Inf;

% initialization of vectors for SDR and termination function during iterations
sdr_iter = NaN(paramsolver.maxit, 1);
obj_val = NaN(paramsolver.maxit, 1);


while cnt <= paramsolver.maxit
    
   % set all but k largest coefficients to zero
   % (complex conjugate pairs are taken into consideration)
   z_bar = hard_thresholding(z_hat + u, k);
      
   objVal = norm(z_hat - z_bar); % update termination function

   if objVal <= bestObj
       data_rec_coef = z_hat;
       bestObj = objVal;
   end
   
   % termination step
   if objVal <= paramsolver.epsilon
       break
   end
   
   % projection onto the set of feasible solutions
   b = z_bar - u; 
   z_hat = proj_parse_frame(b, param.F, data_quantized, param.delta .* fftshift(gana), param.algorithm);
   
   % computation and storing of SDR and termination function if requested
   if paramsolver.comp_sdr
       sdr_iter(cnt) = sdr(data_orig, frsyn(param.F, z_hat));
   end
   if paramsolver.comp_obj
       obj_val(cnt) = objVal;
   end
   
   % dual variable update
   u = u + z_hat - z_bar;
   
   % iteration counter update
   cnt = cnt+1;
   
   % incrementation of variable k (require less sparse signal in next iteration)
   if mod(cnt,paramsolver.r) == 0 
       k = k + paramsolver.s;
   end

end

% final synthesis from reconstructed coefficients and cropping to correct length
data_rec = frsyn(param.F, data_rec_coef);

end
