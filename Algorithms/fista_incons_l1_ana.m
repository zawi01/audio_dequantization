function [data_rec, dsdr_rec, obj_func] = fista_incons_l1_ana(data_quantized, param, paramsolver, data)
% FISTA_INCONS_L1_ANA is the implementation of the FISTA approximating the solution 
% of the analysis version of inconsistent l1 minimization-based audio dequantization.
%
% Pavel Záviška, Brno University of Technology, 2020


% predefined_lambdas are empirical values of lambdas for different level of quantization  
predefined_lambdas = [0.0003, 0.00019, 0.000093, 0.0000057, 0.0000023, 0.000001, 0.0000004];
paramsolver.lambda = predefined_lambdas(param.wordlength - 1);  % regularization parameter (trade-off between sparsity and data fidelity)

paramsolver.mu = 1;  % ensured convergency in case of tight frame for mu <= 1

% dsdr process
dsdr_rec = NaN(paramsolver.maxit, 1);

% objective function process
obj_func = NaN(paramsolver.maxit, 1);

% definition of soft thresholding
soft = @(z, T) sign(z).*max(abs(z)-T, 0);

% initialization of alpha
x = data_quantized;

% initialization of dual variable u
u = x;

% initialization of variable t
t = 1;

% iteration counter
cnt = 0;

while cnt < paramsolver.maxit
    cnt = cnt + 1;
    
    proj = proj_time(u, data_quantized, param.delta, param.algorithm);
      
    x_new = frsyn(param.F, soft(frana(param.F, u - paramsolver.mu * (u - proj)), paramsolver.mu*paramsolver.lambda));
    x_new = postpad(x_new, param.Ls);
    
    t_new = (1 + sqrt(1+4*t^2)) / 2;
    
    u = x_new + ((t - 1)/t_new)*(x_new - x);
    
    x = x_new;
    t = t_new;
    
   
    if paramsolver.comp_obj % computing objective function (lambda * l1 norm of coefficients + distance function from the feasible set)
        obj_func(cnt) = paramsolver.lambda * norm(frana(param.F, x), 1) + norm(x - proj_time(x, data_quantized, param.delta, param.algorithm), 2).^2; 
    end
    
    if paramsolver.comp_dsdr || paramsolver.dsdr_decterm
        data_rec_tmp = x; % reconstructed signal
        dsdr_rec(cnt) = sdr(data, data_rec_tmp) - sdr(data, data_quantized); % computing dSDR
        if paramsolver.dsdr_decterm && cnt > paramsolver.minit && dsdr_rec(cnt) - dsdr_rec(cnt-1) < 0
            break
        end
    end
    
    if paramsolver.verbose
        fprintf(' Iteration number: %u\n', cnt);
        if paramsolver.comp_dsdr
            fprintf(' SDR improvement: %.3f dB\n', dsdr_rec(cnt));
        end
        if paramsolver.comp_obj
            fprintf(' Objective function value: %e \n', obj_func(cnt));
        end
        fprintf('\n')        
    end    
    
    
end

data_rec = proj_time(x, data_quantized, param.delta, param.algorithm); % reconstructed signal
