function [data_rec, dsdr_rec, obj_func, cnt] = cp_incons_l1_ana(data_quantized, param, paramsolver, data)
% CP_INCONS_L1_ANA is the implementation of the Chamboll-Pock algorithm solving 
% the analysis version of inconsistent l1 minimization-based audio dequantization.
%
% Pavel Záviška, Brno University of Technology, 2020


% predefined_lambdas are empirical values of lambdas for different level of quantization  
predefined_lambdas = [0.0003, 0.00019, 0.000093, 0.0000057, 0.0000023, 0.000001, 0.0000004];
paramsolver.lambda = predefined_lambdas(param.wordlength - 1); % regularization parameter (trade-off between sparsity and data fidelity)

% predefined_zetas are empirical values of zetas for different level of quantization
predefined_zetas = [10; 10; 10; 10; 10; 10; 10];
paramsolver.zeta = predefined_zetas(param.wordlength - 1);  % CP parameter; here step for projection

paramsolver.sigma = 1/paramsolver.zeta;     % CP parameter; here threshold for soft thresholding
paramsolver.rho = 1;    % step size for CP algorithm rho = <0,1>


% starting poing
p = data_quantized; 
c = frana(param.F, p);

x = data_quantized; % initialization of reconstructed data vector

% definition of clip function (result of the Fenchel-Rockafellar conjugate of soft thresholding)
clip = @(x, T) (sign(x).*min(abs(x), T));

% dsdr process
dsdr_rec = NaN(paramsolver.maxit, 1);

% objective function process
obj_func = NaN(paramsolver.maxit, 1);

% iteration counter
cnt = 0;

while cnt < paramsolver.maxit
    cnt = cnt+1;
    
    c = clip(c + paramsolver.sigma.*frana(param.F, x), paramsolver.lambda);
    u = p - paramsolver.zeta * postpad(frsyn(param.F, c), param.Ls);
    
    proj = proj_time(u, data_quantized, param.delta, param.algorithm);
    
    p_old = p;
    p = (1/(paramsolver.zeta + 1)) * (paramsolver.zeta * proj + u);
    x = p + paramsolver.rho*(p - p_old);
        
    if paramsolver.comp_obj % computing objective function (lambda * l1 norm of coefficients + distance function from the feasible set)
        obj_func(cnt) = paramsolver.lambda * norm(frana(param.F, x), 1) + norm(x - proj_time(x, data_quantized, param.delta, param.algorithm), 2).^2; 
    end
    
    if paramsolver.comp_dsdr || paramsolver.dsdr_decterm
        dsdr_rec(cnt) = sdr(data, x) - sdr(data, data_quantized); % computing dSDR
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

data_rec = proj_time(x, data_quantized, param.delta, param.algorithm); % final projection into the constraints

end
