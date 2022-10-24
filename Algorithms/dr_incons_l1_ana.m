function [data_rec, dsdr_rec, obj_func] = dr_incons_l1_ana(data_quantized, param, paramsolver, data)
% DR_INCONS_L1_ANA is the implementation of the Douglas-Rachford algorithm approximating the solution 
% of the analysis version of inconsistent l1 minimization-based audio dequantization.
%
% Pavel Záviška, Brno University of Technology, 2020


% predefined_lambdas are empirical values of lambdas for different level of quantization
predefined_lambdas = [0.0003, 0.00019, 0.000093, 0.0000057, 0.0000023, 0.000001, 0.0000004];
paramsolver.lambda = predefined_lambdas(param.wordlength - 1);  % regularization parameter (trade-off between sparsity and data fidelity)

% predefined_gammas are empirical values of gammas for different level of quantization            
predefined_gammas = [15.6; 13.7; 13.1; 16.2; 14.3; 13.4; 13.6];
paramsolver.gamma = predefined_gammas(param.wordlength - 1);  % DR parameter; here the threshold for soft thresholding

% starting point
u = data_quantized;

% definition of soft thresholding
soft = @(z, T) sign(z).*max(abs(z)-T, 0);

% dsdr process
dsdr_rec = NaN(paramsolver.maxit, 1);

% objective function process
obj_func = NaN(paramsolver.maxit, 1);

% iteration counter
cnt = 0;

while cnt < paramsolver.maxit
    cnt = cnt + 1;
    
    proj = proj_time(u, data_quantized, param.delta, param.algorithm);
    x = (1/(paramsolver.gamma + 1)) * (paramsolver.gamma * proj + u);
    
    syn = postpad(frsyn(param.F, soft(frana(param.F, 2*x-u), paramsolver.gamma*paramsolver.lambda)), param.Ls);
    u = u + syn - x;
    
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
    
    clear data_rec_tmp
    
end


data_rec = proj_time(x, data_quantized, param.delta, param.algorithm); % final projection into the constraints

end
