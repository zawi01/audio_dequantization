function [data_rec, dsdr_rec, obj_func] = fista_incons_l1_syn(data_quantized, param, paramsolver, data)
% FISTA_INCONS_l1_SYN is the implementation of the FISTA solving
% the synthesis version of inconsistent l1 minimization-based audio dequantization.
%
% Pavel Záviška, Brno University of Technology, 2020


paramsolver.mu = 1; % FISTA parameter; ensured convergency in case of tight frame for mu <= 1

% predefined_lambdas are empirical values of lambdas for different level of quantization
predefined_lambdas = [0.0003, 0.00019, 0.000093, 0.0000057, 0.0000023, 0.000001, 0.0000004];
paramsolver.lambda = predefined_lambdas(param.wordlength - 1); % regularization parameter (trade-off between sparsity and data fidelity)

% dsdr process
dsdr_rec = NaN(paramsolver.maxit, 1);

% objective function process
obj_func = NaN(paramsolver.maxit, 1);

% definition of soft thresholding
soft = @(z, T) sign(z).*max(abs(z)-T, 0);

% initialization of alpha
c = frana(param.F, data_quantized);

% initialization of dual variable u
z = c;

% initialization of variable t
t = 1;

% iteration counter
cnt = 0;

while cnt < paramsolver.maxit
    cnt = cnt + 1;
    
    syn = postpad(frsyn(param.F, z), param.Ls);
    proj = proj_time(syn, data_quantized, param.delta, param.algorithm);
    
    c_new = soft(z - paramsolver.mu * frana(param.F, syn - proj), paramsolver.mu * paramsolver.lambda);
    
    t_new = (1 + sqrt(1+4*t^2)) / 2;
    
    z = c_new + ((t - 1)/t_new)*(c_new - c);
    
    c = c_new;
    t = t_new;
    
   
    if paramsolver.comp_obj % computing objective function (lambda * l1 norm of coefficients + distance function from the feasible set)
        data_rec_tmp = postpad(frsyn(param.F, c), param.Ls); % reconstructed signal
        obj_func(cnt) = paramsolver.lambda * norm(c, 1) + norm(data_rec_tmp - proj_time(data_rec_tmp, data_quantized, param.delta, param.algorithm), 2).^2; 
    end
    
    if paramsolver.comp_dsdr || paramsolver.dsdr_decterm
        data_rec_tmp = postpad(frsyn(param.F, c), param.Ls); % reconstructed signal
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

data_rec = postpad(frsyn(param.F, c), param.Ls); % reconstructed signal
