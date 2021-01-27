function [data_rec, dsdr_rec, obj_func] = cp_cons_l1_ana(data_quantized, param, paramsolver, data)
% CP_CONS_l1_ANA is the implementation of the Chamboll-Pock algorithm solving
% the analysis version of consistent l1 minimization-based audio dequantization.
% 
% Pavel Záviška, Brno University of Technology, 2020


% predefined_zetas are empirical values of zetas for different level of quantization
predefined_zetas = [0.0038; 0.0023; 0.0012; 0.000094; 0.000032; 0.000013; 0.0000055];
paramsolver.zeta = predefined_zetas(param.wordlength - 1); % CP parameter; here step for projection

paramsolver.sigma = 1/paramsolver.zeta;     % CP parameter; here threshold for soft thresholding
paramsolver.rho = 1;    % step size for CP algorithm rho = <0,1>

% starting poing
p = data_quantized; 
q = frana(param.F, p);

x = data_quantized; % initialization of reconstructed data vector

% definition of clip function (result of the Fenchel-Rockafellar conjugate of soft thresholding)
clip = @(x) (sign(x).*min(abs(x), 1));

% dsdr process
dsdr_rec = NaN(paramsolver.maxit, 1);

% objective function process
obj_func = NaN(paramsolver.maxit, 1);

% iteration counter
cnt = 0;

while cnt < paramsolver.maxit
    cnt = cnt+1;
    
    q = clip(q + paramsolver.sigma.*frana(param.F, x));
    
    p_old = p;
    p = proj_time(p - paramsolver.zeta*postpad(frsyn(param.F, q), param.Ls), data_quantized, param.delta, param.algorithm);
    x = p + paramsolver.rho*(p - p_old);
        
    if paramsolver.comp_obj
        obj_func(cnt) = norm(frana(param.F, x), 1); % computing objective function (l1 norm of coefficients)
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
