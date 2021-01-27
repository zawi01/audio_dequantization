function [data_rec, dsdr_rec, obj_func] = dr_cons_l1_syn(data_quantized, param, paramsolver, data)
% DR_CONS_L1_SYN is the implementation of the Douglas-Rachford algorithm solving
% the synthesis version of consistent l1 minimization-based audio dequantization.
% 
% Pavel Záviška, Brno University of Technology, 2020


% predefined_gammas are empirical values of gammas for different level of quantization 
predefined_gammas = [0.0047; 0.0026; 0.0012; 0.000095; 0.000033; 0.000013; 0.0000055];  
paramsolver.gamma = predefined_gammas(param.wordlength - 1);  % DR parameter; here the threshold for soft thresholding

% starting point
z = frana(param.F, data_quantized);

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
    
    c = proj_parse_frame(z, param.F, data_quantized, param.delta, param.algorithm);
    z = z + (soft(2*c-z, paramsolver.gamma)-c);
    
    if paramsolver.comp_obj
        obj_func(cnt) = norm(c, 1); % computing objective function (l1 norm of coefficients)
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
    
    clear data_rec_tmp
    
end


c = proj_parse_frame(z, param.F, data_quantized, param.delta, param.algorithm); % final projection into the constraints
data_rec = postpad(frsyn(param.F, c), param.Ls); % reconstructed signal

end
