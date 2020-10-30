%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                                %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%        DEQUANTIZATION          %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%        WHOLE DATABASE          %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                                %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Pavel Záviška, Brno University of Technology, 2020

% using toolbox LTFAT
ltfatstart

addpath('Algorithms')
addpath('Sounds')
addpath('Tools')

algorithms = {'DR_cons_l1_syn', 'CP_cons_l1_ana', 'A_SPADQ', 'S_SPADQ', 'S_SPADQ_DR', 'FISTA_incons_l1_syn', ...
        'DR_incons_l1_syn', 'CP_incons_l1_ana', 'DR_incons_l1_ana', 'FISTA_incons_l1_ana'};

sounds = {'a08_violin', 'a16_clarinet', 'a18_bassoon', 'a25_harp', 'a35_glockenspiel', 'a41_celesta', ...
          'a42_accordion', 'a58_guitar_sarasate', 'a60_piano_schubert', 'a66_wind_ensemble_stravinsky'};

% here set algorithms, sounds, and wordlenghts to compute.      
alg_idxs = 1:length(algorithms);
sound_idxs = 1:length(sounds);      
wordlengths = 2:8;    
    
STORE_dSDR_PROCESS = true;
STORE_OBJ_PROCESS = true;
STORE_DEQ_SOUNDS = true;

% initialization of matrices
SDR = NaN(length(sounds), length(wordlengths));
dSDR = NaN(length(sounds), length(wordlengths));
time_mat = NaN(length(sounds), length(wordlengths));

if STORE_dSDR_PROCESS == true
    dSDR_process = cell(length(sounds), length(wordlengths));
end

if STORE_OBJ_PROCESS == true
    objective_process = cell(length(sounds), length(wordlengths));
end

% initialization of counter
cnt = 0;
cases = length(alg_idxs)*length(sound_idxs)*length(wordlengths);

for algorithm = alg_idxs
    for sound = sound_idxs
        for wordlength = wordlengths
            
            cnt = cnt+1;
            %% input file settings
            eval(['[data, fs] = audioread(''Sounds\', sounds{sound}, '.wav'');']);

            % peak-normalization
            maxAbsVal = max(abs(data));
            data = data/maxAbsVal;

            % signal length
            param.Ls = length(data);


            %% General settings
            param.wordlength = wordlength;     % set the wordlength in bits
            param.algorithm = algorithms{algorithm}; % algorithm to compute declipping, options: 'DR', 'CP', 'A-SPADQ', 'S-SPADQ', 'S-SPADQ_DR', 'FISTA'


            %% Settings for l1-minimization algorithms (CP, DR)
            if any(strcmp(param.algorithm, {'DR_cons_l1_syn', 'CP_cons_l1_ana', 'FISTA_incons_l1_syn', ...
                'DR_incons_l1_syn', 'CP_incons_l1_ana', 'DR_incons_l1_ana', 'FISTA_incons_l1_ana'}))
               
                % frame settings
                param.w = 8192;
                param.a = param.w/4;
                param.M = 2*param.w; % M >= w
                param.wtype = 'hann'; % options available on: http://ltfat.github.io/doc/sigproc/firwin.html

                % construction of frame
                param.F = frametight(frame('dgtreal', {param.wtype, param.w}, param.a, param.M));
                param.F = frameaccel(param.F, param.Ls);  % precomputation for a fixed signal length

                % general settings of the l1 minimization algorithms (algorithm parameters are set directly in the respective m-file.)
                paramsolver.maxit = 500;    % maximum number of iterations
                paramsolver.minit = 25 ;    % minimum number of iterations 
                paramsolver.verbose = 0;    % display parameter
                paramsolver.comp_dsdr = STORE_dSDR_PROCESS;  % compute and store dSDR during iterations
                paramsolver.dsdr_decterm = 0;  % terminate algorithm if the SDR value starts to decrease
                paramsolver.comp_obj = STORE_OBJ_PROCESS;   % compute and store objective function values during iterations
    
            end

            %% Settings for SPADQ algorithms
            if any(strcmp(param.algorithm, {'A_SPADQ', 'S_SPADQ', 'S_SPADQ_DR'}))
            
                % window parameters
                param.w = 8192;       % window length
                param.a = param.w/4;  % window shift
                param.wtype = 'hann'; % options available on: http://ltfat.github.io/doc/sigproc/firwin.html

                % DFT parameters
                param.F = frame('dft');
                param.F.redundancy = 2;  %non-native, our own additional parameter
                param.F.frana = @(insig)dft([insig; zeros(length(insig)*(param.F.redundancy-1),1)]);
                param.F.frsyn = @(insig)postpad(idft(insig),length(insig)/param.F.redundancy);
                
                % general settings of the SPADQ algorithms
                paramsolver.verbose = 0;
                paramsolver.comp_sdr = STORE_dSDR_PROCESS;
                paramsolver.comp_obj = STORE_OBJ_PROCESS;
            end


            %% quantization
            [data_quant, param.delta] = quant(data, param.wordlength); % quantizing the original signal and computing the quantization step


            %% Optimization algorithm

            tic;

            switch param.algorithm
                case {'DR_cons_l1_syn'} % consistent l1-minimization using synthesis model of the signal, Douglas-Rachford algorithm      
                    [data_rec, dsdr_iter, obj_iter] = dr_cons_l1_syn(data_quant, param, paramsolver, data);

                case {'CP_cons_l1_ana'} % consistent l1-minimization using analysis model of the signal, Chambolle-Pock algorithm
                    [data_rec, dsdr_iter, obj_iter] = cp_cons_l1_ana(data_quant, param, paramsolver, data);

                case {'A_SPADQ', 'S_SPADQ', 'S_SPADQ_DR'} % non-convex l0-minimization based on ADMM, SPADQ algorithms
                    % paramsolver parameters
                    paramsolver.s = 1;   % increment of k
                    paramsolver.r = 1;   % every r-th iteration increment k by s
                    paramsolver.epsilon = 0.01;  % stopping criterion of termination function
                    paramsolver.maxit = ceil(floor(param.w*param.F.redundancy/2+1)*paramsolver.r/paramsolver.s); % maximum number of iterations

                    [data_rec, dsdr_iter, obj_iter] = spadq_segmentation(data_quant, param, paramsolver, data);

                case {'FISTA_incons_l1_syn'} % inconsistent l1-minimization using synthesis model of the signal, FISTA
                    [data_rec, dsdr_iter, obj_iter] = fista_incons_l1_syn(data_quant, param, paramsolver, data);

                case {'DR_incons_l1_syn'} % inconsistent l1-minimization using synthesis model of the signal, Douglas-Rachford algorithm
                    [data_rec, dsdr_iter, obj_iter] = dr_incons_l1_syn(data_quant, param, paramsolver, data);

                case {'CP_incons_l1_ana'} % inconsistent l1-minimization using analysis model of the signal, Chambolle-Pock algorithm   
                    [data_rec, dsdr_iter, obj_iter] = cp_incons_l1_ana(data_quant, param, paramsolver, data);        

                case {'DR_incons_l1_ana'} % inconsistent l1-minimization using analysis model of the signal, Douglas-Rachford algorithm
                    [data_rec, dsdr_iter, obj_iter] = dr_incons_l1_ana(data_quant, param, paramsolver, data);

                case {'FISTA_incons_l1_ana'} % inconsistent l1-minimization using analysis model of the signal, FISTA
                    [data_rec, dsdr_iter, obj_iter] = fista_incons_l1_ana(data_quant, param, paramsolver, data);

                otherwise
                    error('Invalid algorithm is set!');
            end

            time = toc;
            
            %% Time & SDR evaluation
            
            % rename and save dequantized sound
            if STORE_DEQ_SOUNDS == true, eval([sounds{sound} '_rec_' algorithms{algorithm} '_0' num2str(wordlength) ' = data_rec;']); end
            
            % store dSDR course through iterations
            if STORE_dSDR_PROCESS == true, dSDR_process{sound, wordlength-1} = dsdr_iter; end
            
            % store course of the objective function through iterations
            if STORE_OBJ_PROCESS == true, objective_process{sound, wordlength-1} = obj_iter; end
            
            % store computational time
            time_mat(sound, wordlength-1) = time;
            
            % compute and store the SDR and dSDR values of the
            % reconstructed (dequantized) signal
            sdr_quant = sdr(data, data_quant);
            sdr_rec = sdr(data, data_rec);
            SDR(sound, wordlength-1) = sdr_rec;
            dSDR(sound, wordlength-1) = sdr_rec - sdr_quant;
           
            
            disp(['Done: ', num2str(cnt), ' / ', num2str(cases)]);
            

        end
    end
    
    eval(['SDR_' algorithms{algorithm} ' = SDR;' ]);
    eval(['dSDR_' algorithms{algorithm} ' = dSDR;' ]);
    eval(['TIME_' algorithms{algorithm} ' = time_mat;' ]);
    
    if STORE_dSDR_PROCESS == true, eval(['dSDR_process_' algorithms{algorithm} ' = dSDR_process;' ]); end
    
    if STORE_OBJ_PROCESS == true, eval(['objective_process_' algorithms{algorithm} ' = objective_process;' ]); end
 
end

% clean-up
clear alg_idxs algorithm algorithms cases cnt data data_quant data_rec dSDR dsdr_iter ...
    dSDR_process fs maxAbsVal obj_iter objective_process SDR sdr_quant sdr_rec sound sound_idxs sounds ...
    STORE_DEQ_SOUNDS STORE_dSDR_PROCESS STORE_OBJ_PROCESS time time_mat wordlengths wordlength
