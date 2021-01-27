%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                                %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%    DEQUANTIZATION MAIN FILE    %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                                %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Pavel Záviška, Brno University of Technology, 2020
%
% using toolbox LTFAT
ltfatstart

addpath('Algorithms')
addpath('Sounds')
addpath('Tools')

plot_results = true; % Plot the results after the dequantization process

%% input file settings
% audio file path
fprintf('Loading audio.\n')
audio_path = 'Sounds'; 

audio_file = 'a08_violin';  %  'a08_violin'
                            %  'a16_clarinet'
                            %  'a18_bassoon'
                            %  'a25_harp'
                            %  'a35_glockenspiel'
                            %  'a41_celesta'
                            %  'a42_accordion'
                            %  'a58_guitar_sarasate'
                            %  'a60_piano_schubert'
                            %  'a66_wind_ensemble_stravinsky'

% load audio-file
[data, fs] = audioread([audio_path '\' audio_file '.wav']);

% peak-normalization
fprintf('Normalizing audio. ')
maxAbsVal = max(abs(data));
data = data/maxAbsVal;
fprintf('(The original maximum was %f)\n', maxAbsVal)

% signal length
param.Ls = length(data);


%% General settings
fprintf('Setting up the frame parameters.\n')
param.wordlength = 8;     % set the wordlength in bits

param.algorithm = 'DR_cons_l1_syn'; % algorithm to compute declipping, options: 
                                         % 'DR_cons_l1_syn', 
                                         % 'CP_cons_l1_ana', 
                                         % 'A_SPADQ', 
                                         % 'S_SPADQ', 
                                         % 'S_SPADQ_DR', 
                                         % 'FISTA_incons_l1_syn'
                                         % 'DR_incons_l1_syn'
                                         % 'CP_incons_l1_ana'
                                         % 'DR_incons_l1_ana'
                                         % 'FISTA_incons_l1_ana'

fprintf('Setting up algorithm: %s.\n', param.algorithm)


%% Settings for l1-minimization algorithms (CP, DR, FISTA)
if any(strcmp(param.algorithm, {'DR_cons_l1_syn', 'CP_cons_l1_ana', 'FISTA_incons_l1_syn', ...
        'DR_incons_l1_syn', 'CP_incons_l1_ana', 'DR_incons_l1_ana', 'FISTA_incons_l1_ana'}))
    
    % frame settings
    param.w = 8192;       % window length           
    param.a = param.w/4;  % window shift
    param.M = 2*param.w;  % number of frequency channels
    param.wtype = 'hann'; % window type, options available on: http://ltfat.github.io/doc/sigproc/firwin.html
    
    % construction of the frame
    fprintf('Creating the frame.\n')
    param.F = frametight(frame('dgtreal', {param.wtype, param.w}, param.a, param.M));
    param.F = frameaccel(param.F, param.Ls);  % precomputation for a fixed signal length
    
    % general settings of the l1 minimization algorithms (algorithm parameters are set directly in the respective m-file.)
    paramsolver.maxit = 500;    % maximum number of iterations
    paramsolver.minit = 25 ;    % minimum number of iterations 
    paramsolver.verbose = 0;    % display parameter
    paramsolver.comp_dsdr = 1;  % compute and store dSDR during iterations
    paramsolver.dsdr_decterm = 0;  % terminate algorithm if the SDR value starts to decrease
    paramsolver.comp_obj = 1;   % compute and store objective function values during iterations
    
end


%% Settings for SPADQ algorithms
if any(strcmp(param.algorithm, {'A_SPADQ', 'S_SPADQ', 'S_SPADQ_DR'}))

    % window parameters
    param.w = 8192;       % window length
    param.a = param.w/4;  % window shift
    param.wtype = 'hann'; % options available on: http://ltfat.github.io/doc/sigproc/firwin.html
    
    % DFT parameters 
    fprintf('Creating the frame.\n')
    param.F = frame('dft');
    param.F.redundancy = 2;  % redundancy of the DFT transform
    param.F.frana = @(insig)dft([insig; zeros(length(insig)*(param.F.redundancy-1),1)]); % redundant analysis
    param.F.frsyn = @(insig)postpad(idft(insig),length(insig)/param.F.redundancy);       % redundant synthesis   
   
    % general settings of the SPADQ algorithms
    paramsolver.verbose = 0;
    paramsolver.comp_sdr = 1;
    paramsolver.comp_obj = 1;
end


%% quantization
fprintf('Generating quantized signal.\n')
[data_quant, param.delta] = quant(data, param.wordlength); % quantizing the original signal and computing the quantization step


%% Optimization algorithm
fprintf('Starting the optimization algorithm.\n')

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


%% Evaluation
fprintf('Dequantization finished.\n')

% time
fprintf('Result obtained in %4.3f seconds.\n', time);

% SDR
sdr_quant = sdr(data, data_quant);
sdr_rec = sdr(data, data_rec);
dsdr = sdr_rec - sdr_quant;
fprintf('SDR of the quantized signal is %4.3f dB.\n', sdr_quant);
fprintf('SDR of the reconstructed signal is %4.3f dB.\n', sdr_rec);
fprintf('SDR improvement is %4.3f dB.\n', dsdr);


%% Plot results

if plot_results
    
    % Plot of signal waveforms
    t = linspace(0, length(data)/fs, length(data));
    figure
    subplot(3,1,1);
    plot(t, data);
    xlabel('time (s)')
    ylabel('Amplitude')
    title('Original')
    subplot(3,1,2);
    plot(t, data_quant);
    xlabel('time (s)')
    ylabel('Amplitude')
    title('Quantized')
    subplot(3,1,3);
    plot(t, data_rec);
    xlabel('time (s)')
    ylabel('Amplitude')
    title('Restored')
    
    % Plot signal spectrograms
    figure
    subplot(1,3,1)
    sgram(data, 'fs', fs, 'dynrange', 80);
    title('Original')
    subplot(1,3,2)
    sgram(data_quant, 'fs', fs, 'dynrange', 80);
    title('Quantized')
    subplot(1,3,3)
    sgram(data_rec, 'fs', fs, 'dynrange', 80);
    title('Restored')
    
    % Plot of dSDR and objective function
    if paramsolver.comp_dsdr && paramsolver.comp_obj
        figure
        len = length(dsdr_iter(~isnan(dsdr_iter)));
        t = linspace(0, time, len);
        yyaxis left
        p1 = plot(t.',dsdr_iter(~isnan(dsdr_iter)));
        hold on
        ylabel('{\Delta}SDR (dB)');
        
        yyaxis right
        p2 = plot(t.', obj_iter(~isnan(obj_iter)));
        title('\Delta{}SDR and objective function over time');
        xlabel('time (s)');
        ylabel('Objective function');
        grid on
        
    elseif paramsolver.comp_dsdr
        figure
        len = length(dsdr_iter(~isnan(dsdr_iter)));
        t = linspace(0, time, len);
        p1 = plot(t.',dsdr_iter(~isnan(dsdr_iter)));
        title('\Delta{}SDR over time')
        xlabel('time (s)');
        ylabel('{\Delta}SDR (dB)');
        
    elseif paramsolver.comp_obj
        figure
        len = length(dsdr_iter(~isnan(dsdr_iter)));
        t = linspace(0, time, len);
        p2 = plot(t.', obj_iter(~isnan(obj_iter)));
        title('Objective function over time')
        xlabel('time (s)');
        ylabel('Objective function');
        
    end
end
