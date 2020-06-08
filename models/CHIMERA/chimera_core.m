%% CHIMERA: Combining Mechanistic Models and Machine Learning for Personalized Chemotherapy and Surgery Sequencing in Breast Cancer
%% PREPARE ENVIRONMENT
clear all; clc; close all; format long; pause(2);
%% INIT SIMULATION
% enables dynamic visualization on network runtime
DYN_VISUAL = 1;
% number of populations in the network
N_SOM      = 2;
% number of neurons in each population
N_NEURONS  = 50;
% MAX_EPOCHS for SOM relaxation
MAX_EPOCHS_IN_LRN = 100;
MAX_EPOCHS_XMOD_LRN = 100;
% decay factors
ETA = 1.0; % activity decay
XI = 1e-3; % weights decay
% enable population wrap-up to cancel out boundary effects
WRAP_ON = 0;
% init data
sensory_data.x = [];
sensory_data.y = [];
%% SELECT DATA SOURCE (arbitrary function or dataset)
DATASET = 0; % if dataset is 1 load cancer dataset, otherwise drug concentration
if DATASET == 0
    % Data from Kuh, Hyo-Jeong, et al.
    % "Computational model of intracellular pharmacokinetics of paclitaxel."
    % Journal of Pharmacology and Experimental Therapeutics 293.3 (2000): 761-770.
    
    experiment_dataset = 2; % {1, 2, 3}
    % read from sample datasets
    switch experiment_dataset
        case 1 % learning the cellular concentration of chemotoxic drug
            % reproducing Figure 1 left from the Kuh et al. 2000 paper
            opts = delimitedTextImportOptions("NumVariables", 2);
            
            % Specify range and delimiter
            opts.DataLines = [2, Inf];
            opts.Delimiter = ",";
            
            % Specify column names and types
            opts.VariableNames = ["time", "concentration"];
            opts.VariableTypes = ["double", "double"];
            
            % Specify file level properties
            opts.ExtraColumnsRule = "ignore";
            opts.EmptyLineRule = "read";
            
            % Import the data
            filename = ['..' filesep '..' filesep 'datasets' filesep '0' filesep 'Paclitaxel_10nm_cell_conc.csv'];
            Paclitaxel10nmcellconc = readtable(filename, opts);
            % Clear temporary variables
            clear opts
            
            % populate data structure
            sensory_data.x =  Paclitaxel10nmcellconc.time;
            sensory_data.y =  Paclitaxel10nmcellconc.concentration;
            
                        % save the original dataset
            sensory_data_orig = sensory_data;
            DATASET_LEN_ORIG = length(sensory_data_orig.x);
            % change range
            sensory_data.range  = 1.0;
            % convert x axis data to [-sensory_data.range, +sensory_data.range]
            minVal = min(sensory_data.x);
            maxVal = max(sensory_data.x);
            sensory_data.x = (((sensory_data.x - minVal) * (sensory_data.range - (-sensory_data.range))) / (maxVal - minVal)) + (-sensory_data.range);
            % convert y axis data to [-sensory_data.range, +sensory_data.range]
            minVal = min(sensory_data.y);
            maxVal = max(sensory_data.y);
            sensory_data.y = (((sensory_data.y - minVal) * (sensory_data.range - (-sensory_data.range))) / (maxVal - minVal)) + (-sensory_data.range);
            % load the data and extrapolate for more density in x axis
            upsample_factor = 20;
            datax = sensory_data.x';
            idx_data = 1:length(datax);
            idx_upsampled_data = 1:1/upsample_factor:length(datax);
            datax_extrapolated = interp1(idx_data, datax, idx_upsampled_data, 'linear');
            % load the data and extrapolate for more density in y axis
            datay = sensory_data.y';
            idx_data = 1:length(datay);
            idx_upsampled_data = 1:1/upsample_factor:length(datay);
            datay_extrapolated = interp1(idx_data, datay, idx_upsampled_data, 'linear');
            
            sensory_data.x = datax_extrapolated;
            sensory_data.y = datay_extrapolated;
            
            
        case 2 % learning the extracellular concentration of chemotoxic drug
            % reproducing Figure 1 right from the Kuh et al. 2000 paper
            opts = delimitedTextImportOptions("NumVariables", 2);
            
            % Specify range and delimiter
            opts.DataLines = [2, Inf];
            opts.Delimiter = ",";
            
            % Specify column names and types
            opts.VariableNames = ["time", "concentration"];
            opts.VariableTypes = ["double", "double"];
            
            % Specify file level properties
            opts.ExtraColumnsRule = "ignore";
            opts.EmptyLineRule = "read";
            
            % Import the data
            filename = ['..' filesep '..' filesep 'datasets' filesep '0' filesep 'Paclitaxel_10nm_medium_conc.csv'];
            Paclitaxel10nmmediumconc = readtable(filename, opts);
            % Clear temporary variables
            clear opts
            
            % populate data structure
            sensory_data.x =  Paclitaxel10nmmediumconc.time;
            sensory_data.y =  Paclitaxel10nmmediumconc.concentration;
            
            % save the original dataset
            sensory_data_orig = sensory_data;
            DATASET_LEN_ORIG = length(sensory_data_orig.x);
            % change range
            sensory_data.range  = 1.0;
            % convert x axis data to [-sensory_data.range, +sensory_data.range]
            minVal = min(sensory_data.x);
            maxVal = max(sensory_data.x);
            sensory_data.x = (((sensory_data.x - minVal) * (sensory_data.range - (-sensory_data.range))) / (maxVal - minVal)) + (-sensory_data.range);
            % convert y axis data to [-sensory_data.range, +sensory_data.range]
            minVal = min(sensory_data.y);
            maxVal = max(sensory_data.y);
            sensory_data.y = (((sensory_data.y - minVal) * (sensory_data.range - (-sensory_data.range))) / (maxVal - minVal)) + (-sensory_data.range);
            % load the data and extrapolate for more density in x axis
            upsample_factor = 20;
            datax = sensory_data.x';
            idx_data = 1:length(datax);
            idx_upsampled_data = 1:1/upsample_factor:length(datax);
            datax_extrapolated = interp1(idx_data, datax, idx_upsampled_data, 'linear');
            % load the data and extrapolate for more density in y axis
            datay = sensory_data.y';
            idx_data = 1:length(datay);
            idx_upsampled_data = 1:1/upsample_factor:length(datay);
            datay_extrapolated = interp1(idx_data, datay, idx_upsampled_data, 'linear');
            
            sensory_data.x = datax_extrapolated;
            sensory_data.y = datay_extrapolated;
            
        case 3 % learning the underlying volume dynamics in the concentration equations
            filename=['..' filesep '..' filesep 'datasets' filesep '0' filesep 'Paclitaxel_conc_volume_dynamics'];
            MIN_VAL         = 0.0;
            MAX_VAL         = 1.0;
            sensory_data.range  = MAX_VAL;
            
            % setup the number of random input samples to generate
            NUM_VALS        = 450;
            % generate NUM_VALS random samples in the given interval
            sensory_data.x  = MIN_VAL + rand(NUM_VALS, 1)*(MAX_VAL - MIN_VAL);
            % volume relation given paper parameters in Table 3
            % from the Kuh et al. 2000 paper
            Vonecell = 0.00002; % average volume of one cell
            ICN = 50; % initial cells number
            Kcellnum = 6; % rate constant for changes in cell number due to drug
            
            % Kuh et al. 2000 propose Vc = Vonecell * ICN * exp(Kcellnum*t)
            % we generate the data as suggested and fed the data to CHIMERA
            % to learn the underlying volume dynamics
            sensory_data.y  = Vonecell * ICN * exp(sensory_data.x * Kcellnum);
            sensory_data_orig = sensory_data;
            DATASET_LEN_ORIG = length(sensory_data_orig.x);
            
    end
    
else
    % select the dataset of interest
    experiment_dataset = 3; % {1, 2(a,b), 3}
    % read from sample datasets
    switch experiment_dataset
        case 1
            
            % Rodallec, Anne, Giacometti, Sarah, Ciccolini, Joseph, & Fanciullino, Raphaelle. (2019).
            % Tumor growth kinetics of human MDA-MB-231 cells transfected with dTomato lentivirus [Data set].
            % Zenodo. http://doi.org/10.5281/zenodo.3593919
            
            filename = ['..' filesep '..' filesep 'datasets' filesep '1' filesep 'MDA-MB-231dTomato.csv'];
            delimiter = ',';
            startRow = 2;
            formatSpec = '%f%f%f%[^\n\r]';
            % Open the text file.
            fileID = fopen(filename,'r');
            % Read columns of data according to the format.
            dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
            % Close the text file.
            fclose(fileID);
            
            % Create output variable as table
            MDAMB231dTomato = table(dataArray{1:end-1}, 'VariableNames', {'ID','Time','Observation'});
            % or as a simple matrix
            % MDAMB231dTomato = [dataArray{1:end-1}];
            
            % Clear temporary variables
            clearvars delimiter startRow formatSpec fileID dataArray ans;
            % check which ID one needs
            ID = 0; % ID is one of {0, 1, 2, 3, 4, 5, 6, 7}
            sensory_data.x =  MDAMB231dTomato.Time(MDAMB231dTomato.ID == ID);
            sensory_data.y =  MDAMB231dTomato.Observation(MDAMB231dTomato.ID == ID);
            
        case 2
            
            % Gaddy, Thomas D.; Wu, Qianhui; Arnheim, Alyssa D.; D. Finley, Stacey (2017)
            % Mechanistic modeling quantifies the influence of tumor growth kinetics on the response to anti-angiogenic treatment.
            % PLOS Computational Biology. Dataset. https://doi.org/10.1371/journal.pcbi.1005874
            
            % Import the data
            filename = ['..' filesep '..' filesep 'datasets' filesep '2' filesep 'angio-genesis.csv'];
            delimiter = ',';
            startRow = 2;
            
            formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
            
            fileID = fopen(filename,'r');
            
            dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
            fclose(fileID);
            
            S1Table = table(dataArray{1:end-1}, 'VariableNames', {'RolandTimedays','RolandVolumemm3','ZibaraTimedays','ZibaraVolumemm3','Volk2008Timedays','Volk2008Volumemm3','TanTimedays','TanVolumemm3','Volk2011aTimedays','Volk2011aVolumemm3','Volk2011bTimedays','Volk2011bVolumemm3'});
            
            clearvars delimiter startRow formatSpec fileID dataArray ans;
            
            % Add filtering for sub-dataset
            study_id = 'Volk11a'; % {Tan, Volk11a}
            switch study_id
                case 'Tan'
                    sensory_data.x = S1Table.TanTimedays(~isnan(S1Table.TanTimedays));
                    sensory_data.y = S1Table.TanVolumemm3(~isnan(S1Table.TanVolumemm3));
                case 'Volk11a'
                    sensory_data.x = S1Table.Volk2011aTimedays(~isnan(S1Table.Volk2011aTimedays));
                    sensory_data.y = S1Table.Volk2011aVolumemm3(~isnan(S1Table.Volk2011aVolumemm3));
            end
            % udpate the filename to contain also sub-dataset and study id
            filename = ['..' filesep '..' filesep 'datasets' filesep '2' filesep 'angio-genesis.csv' study_id];
            
        case 3
            
            % Mastri, Michalis, Tracz, Amanda, & Ebos, John ML. (2019).
            % Tumor growth kinetics of human LM2-4LUC+ triple negative breast carcinoma cells [Data set].
            % Zenodo. http://doi.org/10.5281/zenodo.3574531
            
            filename = ['..' filesep '..' filesep 'datasets' filesep '3' filesep 'LM2-4LUC.csv'];
            delimiter = ',';
            startRow = 2;
            formatSpec = '%f%f%f%[^\n\r]';
            % Open the text file.
            fileID = fopen(filename,'r');
            dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
            % Close the text file.
            fclose(fileID);
            
            % Create output variable for table import
            LM24LUC = table(dataArray{1:end-1}, 'VariableNames', {'ID','Time','Observation'});
            % for numeric array import
            % LM24LUC = [dataArray{1:end-1}];
            
            % Clear temporary variables
            clearvars delimiter startRow formatSpec fileID dataArray ans;
            
            % check which ID one needs
            ID = 60; % ID is one of {0, 1, 2, 3, 4, 5, ..., 65}
            sensory_data.x =  LM24LUC.Time(LM24LUC.ID == ID);
            sensory_data.y =  LM24LUC.Observation(LM24LUC.ID == ID);
    end
    % save the original dataset
    sensory_data_orig = sensory_data;
    DATASET_LEN_ORIG = length(sensory_data_orig.x);
    % change range
    sensory_data.range  = 1.0;
    % convert x axis data to [-sensory_data.range, +sensory_data.range]
    minVal = min(sensory_data.x);
    maxVal = max(sensory_data.x);
    sensory_data.x = (((sensory_data.x - minVal) * (sensory_data.range - (-sensory_data.range))) / (maxVal - minVal)) + (-sensory_data.range);
    % convert y axis data to [-sensory_data.range, +sensory_data.range]
    minVal = min(sensory_data.y);
    maxVal = max(sensory_data.y);
    sensory_data.y = (((sensory_data.y - minVal) * (sensory_data.range - (-sensory_data.range))) / (maxVal - minVal)) + (-sensory_data.range);
    % load the data and extrapolate for more density in x axis
    upsample_factor = 10;
    datax = sensory_data.x';
    idx_data = 1:length(datax);
    idx_upsampled_data = 1:1/upsample_factor:length(datax);
    datax_extrapolated = interp1(idx_data, datax, idx_upsampled_data, 'linear');
    % load the data and extrapolate for more density in y axis
    datay = sensory_data.y';
    idx_data = 1:length(datay);
    idx_upsampled_data = 1:1/upsample_factor:length(datay);
    datay_extrapolated = interp1(idx_data, datay, idx_upsampled_data, 'linear');
end
% re-assign data
if DATASET == 1
    sensory_data.x = datax_extrapolated;
    sensory_data.y = datay_extrapolated;
end
DATASET_LEN     = length(sensory_data.x);
%% CREATE NETWORK AND INITIALIZE PARAMS
% create a network of SOMs given the simulation constants
populations = create_init_network(N_SOM, N_NEURONS);
% init activity vector
act_cur = zeros(N_NEURONS, 1);
% init neighborhood function
hwi = zeros(N_NEURONS, 1);
% learning params
learning_params.t0 = 1;
learning_params.tf_learn_in = MAX_EPOCHS_IN_LRN;
learning_params.tf_learn_cross = MAX_EPOCHS_XMOD_LRN;
% init width of neighborhood kernel
sigma0 = N_NEURONS/2;
sigmaf = 1;
learning_params.sigmat = parametrize_learning_law(sigma0, sigmaf, learning_params.t0, learning_params.tf_learn_in, 'invtime');
% init learning rate
alpha0 = 0.01;
alphaf = 0.001;
learning_params.alphat = parametrize_learning_law(alpha0, alphaf, learning_params.t0, learning_params.tf_learn_in, 'invtime');
% cross-modal learning rule type
cross_learning = 'covariance';    % {hebb - Hebbian, covariance - Covariance, oja - Oja's Local PCA}
% mean activities for covariance learning
avg_act=zeros(N_NEURONS, N_SOM);
%% NETWORK SIMULATION LOOP
fprintf('Started training sequence ...\n');
% present each entry in the dataset for MAX_EPOCHS epochs to train the net
for t = 1:learning_params.tf_learn_cross
    % update visualization of the Hebbian links
    if DYN_VISUAL==1
        if DATASET == 1
            visualize_runtime(populations, t);
        end
    end    % learn the sensory space data distribution
    if(t<learning_params.tf_learn_in)
        for didx = 1:DATASET_LEN
            % loop through populations
            for pidx = 1:N_SOM
                % reinit population activity
                act_cur = zeros(populations(pidx).lsize, 1);
                % pick a new sample from the dataset and feed it to the current layer
                switch pidx
                    case 1
                        input_sample = sensory_data.x(didx);
                    case 2
                        input_sample = sensory_data.y(didx);
                end
                % compute new activity given the current input sample
                for idx = 1:populations(pidx).lsize
                    act_cur(idx) = (1/(sqrt(2*pi)*populations(pidx).s(idx)))*...
                        exp(-(input_sample - populations(pidx).Winput(idx))^2/(2*populations(pidx).s(idx)^2));
                end
                % normalize the activity vector of the population
                act_cur = act_cur./sum(act_cur);
                % update the activity for the next iteration
                populations(pidx).a = (1-ETA)*populations(pidx).a + ETA*act_cur;
                % competition step: find the winner in the population given the input data
                % the winner is the neuron with the highest activity elicited
                % by the input sample
                [win_act, win_pos] = max(populations(pidx).a);
                for idx = 1:populations(pidx).lsize % go through neurons in the population
                    % cooperation step: compute the neighborhood kernel
                    if(WRAP_ON==1)
                        % wrap up the population to avoid boundary effects
                        % dist = min{|i-j|, N - |i-j|}
                        hwi(idx) = exp(-norm(min([norm(idx-win_pos), N_NEURONS - norm(idx-win_pos)]))^2/(2*learning_params.sigmat(t)^2));
                    else
                        % simple Gaussian kernel with no boundary compensation
                        hwi(idx) = exp(-norm(idx-win_pos)^2/(2*learning_params.sigmat(t)^2));
                    end
                    % learning step: compute the weight update
                    populations(pidx).Winput(idx) = populations(pidx).Winput(idx) + ...
                        learning_params.alphat(t)*...
                        hwi(idx)*...
                        (input_sample - populations(pidx).Winput(idx));
                    % update the shape of the tuning curve for current neuron
                    populations(pidx).s(idx) = populations(pidx).s(idx) + ...
                        learning_params.alphat(t)*...
                        (1/(sqrt(2*pi)*learning_params.sigmat(t)))*...
                        hwi(idx)*...
                        ((input_sample - populations(pidx).Winput(idx))^2 - populations(pidx).s(idx)^2);
                end
            end % end for population pidx
        end % end samples in the dataset
    end % allow the som to learn the sensory space data distribution
    % % learn the cross-modal correlation
    for didx = 1:DATASET_LEN
        % use the learned weights and compute activation
        % loop through populations
        for pidx = 1:N_SOM
            % pick a new sample from the dataset and feed it to the current layer
            switch pidx
                case 1
                    input_sample = sensory_data.x(didx);
                case 2
                    input_sample = sensory_data.y(didx);
            end
            % compute new activity given the current input sample
            for idx = 1:populations(pidx).lsize
                act_cur(idx) = (1/(sqrt(2*pi)*populations(pidx).s(idx)))*...
                    exp(-(input_sample - populations(pidx).Winput(idx))^2/(2*populations(pidx).s(idx)^2));
            end
            % normalize the activity vector of the population
            act_cur = act_cur./sum(act_cur);
            % update the activity for the next iteration
            populations(pidx).a = (1-ETA)*populations(pidx).a + ETA*act_cur;
        end
        % check which learning rule we employ
        switch(cross_learning)
            case 'hebb'
                % cross-modal Hebbian learning rule
                populations(1).Wcross = (1-XI)*populations(1).Wcross + XI*populations(1).a*populations(2).a';
                populations(2).Wcross = (1-XI)*populations(2).Wcross + XI*populations(2).a*populations(1).a';
            case 'covariance'
                % compute the mean value computation decay
                OMEGA = 0.002 + 0.998/(t+2);
                % compute the average activity for Hebbian covariance rule
                for pidx = 1:N_SOM
                    avg_act(:, pidx) = (1-OMEGA)*avg_act(:, pidx) + OMEGA*populations(pidx).a;
                end
                % cross-modal Hebbian covariance learning rule: update the synaptic weights
                populations(1).Wcross = (1-XI)*populations(1).Wcross + XI*(populations(1).a - avg_act(:, 1))*(populations(2).a - avg_act(:, 2))';
                populations(2).Wcross = (1-XI)*populations(2).Wcross + XI*(populations(2).a - avg_act(:, 2))*(populations(1).a - avg_act(:, 1))';
            case 'oja'
                % Oja's local PCA learning rule
                populations(1).Wcross = ((1-XI)*populations(1).Wcross + XI*populations(1).a*populations(2).a')/...
                    sqrt(sum(sum((1-XI)*populations(1).Wcross + XI*populations(1).a*populations(2).a')));
                populations(2).Wcross = ((1-XI)*populations(2).Wcross + XI*populations(2).a*populations(1).a')/...
                    sqrt(sum(sum((1-XI)*populations(2).Wcross + XI*populations(2).a*populations(1).a')));
        end
    end % end for values in dataset
end % end for training epochs
fprintf('Ended training sequence. Presenting results ...\n');
%% VISUALIZATION
present_tuning_curves(populations(1), sensory_data);
present_tuning_curves(populations(2), sensory_data);
% normalize weights between [0,1] for display
populations(1).Wcross = populations(1).Wcross ./ max(populations(1).Wcross(:));
populations(2).Wcross = populations(2).Wcross ./ max(populations(2).Wcross(:));
if DATASET == 0 && experiment_dataset == 2
    populations(1).Wcross = rot90(rot90(populations(1).Wcross));
    populations(2).Wcross = rot90(rot90(populations(2).Wcross));
end
% visualize post-simulation weight matrices encoding learned relation
[sensory_data, neural_model] = visualize_results(sensory_data, populations, learning_params, DATASET, experiment_dataset);
% denormalize the neural model fit to match original input data range
minVal = min(neural_model);
maxVal = max(neural_model);
neural_model = (((neural_model - minVal) * (max(sensory_data_orig.y) - min(sensory_data_orig.y))) / (maxVal - minVal)) + min(sensory_data_orig.y);
% downsample to match input data size, on various datasets
neural_model = interp1(1:length(neural_model), neural_model, linspace(1,length(neural_model),DATASET_LEN_ORIG));
% save runtime data in a file for later analysis and evaluation against
% other models - imported in evaluation script
runtime_data_file = sprintf('Experiment_dataset_%s_runtime.mat',...
    filename(18:end)); % get only the name of the file remove path
save(runtime_data_file);

