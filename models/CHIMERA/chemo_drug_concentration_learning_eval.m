% CHIMERA: Combining Mechanistic Models and Machine Learning for Personalized Chemotherapy and Surgery Sequencing in Breast Cancer
% prepare environment
clc; close all;
intracell = 1; % intracell or extracellular experiment {1, 0}
% dataset name
if intracell == 1
    dataset = 'Experiment_dataset_Paclitaxel_10nm_cell_conc.csv_runtime.mat';
else
    dataset = 'Experiment_dataset_Paclitaxel_10nm_medium_conc.csv_runtime.mat';
end
load(dataset);
global T M
% data points and time indices
T = 1:DATASET_LEN_ORIG; T = T';
M = sensory_data_orig.y;
figure();
set(gcf, 'color', 'w');hold on; box off;
% dataset points
semilogy(T,M, 'ro', 'MarkerEdgeColor','r',...
    'MarkerFaceColor','r',...
    'MarkerSize',10);
% neural model
if DATASET == 0 && experiment_dataset == 2
    neural_model = fliplr(neural_model);
end
semilogy(neural_model, 'LineWidth',3);
set(gca, 'YScale', 'log');
set(gca, 'XTickLabel', {'0','6','','12','','18','','24'});
ylabel('Drug concentration \muM'); xlabel('Time (h)');

if intracell == 1
    title('Chemotoxic drug pharmacokinetics P(t,V) intracellular c(t)');
    legend('Data points', 'Learnt P(t,V) intracellular c(t)');
else
    title('Chemotoxic drug pharmacokinetics P(t,V) extracellular c_m(t)');
    legend('Data points', 'Learnt P(t,V) extracellular c_m(t)');
end
legend boxoff


% Evaluate SSE, RMSE, MAPE
% params as from [Benzekry et al., 2014c]
alfa = 0.84;
sigma = 0.21;

% SSE
SSEn = model_sse(alfa, sigma, M, neural_model);
% RMSE
RMSEn = model_rmse(alfa, sigma, 0, M, neural_model);
% sMAPE
sMAPEn = mean(2*abs((M-neural_model'))./(abs(M) + abs(neural_model')));

annotation('textbox', [0.6, 0.6, 0.3, 0.3], ...
    'String', sprintf('SSE %.4f\nRMSE %.4f\nsMAPE %.4f', SSEn, RMSEn, sMAPEn),...
    'LineStyle','none',...
    'FitBoxToText','on');

% sum of concentration for sequencing
c = sum(M);