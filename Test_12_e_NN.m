close all; clear all; clc
addpath('/MatlabFunctions/')
addpath '/matlab'


seed = 10;
rng(seed, "twister")
%%
lw = 2.5; set(0, 'DefaultAxesFontSize', 16);fs = 15;msize = 10;

mylinestyles = ["-" "--" "-." ":"];

%%
% Parameters for Wavelet decomposition
n = 2^16;               % Signal length
J = log2(n);  L = 1;  % Wavelet decomposition levels
        
family = 'Symmlet';               % Wavelet filter family
filt = MakeONFilter(family, 6);   % Wavelet filter with size 6

%%
% Parameters for generating synthetic data signals
Hval = [.1 .2 .3 .4 .5 .6 .7 .8];    % True Hurst exponent
sigma = [ 0.0 0.25 .5 .75 1 1.5];    % Noise Level

%%
% Selection of level pairs for the estimation process
a = 3; b = 15;   % selected min and max levels 
d = b-a;         % difference between the level pair  j1 - j2

%% Number of pairs to be used in ALPHEE and NC-ALPHEE
pairs = nchoosek(a:b, 2);                      % all pairs between level a and b 
pairs = select_pairs_with_difference(pairs,d); % pairs with  d difference
num_pairs = size(pairs, 1); 

%% Paramerters for level-pairwise methods 
%Averaging methods 
isPairMean = 1;  % level-pairwise estimation aggregation: 0 - arithmetic mean, 1 - weighted median
isEnergyMean = 0; % signal energy estimation: 0 - arithmetic mean, 1 - median

%%
% Number of repeated evaluations per H and noise level
nrep = 100;

% Noise level to be tested; 
sig_id = 1;

% Store candidate estimators
H_est = zeros(length(Hval), length(sigma), nrep);

X = []; X_var = [];
for i = 1:length(Hval)
    
     % True Hurst exponent 
    H_true = Hval(i);          
    fprintf('True H value %.2f\n', H_true )

    % Generate fractional Brownian motion with H = H_true(i)
    y_true = MakeFBMNew(n, H_true);
    
    H_sigma =zeros(length(sigma), nrep);
    
    k = 1;
    for  s = sig_id%:length(sigma) 
        
        
        for n_rep = 1:nrep
            
            if sigma(s) == 0 % when noise-free, let the signal change
                isEnergyMean = 1;
                y = MakeFBMNew(n, H_true) + sigma(s)*randn(n, 1);       % Add noise
            else  % when noisy, let only the noise change 

                y = y_true + sigma(s)*randn(n, 1);  % Add noise
            end
            
            % Wavelet decomposition
            wddata = dwtr(y, J-1, filt);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            %%%%%%%      NCALPHEE_Update function    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %noise_variance = var(wddata(2^(J)+1:end));%
            W_finest = wddata(2^(J-1)+1:end);

            Wd = sort(abs(wddata(2^(1)+1:end)));
            n_smallest = round(0.3 * length(Wd));
            noise_coeffs = Wd(1:n_smallest);

            if s == 1 
                noise_sigma =  mad(noise_coeffs) / 0.6745;%
            else
                noise_sigma = std(W_finest); 
            end

            noise_variance = noise_sigma^2;

            % Compute robust means (MAD of absolute coefficients
            N = length(y);
            energy = zeros(1, J-1);
            n_coeffs = zeros(1, J-1);

            for j = 1:J-1
                dj = wddata(2^(j) + 1 : 2^(j + 1));
                mn = median(dj);

                if isEnergyMean == 0
                    energy(j) = mean( (dj- 0).^2);

                elseif isEnergyMean == 1 
                    energy(j) = median( (dj - 0).^2);
                end 
                n_coeffs(j) = length(dj);
            end
            
            H_ests = zeros(num_pairs,1); var_H_ests = zeros(num_pairs,1);
            for p = 1:num_pairs
                j1 = pairs(p, 1); j2 = pairs(p, 2);

                n1 = 2^j1; n2 = 2^j2;

                % Compute digamma terms
                psi_n1 = psi(n1 / 2); 
                psi_n2 = psi(n2 / 2);

                delta_psi = (psi_n1 - psi_n2) / log(2);

                % Log energies
                log_dj1 = n1*energy(j1) - 2*noise_variance*exp(psi_n1);
                log_dj2 = n2*energy(j2) - 2*noise_variance*exp(psi_n2);

                delta_log = log2(log_dj1/log_dj2);

                % Denominator
                denom = 2 * (j1 - j2);

                 % Compute H
                H =  (delta_psi - delta_log) / denom - 1/2;

                H_ests( p) = real(H);

                % Compute variance
                % term1
                trigamma_n1 = psi(1, n1 / 2); trigamma_n2 = psi(1, n2 / 2);
                term1 =  trigamma_n1 + trigamma_n2;

                % term2
                term2_1 = ( 8*noise_variance^2*exp(2*psi_n1) / (energy(j1)^2*(n1 -2)^2*(n1-4)) );
                term2_2 = ( 8*noise_variance^2*exp(2*psi_n2) / (energy(j2)^2*(n2 -2)^2*(n2-4)) );

                term2 = term2_1 + term2_2;

                % term3
                term3_1 = ( 8*noise_variance*exp(psi_n1) / (energy(j1)*(n1 -2)*n1) ); 
                term3_2 = ( 8*noise_variance*exp(psi_n2) / (energy(j2)*(n2 -2)*n2) ); 

                term3 = term3_1 + term3_2;

                % Variance Denominator
                v_denom = denom*log(2);

                % Variance
                var_H = (1/v_denom)^2 * (term1 + term2 + term3);

                var_H_ests(p) = var_H;

            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

            % Collect Candidate estimators for training a NN
            weights = (1 ./ var_H_ests)/ sum((1 ./ var_H_ests));
           X = [X; [H_ests' H_true ]];  X_var = [X_var; weights'];
        end
        
    k = k+1;
    end    
end

%% Neural Network Based Aggregation of Candidate Estimators

%% Suffle the data matrix
rand_id = randperm(size(X,1),size(X,1)); 

x = X(rand_id,1:end-1); y = X(rand_id,end); W =  X_var(rand_id,:);

% %%
% [numSamples inputSize] = size(x); 
% X = x';%rand(inputSize, numSamples); % Inputs (78 features, 1000 samples)
% Y = y';%rand(1, numSamples);         % Targets (1 output, 1000 samples)
% 
% % Split data into training, validation, and testing
% trainRatio = 0.7;
% valRatio = 0.15;
% testRatio = 0.15;
% 
% [trainInd, valInd, testInd] = dividerand(numSamples, trainRatio, valRatio, testRatio);
% 
% XTrain = X(:, trainInd);
% YTrain = Y(:, trainInd);
% 
% XVal = X(:, valInd);
% YVal = Y(:, valInd);
% 
% XTest = X(:, testInd);
% YTest = Y(:, testInd);
% 
% %% NN model 
% % tune the number of nodes per hidden layer
% optimVars = [
%     optimizableVariable('Layer1',[10,100],'Type','integer')
%     optimizableVariable('Layer2',[5,50],'Type','integer')
%     optimizableVariable('Layer3',[5,50],'Type','integer')
%     optimizableVariable('Layer4',[5,50],'Type','integer')
%     optimizableVariable('Layer5',[5,30],'Type','integer')
%     ];
% 
% ObjFcn = @(optVars) nnObjectiveFcn(...
%     struct('hiddenLayerSizes',[optVars.Layer1,optVars.Layer2,optVars.Layer3,optVars.Layer4,optVars.Layer5]),...
%     XTrain,YTrain,XVal,YVal);
% 
% % Bayes optimization
% results = bayesopt(ObjFcn,optimVars,...
%     'MaxObjectiveEvaluations',10,...
%     'AcquisitionFunctionName','expected-improvement-plus');
% 
% % Optimized hidden layers 
% bestHiddenLayers = [...
%     results.XAtMinEstimatedObjective.Layer1,...
%     results.XAtMinEstimatedObjective.Layer2,...
%     results.XAtMinEstimatedObjective.Layer3,...
%     results.XAtMinEstimatedObjective.Layer4,...
%     results.XAtMinEstimatedObjective.Layer5];
% 
% % Define the NN
% finalNet = feedforwardnet(bestHiddenLayers);
% finalNet.trainFcn = 'trainlm';
% finalNet.divideFcn = 'divideind';
% finalNet.divideParam.trainInd = trainInd;
% finalNet.divideParam.valInd = valInd;
% finalNet.divideParam.testInd = testInd;
% 
% %  if size(W, 1) == 800 && size(W, 2) == 78
% %     % Check if the first hidden layer size matches
% %     if bestHiddenLayers(1) == 800
% %         finalNet.IW{1,1} = W(1:bestHiddenLayers(1), :);  % Direct assignment if sizes match
% %     else
% %         % If sizes don't match, you might need to transform W
% %         % Here's one approach - take a subset or interpolate
% %         % This is just an example - adjust according to your needs
% %         %finalNet.IW{1,1} = W(1:bestHiddenLayers(1), :);
% %         warning('Input weight matrix was resized to match first hidden layer dimensions');
% %     end
% % else
% %     error('Weight matrix W must be of size 800x55');
% % end
% 
% 
% finalNet = train(finalNet, X, Y);
% 
% % Evaluate performance on test data
% YPredTest = finalNet(XTest);
% testError = perform(finalNet,YTest,YPredTest);

[YPredTest, testError, testInd]  = NN_H(x,y);

fprintf('Final Test MSE: %.5f\n', testError);

%%
h = figure('Renderer', 'painters', 'Position', [5 18 1800 600]);
YTest = y(testInd);
plot(YTest, 'o', ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', 'b', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Actual H ')
hold on; 
plot(YPredTest,'*', ...
    'MarkerSize', 4, ...
    'MarkerFaceColor', 'r', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Predicted H')
grid on
legend('Location', 'northwest');
ylim([0 1.0])
xlabel('Sample index');
ylabel(' H Values', 'Interpreter','latex');
title('NN Prediction Accuracy');
%axis tight

%% NC-ALPHEE & NN performance 

% Select the data corresponds to test indexes
X = x';
Y = y';      
Xk = X'; 
HxTest = Xk(testInd,:); WxTest = X_var(testInd,:); HyTests = Y(testInd)';

HwEstimates = zeros(1,size(HxTest,1));    % weigted median
WeightedMean_H = zeros(1,size(HxTest,1)); % weighted mean
Median_H = zeros(1,size(HxTest,1));       % arimetic median
Mean_H = zeros(1,size(HxTest,1));         % arithmetic mean

for i = 1 : size(HxTest,1)
    HwEstimates(i) = weighted_median(HxTest(i,:), WxTest(i,:));
    WeightedMean_H(i) =  [sum(HxTest(i,:).* WxTest(i,:))];

    Median_H(i) = median(HxTest(i,:));
    Mean_H(i) = mean(HxTest(i,:));
end

% Get unique values from YTest
unique_vals = unique(YTest);

%% Neural Network
% Calculate mean predictions for each unique value
NN_mean_preds = arrayfun(@(x) mean(YPredTest(YTest == x)), unique_vals);

% Calculate standard deviation of predictions for error bars
NN_std_preds = arrayfun(@(x) std(YPredTest(YTest == x)), unique_vals);

%% Weighted median
% Calculate mean predictions for each unique value
Wmed_mean_preds = arrayfun(@(x) mean(HwEstimates(YTest == x)), unique_vals);

% Calculate standard deviation of predictions for error bars
Wemd_std_preds = arrayfun(@(x) std(HwEstimates(YTest == x)), unique_vals);

%% Weighted mean
% Calculate standard deviation of predictions for error bars
Wemn_mean_preds = arrayfun(@(x) mean(WeightedMean_H(YTest == x)), unique_vals);

% Calculate standard deviation of predictions for error bars
Wemn_std_preds = arrayfun(@(x) std(WeightedMean_H(YTest == x)), unique_vals);

%% Arimthetic median
% Calculate standard deviation of predictions for error bars
mdn_mean_preds = arrayfun(@(x) mean(Median_H(YTest == x)), unique_vals);

% Calculate standard deviation of predictions for error bars
mdn_std_preds = arrayfun(@(x) std(Median_H(YTest == x)), unique_vals);

%% Arimthetic mean
% Calculate standard deviation of predictions for error bars
mn_mean_preds = arrayfun(@(x) mean(Mean_H(YTest == x)), unique_vals);

% Calculate standard deviation of predictions for error bars
mn_std_preds = arrayfun(@(x) std(Mean_H(YTest == x)), unique_vals);

%% 
h = figure('Renderer', 'painters', 'Position', [5 18 1800 600]);

hold on;

% Plot perfect prediction line (y = x)
plot([0 0.9], [0 0.9], 'k--', 'DisplayName', 'Actual H');

% NN: Plot actual vs mean predicted values with error bars
errorbar(unique_vals, NN_mean_preds, NN_std_preds, 'o-', ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', 'b', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'NN');

% % Weighted Median: Plot actual vs mean predicted values with error bars
% 
% errorbar(unique_vals, Wmed_mean_preds, Wemd_std_preds, 'o-', ...
%     'MarkerSize', 8, ...
%     'MarkerFaceColor', 'r', ...
%     'LineWidth', 1.5, ...
%     'DisplayName', 'Weighted Median');
% 
% % Weighted Mean: Plot actual vs mean predicted values with error bars
% 
% errorbar(unique_vals, Wemn_mean_preds, Wemn_std_preds, 'o-', ...
%     'MarkerSize', 8, ...
%     'MarkerFaceColor', 'g', ...
%     'LineWidth', 1.5, ...
%     'DisplayName', 'Weighted Mean');
% 
% % Median: Plot actual vs mean predicted values with error bars
% 
% errorbar(unique_vals, mdn_mean_preds, mdn_std_preds, 'o-', ...
%     'MarkerSize', 8, ...
%     'MarkerFaceColor', 'y', ...
%     'LineWidth', 1.5, ...
%     'DisplayName', 'Arithmatic Median');
% 
% % Median: Plot actual vs mean predicted values with error bars
% 
% errorbar(unique_vals, mn_mean_preds, mn_std_preds, 'o-', ...
%     'MarkerSize', 8, ...
%     'MarkerFaceColor', 'k', ...
%     'LineWidth', 1.5, ...
%     'DisplayName', 'Arithmatic Mean');
% 
% 
% % Customize plot
xlabel('Actual Hurst Exponent Values (H)');
ylabel('Predicted H Values ($\hat{H}$)', 'Interpreter','latex');
if sig_id == 1
    title(' Estimations in Noise-free Conditions');
else
    title(' Estimation Accuracy in Noisy Conditions');
end

legend('Location', 'northwest');
grid on;
%axis equal;
xlim([0 0.9]);
%ylim([0 0.9]);
hold off;
if sig_id == 1
    filename = sprintf('./NewFigs/Test_12_e_NN_NoNoise.png');
else
   filename = sprintf('./NewFigs/Test_12_e_NN_Noise.png');
end
%saveas(gcf, filename)
