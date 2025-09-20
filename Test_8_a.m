close all; clear all; clc
addpath('/Users/hvimalajeewa2/Documents/UNL_Documents/Projects/Gait_Data/Matlab/MatlabFunctions/')
addpath '/Users/hvimalajeewa2/Documents/TAMU_Documents/TAMU/WaveletSrinkage/LPM/matlab'


seed = 10;
rng(seed, "twister")
%%
lw = 2.5; set(0, 'DefaultAxesFontSize', 16);fs = 15;msize = 10;

mylinestyles = ["-" "--" "-." ":"];
%%
n = 2^16;               % Signal length
J = log2(n)-1;  L = 1;
        
family = 'Symmlet';
filt = MakeONFilter(family, 6);

% Generate synthetic data
Hval = [.1 .2 .3 .4 .5 .6 .7 .8];        % True Hurst exponent
%Hval = [.6];        % True Hurst exponent
sigma = [ 0.0 0.25 .5 .75 1 1.5];    % Noise Level


a = 3; b = J-1;
d = b-a; % difference between the level pair  j1 - j2


nrep = 100;

H_est = zeros(length(Hval), length(sigma), nrep);

X = []; 
for i = 1:length(Hval)

    H_true = Hval(i);           % True Hurst exponen
    
    % if H_true > .6 
    %     a=1; b = 3;
    % else
    %    a= 1;b = 10;
    % end 
    
    y_true = MakeFBMNew(n, H_true,seed);
    fprintf('True H value %.2f\n', H_true )
    
    H_sigma =zeros(length(sigma), nrep);
    k = 1;
    for  s = 6%:length(sigma) 
        
        for n_rep = 1:nrep
            %y_true = MakeFBMNew(n, H_true);
            y = y_true + sigma(s)*randn(n, 1);  % Add noise
            
            % Wavelet decomposition
            wddata = dwtr(y, J, filt);
            %if H_true <= 1 
            %ybams= recmybams(y,filt);

            %noise_variance = var(y - ybams)
            noise_variance = var(wddata(2^(J)+1:end));%
            
            %else
            %noise_variance = mad(wddata(2^(J)+1:end))^2/0.67456^2
            %end
    
            % if noise_variance >= .6 
            % 
            %     if H_true <= .6 
            %         b = 9;
            %     else
            %         b = 3;
            %     end
            % end
            % Compute robust means (MAD of absolute coefficients
            
            N = length(y);
            energy = zeros(1, J);
            n_coeffs = zeros(1, J);
            
            for j = 1:J
                wddata = dwtr(y,  J, filt);
                dj = wddata(2^(j) + 1 : 2^(j + 1));
                mn = median(dj);
                if noise_variance <= 0.00
                    energy(j) = mean( (dj- 0).^2);
                else 
                    energy(j) = mean( (dj - 0).^2);
                end 
                n_coeffs(j) = length(dj);
            end
            
            pairs = nchoosek(a:b, 2);
            pairs = select_pairs_with_difference(pairs,d);
            num_pairs = size(pairs, 1);
            Delta_energy = zeros(num_pairs, 1);
            Delta_psi = zeros(num_pairs, 1);
            H_ests = zeros(num_pairs, 1);
            var_H_ests = zeros(num_pairs, 1);
            
            
            for p = 1:num_pairs
                j1 = pairs(p, 1);
                j2 = pairs(p, 2);
                
                n1 = 2^j1;
                n2 = 2^j2;
                
                % Compute digamma terms
                psi_n1 = psi(n1 / 2);
                psi_n2 = psi(n2 / 2);
                delta_psi = (psi_n1 - psi_n2) / log(2);
        
                Delta_psi(p) = delta_psi;
                
                % Log energies
                log_dj1 = n1*energy(j1) - 2*noise_variance*exp(psi_n1);
                log_dj2 = n2*energy(j2) - 2*noise_variance*exp(psi_n2);
                delta_log = log2(log_dj1/log_dj2);
                Delta_energy(p) = delta_log;
                
                % Compute H
                denom = 2 * (j1 - j2);
        
                H =  (delta_psi - delta_log) / denom - 1/2;
                
                H_ests(p) = H;
                
                % Compute variance
                % trigamma_n1 = psi(1, n1 / 2);
                % trigamma_n2 = psi(1, n2 / 2);
                % 
                % term_1 = (n1*energy(j1))/(log_dj1)^2; 
                % term_2 = (n2*energy(j2))/(log_dj2)^2;
                % 
                % var_H = 2*(term_1 *trigamma_n1 + term_2 *trigamma_n2) / (denom*log(2))^2;
                % var_H_ests(p) = var_H;
    
                %fprintf('j1 = %d, j2 = %d, delta_psi =%.2f,delta_log = %.2f, H= %.2f, var H = %.2f \n', j1, j2,delta_psi,delta_log, H, var_H)
               
            end
            
            H_ests = real(H_ests(find(H_ests)));

            X = [X; [H_ests' H_true ]];
            %mean( H_ests)
            % var_H_ests = var_H_ests(find(var_H_ests));
            % 
            % 
            % weights = (1 ./ var_H_ests)/ sum((1 ./ var_H_ests));
            % valid = isfinite(weights) & isfinite(H_ests);
            % 
            % if  H_true < 1 
            %     H_median_all = weighted_median(H_ests(valid), weights(valid));
            % else 
            %      H_median_all = mean( H_ests);
            %      %H_median_all = weighted_median(H_ests(valid), weights(valid));
            % end
            % 
            % H_sigma(s, n_rep) =  H_median_all;
            % H_est(i,s,n_rep) = H_median_all;
        
        end
    k = k+1;
    end    
end
%%
rand_id = randperm(size(X,1),size(X,1));

x = X(rand_id,1:end-1); y = X(rand_id,end);

%%
[numSamples inputSize] = size(x); 
X = x';%rand(inputSize, numSamples); % Inputs (78 features, 1000 samples)
Y = y';%rand(1, numSamples);         % Targets (1 output, 1000 samples)

% Split data into training, validation, and testing
trainRatio = 0.7;
valRatio = 0.15;
testRatio = 0.15;

[trainInd, valInd, testInd] = dividerand(numSamples, trainRatio, valRatio, testRatio);

XTrain = X(:, trainInd);
YTrain = Y(:, trainInd);

XVal = X(:, valInd);
YVal = Y(:, valInd);

XTest = X(:, testInd);
YTest = Y(:, testInd);
%%
hiddenLayerSizes = [50, 40, 30, 20, 10]; % example hidden layers sizes

net = feedforwardnet(hiddenLayerSizes);

% Setup data division explicitly
net.divideFcn = 'divideind';
net.divideParam.trainInd = trainInd;
net.divideParam.valInd = valInd;
net.divideParam.testInd = testInd;

% Set Training function and performance
net.trainFcn = 'trainlm'; % Levenberg-Marquardt (fast for smaller datasets)
net.performFcn = 'mse';   % Mean squared error

% Optionally configure training parameters
net.trainParam.epochs = 500;
net.trainParam.goal = 1e-6;

%%
[net, tr] = train(net, X, Y);

%%
% Predictions
YPredTrain = net(XTrain);
YPredVal   = net(XVal);
YPredTest  = net(XTest);

% Performance
trainPerf = perform(net, YTrain, YPredTrain);
valPerf   = perform(net, YVal, YPredVal);
testPerf  = perform(net, YTest, YPredTest);

fprintf('Training MSE: %.5f\n', trainPerf);
fprintf('Validation MSE: %.5f\n', valPerf);
fprintf('Testing MSE: %.5f\n', testPerf);

%%
figure;
plotperform(tr); % Plot performance during training (training, validation, test)

%%

optimVars = [
    optimizableVariable('Layer1',[10,100],'Type','integer')
    optimizableVariable('Layer2',[5,50],'Type','integer')
    optimizableVariable('Layer3',[5,50],'Type','integer')
    optimizableVariable('Layer4',[5,50],'Type','integer')
    optimizableVariable('Layer5',[5,30],'Type','integer')
    ];

ObjFcn = @(optVars) nnObjectiveFcn(...
    struct('hiddenLayerSizes',[optVars.Layer1,optVars.Layer2,optVars.Layer3,optVars.Layer4,optVars.Layer5]),...
    XTrain,YTrain,XVal,YVal);

results = bayesopt(ObjFcn,optimVars,...
    'MaxObjectiveEvaluations',10,...
    'AcquisitionFunctionName','expected-improvement-plus');

%

bestHiddenLayers = [...
    results.XAtMinEstimatedObjective.Layer1,...
    results.XAtMinEstimatedObjective.Layer2,...
    results.XAtMinEstimatedObjective.Layer3,...
    results.XAtMinEstimatedObjective.Layer4,...
    results.XAtMinEstimatedObjective.Layer5];

%
finalNet = feedforwardnet(bestHiddenLayers);
finalNet.trainFcn = 'trainlm';
finalNet.divideFcn = 'divideind';
finalNet.divideParam.trainInd = trainInd;
finalNet.divideParam.valInd = valInd;
finalNet.divideParam.testInd = testInd;
finalNet = train(finalNet,X,Y);

% Evaluate performance on test data
YPredTest = finalNet(XTest);
testError = perform(finalNet,YTest,YPredTest);
fprintf('Final Test MSE: %.5f\n', testError);

figure
plot(YTest, 'r.-')
hold on; 
plot(YPredTest,'b.-')

%%


% Sample data (replace with your actual vectors)
%YTest = 0.1:0.1:0.8;
%YPredTest = [0.12, 0.09, 0.22, 0.19, 0.31, 0.28, 0.42, 0.39]; % example predictions

% Get unique values from YTest
unique_vals = unique(YTest);

% Calculate mean predictions for each unique value
mean_preds = arrayfun(@(x) mean(YPredTest(YTest == x)), unique_vals);

% Calculate standard deviation of predictions for error bars
std_preds = arrayfun(@(x) std(YPredTest(YTest == x)), unique_vals);

h = figure;
hold on;

% Plot perfect prediction line (y = x)
plot([0 0.9], [0 0.9], 'k--', 'DisplayName', 'Actual H');

% Plot actual vs mean predicted values with error bars
errorbar(unique_vals, mean_preds, std_preds, 'o-', ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', 'b', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Predicted H ± std(H)');

% Customize plot
xlabel('Actual H Values');
ylabel('Predicted H Values ($\hat{H}$)', 'Interpreter','latex');
title('NN Prediction Accuracy');
legend('Location', 'northwest');
grid on;
axis equal;
xlim([0 0.9]);
ylim([0 0.9]);
hold off;
filename = sprintf('./NewFigs/Test_8_a.png');
%saveas(gcf, filename)