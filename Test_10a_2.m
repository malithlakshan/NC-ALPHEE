close all; clear all; clc
%addpath('./MatlabFunctions/')
addpath('/Users/hvimalajeewa2/Documents/UNL_Documents/Projects/Gait_Data/Matlab/MatlabFunctions/')

%% Set parameters for plotting
lw = 2.5; set(0, 'DefaultAxesFontSize', 16);fs = 15;msize = 10;
seed = 10;
rng(seed, "twister")

%% Wavelet parameters 
 
J = 16; n = 2^J; % signal length
L = 1; % decompsition levels 

% wavelet filters 
family = 'Symmlet';
filt = MakeONFilter(family,6); 

% Generate all possible level pairs between 1 and J - L
pairs = nchoosek(1 :J-1, 2); 

% Selected level pairs for the H estimation
a = 3; b = 9;
pairs = pairs(find( pairs(:,1) >= a & pairs(:,2 ) <=b ),:);

% Maximum level pair difference used for the H estimation
D = 1:1:b-a; % difference between the level pair  j1 - j2

% Level-wise energy calculation: 0 - median, 1 - mean
ismean_pair = 1;

% Candidate estimators aggregation method: 0 - weigted median, 1 = weighted mean, 2 - Median, 3 - Mean 
ismean_Hhat = 0; 

% True H values
%Hval     = .10:.1:.8;
H_true = 0.6;

% Noise levels 
Noise_level = [ 0.00 0.25 .5 .75 1]; 

% Number of repetited estimations
nrep = 1000;

% Select the estimation method
Est_Method = ["Standard Wavelet Spectrum", "ALPHEE", "NC-ALPHEE"];

for m = 1:length(Est_Method)
    
    %m = 1; % 1 - standard, 2 - ALPHEE, 3 - NC-ALPHEE
    h= figure('Renderer', 'painters', 'Position', [5 21 2100 500]);
    
    % Estimation of H
    for j = 1: length(D) % for each H value

        H_est  = zeros( length(Noise_level), nrep);
        d = D(j); 

        for s = 1:length(Noise_level) % for each noise level
            y_true = MakeFBMNew(n, H_true,seed);
            
            for i = 1: nrep % for each repetition

                % Generate a fBM with known Hurst exponent H(j)
                if s == 1 
                    y_true = MakeFBMNew(n, H_true);
                end
                
                data = y_true + Noise_level(s)*randn(n, 1); 
                %data = MakeFBMNew(n, H_true) + Noise_level(s)*randn(n, 1); 
                
                % Perform DWT with 15 (= J - L) decompostions; J = 16, L = 1
                wddata = dwtr(data, J - L, filt);
                
                if m == 1
                    % Standard wavelet spectrum based H estimation: 
                    ismean_pair = 0; % mean - 0, median - 1, distance cov - 2
                    isplot = 0;      % 0 - don't plot spectrum, 1 - plot spectrum
                    [slope, levels, log2spec] = waveletspectra_new(data, L, filt, a, b, ismean_pair, 0);
                    H_est(s, i) = (slope * (-1) - 1)/2;
                
                elseif m == 2
                    % ALPHEE based H estimation
                    h_hat = MomentMatchHurst_new(wddata, pairs, L, ismean_Hhat);
                    H_est(s, i) = h_hat;
                
                elseif m == 3
                    % NC-ALPHEE based H estimation when noise variance is zero
                    [H_hat, var_H, Var_pair] = NC_ALPHEE(data, L, filt, a, b, d, ismean_pair, ismean_Hhat);
                    H_est(s,i) = H_hat;

                end  % end of estimation method
            
            end % end of repeated estimation of H
        end % end of noise level
        
        % Plotting the outcomes per method
        subplot(1,6,j);
        groups = repelem(1:length(Noise_level), nrep)';
        data_column = reshape(H_est', [], 1);
        boxchart(groups, data_column); 
        hold on
        yline(H_true,'r--', 'LineWidth', 2.0); grid on
        %ylim([0 1])
        xticks(1:5); % Positions at 1 through 5
        xticklabels(Noise_level); % Custom labels
        
        xlabel('Noise Level ($\sigma$)', 'Interpreter','latex');
        ylabel('Estimated Hurst Exponent($\hat{H}$)', 'Interpreter','latex');
        %yscale('log')
        title(sprintf('|j_1 - j_2| %s %d', char(8804), D(j)))
        grid on

    end % end of H val
    filename = sprintf('./NewFigs/Test_10a_2_%s.png', Est_Method(m));
    saveas(gcf, filename)
end 


