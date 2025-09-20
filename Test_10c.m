close all; clear all; clc
addpath('/Users/hvimalajeewa2/Documents/UNL_Documents/Projects/Gait_Data/Matlab/MatlabFunctions/')
addpath '/Users/hvimalajeewa2/Documents/TAMU_Documents/TAMU/WaveletSrinkage/LPM/matlab'


seed = 1;
rng(seed, "twister")
%%
lw = 2.5; set(0, 'DefaultAxesFontSize', 16);fs = 15;msize = 10;

mylinestyles = ["-" "--" "-." ":"];
%%
% Generate synthetic data
Hval = [.1 .2 .3 .4 .5 .6 .7 .8];        % True Hurst exponent
%Hval = [.1 .2 .3 .4 .5 .6 ];        % True Hurst exponent
%Hval = [.7];
sigma = [ 0.0 0.25 .5 .75 1];    % Noise Level

n = 2^16;               % Signal length
a = 3; b = 9;
d = 3; % difference between the level pair  j1 - j2
J = log2(n)-1;  L = 1;
        
family = 'Symmlet';
filt = MakeONFilter(family, 6);

nrep = 1000;

% energy calculation: 0 - median, 1 - mean
ismean_pair = 1; 

% Aggregation of candidate estimators: 0 - w_median, 1 - w_mean, 2 - median, 3 - mean 
ismean_hat = 0; 

H_est = zeros(length(Hval), length(sigma), nrep);

h= figure('Renderer', 'painters', 'Position', [5 18 1800 600]);

for i = 1:length(Hval)

    H_true = Hval(i);           % True Hurst exponen
    
    %y_true = MakeFBMNew(n, H_true,seed);
    fprintf('True H value %.2f\n', H_true )
    
    H_sigma =zeros(length(sigma), nrep);
    k = 1;
    for  s = 1:length(sigma)
        
        for n_rep = 1:nrep
            
            y =  MakeFBMNew(n, H_true) + sigma(s)*randn(n, 1);  % Add noise
            
            % Wavelet decomposition
            wddata = dwtr(y, J, filt);
            
            % Compute noise variance using finest level of details
            W_finest = wddata(2^(J-1)+1:end);
            noise_variance = var(wddata(2^(J)+1:end));%
            %noise_variance = (mad(W_finest)^2/0.67456^2)^2;

            % Compute wavelet energy
            energy = zeros(1, J); n_coeffs = zeros(1, J);
            
            for j = 1:J
                wddata = dwtr(y, J, filt);
                dj = wddata(2^(j) + 1 : 2^(j + 1));
                
                if ismean_pair == 1
                    energy(j) = mean( dj.^2);

                elseif ismean_pair == 0 
                    energy(j) = median( dj.^2);
                end 

                n_coeffs(j) = length(dj);
            end
            
            % Generate a set of level-pairs between a and b, 1 <= a <= b <= J
            pairs = nchoosek(a:b, 2);

            % Select level-pairs of which the distance is d, |j1-j2| <= d
            pairs = select_pairs_with_difference(pairs,d);
           
            % Estimate H using NC-ALPHEE
            num_pairs = size(pairs,1);
            H_ests = zeros(num_pairs, 1); var_H_ests = zeros(num_pairs, 1);
            
            % Compute H for each selected pair
            for p = 1:num_pairs
                % Levels of the selected pair
                j1 = pairs(p, 1); j2 = pairs(p, 2);
                
                % Number of coefficients in the selected levels
                n1 = 2^j1; n2 = 2^j2;
                
                % Compute digamma terms
                psi_n1 = psi(n1 / 2);
                psi_n2 = psi(n2 / 2);
                delta_psi = log2(exp(1))*(psi_n1 - psi_n2);
                
                % Log energies
                log_dj1 = n1*energy(j1) - 2*noise_variance*exp(psi_n1);
                log_dj2 = n2*energy(j2) - 2*noise_variance*exp(psi_n2);
                delta_log = log2(log_dj1/log_dj2);
                
                % Denominator
                denom = 2 * (j1 - j2);
                
                % Compute H 
                H =  (delta_psi - delta_log) / denom - 1/2;
                
                H_ests(p) = H;
                
                % Compute variance
                % term1
                trigamma_n1 = psi(1, n1 / 2); trigamma_n2 = psi(1, n2 / 2);
                term1 =  trigamma_n1 + trigamma_n2;

                % term2
                term2_1 = ( 8*noise_variance^4*exp(2*psi_n1) / (energy(j1)^4*(n1 -2)^2*(n1-4)) );
                term2_2 = ( 8*noise_variance^4*exp(2*psi_n2) / (energy(j2)^4*(n2 -2)^2*(n2-4)) );

                term2 = term2_1 + term2_2;

                % term3
                term3_1 = ( 8*noise_variance^2*exp(psi_n1) / (energy(j1)^2*(n1 -2)*n1) ); 
                term3_2 = ( 8*noise_variance^2*exp(psi_n2) / (energy(j2)^2*(n2 -2)*n2) ); 

                term3 = term3_1 + term3_2;

                v_denom = denom*log(2);

                var_H = (1/v_denom)^2 * (term1 + term2 - term3);

                var_H_ests(p) = var_H;
                
               %fprintf('j1 = %d, j2 = %d, delta_psi =%.2f,delta_log = %.2f, H= %.2f, var H = %.2f \n', j1, j2,delta_psi,delta_log, H, var_H);
                 
            end
            
            H_ests = real(H_ests(find(H_ests))); var_H_ests = var_H_ests(find(var_H_ests));
            
            % Aggregation of candidate estimators
            weights = (1 ./ var_H_ests)/ sum((1 ./ var_H_ests));
            valid = isfinite(weights) & isfinite(H_ests);

            if  ismean_hat == 0
                H_hat = weighted_median(H_ests(valid), weights(valid));
    
            elseif  ismean_hat == 1                
                H_hat =  [sum(H_ests(valid).* weights(valid))]; 
            
            elseif ismean_hat == 2
                H_hat = median( H_ests(valid));

            elseif ismean_hat == 3
                H_hat = mean( H_ests(valid));
            end
            H_hat;
            H_sigma(s, n_rep) =  H_hat;
            H_est(i,s,n_rep) = H_hat;
        
        end
    k = k+1;
    end
    
    % 
    subplot(2,4,i);
    groups = repelem(1:length(sigma), nrep)';
    data_column = reshape(H_sigma', [], 1);
    boxchart(groups, data_column); 
    hold on
    yline(H_true,'r--', 'LineWidth', 2.0); grid on
    ylim([-1.5 1.5])
    xticks(1:5); % Positions at 1 through 5
    xticklabels(sigma); % Custom labels

    xlabel('Noise Level ($\sigma$)', 'Interpreter','latex');
    ylabel('Estimated Hurst Exponent($\hat{H}$)', 'Interpreter','latex');

    title(sprintf('H = %.2f', H_true))
    % 
    
end
filename = sprintf('./NewFigs/Test_10c_NCALPHEE.png');
saveas(gcf, filename)


%%
% h= figure('Renderer', 'painters', 'Position', [5 18 1800 600]);
% 
% trial1 = zeros(length(Hval),nrep); trial1(:) = H_est(:, 1, :); trial1 = trial1';
% trial2 = zeros(length(Hval),nrep); trial2(:) = H_est(:, 2, :); trial2 = trial2';
% trial3 = zeros(length(Hval),nrep); trial3(:) = H_est(:, 3, :); trial3 = trial3';
% trial4 = zeros(length(Hval),nrep); trial4(:) = H_est(:, 4, :); trial4 = trial4';
% trial5 = zeros(length(Hval),nrep); trial5(:) = H_est(:, 5, :); trial5 = trial5';
% 
% % These grouping matrices label the columns:
% grp1 = repmat(1:length(Hval),size(trial1,1),1);
% grp2 = repmat(1:length(Hval),size(trial2,1),1);
% grp3 = repmat(1:length(Hval),size(trial3,1),1);
% grp4 = repmat(1:length(Hval),size(trial4,1),1);
% grp5 = repmat(1:length(Hval),size(trial5,1),1);
% 
% % These color matrices label the matrix id:
% clr1 = repmat(1,size(trial1));
% clr2 = repmat(2,size(trial2));
% clr3 = repmat(3,size(trial3));
% clr4 = repmat(4,size(trial4));
% clr5 = repmat(5,size(trial5));
% 
% % Combine the above matrices into one for x, y, and c:
% x = [grp1;grp2;grp3;grp4;grp5];
% y = [trial1;trial2;trial3;trial4;trial5];
% c = [clr1;clr2;clr3;clr4;clr5];
% 
% % Convert those matrices to vectors:
% x = x(:);
% y = y(:);
% c = c(:);
% 
% % Multiply x by 2 so that they're spread out:
% x = x*2;
% 
% % Make the boxchart, 
% boxchart(x(:),y(:),'GroupByColor',c(:))
% 
% % Set the x ticks and labels, and add a legend
% xticks(2:2:50);
% xticklabels(Hval)
% legend('\sigma_{\epsilon}=0','\sigma_{\epsilon} = 0.25', '\sigma_{\epsilon} = 0.50', '\sigma_{\epsilon} = 0.75', '\sigma_{\epsilon} = 1.00', 'NumColumns', 1,'Location','best')
% xlabel('Actual Hurst Exponent(H)', 'Interpreter','latex'); 
% ylabel('Estimated Hurst Exponent($\hat{H}$)', 'Interpreter','latex');
% %title('Estimated Hurst Exponent(H)');
% ylim([0.0 1.00]);
% grid on 
% filename = sprintf('./NewFigs/Test_5c3_NCALPHEE_all.png');
% %saveas(gcf, filename)
% 
