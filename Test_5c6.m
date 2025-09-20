close all; clear all; clc
addpath('/Users/hvimalajeewa2/Documents/UNL_Documents/Projects/Gait_Data/Matlab/MatlabFunctions/')
seed = 1;
rng(seed, "twister")
%%
lw = 2.5; set(0, 'DefaultAxesFontSize', 16);fs = 15;msize = 10;

mylinestyles = ["-" "--" "-." ":"];
%%
% Generate synthetic data
Hval = [.1 .2 .3 .4 .5 .6 .7 .8];        % True Hurst exponent
Hval = [.6];        % True Hurst exponent
sigma = [ 0.0 0.25 .5 .75 1];    % Noise Level

n = 2^16;               % Signal length
a = 2; b = 12;
D = 1:1:10; % difference between the level pair  j1 - j2
J = log2(n)-1;  L = 1;
        
family = 'Symmlet';
filt = MakeONFilter(family, 6);

nrep = 10;

H_est = zeros(length(D), length(sigma), nrep);

h= figure('Renderer', 'painters', 'Position', [5 18 1800 600]);

for i = 1:length(D)
    d = D(i);
    H_true = Hval;           % True Hurst exponen
    
    % if H_true > .6 
    %     a=1; b = 10;
    % else
    %    a= 2;b = 10;
    % end 
    
    y_true = MakeFBMNew(n, H_true, seed);
    fprintf('True H value %.2f\n', H_true )
    
    H_sigma =zeros(length(sigma), nrep);
    k = 1;
    for  s = 1:length(sigma)

        for n_rep = 1:nrep
            
            y = y_true + sigma(s)*randn(n, 1);  % Add noise
            
            % Wavelet decomposition
            wddata = dwtr(y, J, filt);
            noise_variance = var(wddata(2^(J)+1:end));%
           
    
            % if noise_variance >= .6 
            % 
            %     if H_true <= .6 
            %         a=2;b = 10;
            %     else
            %         a=2;b = 10;
            %     end
            % end
            
            N = length(y);
            energy = zeros(1, J);
            n_coeffs = zeros(1, J);
            
            for j = 1:J
                wddata = dwtr(y,  J, filt);
                dj = wddata(2^(j) + 1 : 2^(j + 1));

                energy(j) = mean( dj.^2);
 
                n_coeffs(j) = length(dj);
            end
            
            pairs = nchoosek(a:b, 2);
            pairs = select_pairs_with_difference_new(pairs,d);
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
                trigamma_n1 = psi(1, n1 / 2);
                trigamma_n2 = psi(1, n2 / 2);
                
                term_1 = (n1*energy(j1))/(log_dj1)^2; 
                term_2 = (n2*energy(j2))/(log_dj2)^2;
    
                var_H = 2*(term_1 *trigamma_n1 + term_2 *trigamma_n2) / (denom*log(2))^2;
                var_H_ests(p) = var_H;
    
                %fprintf('j1 = %d, j2 = %d, delta_psi =%.2f,delta_log = %.2f, H= %.2f, var H = %.2f \n', j1, j2,delta_psi,delta_log, H, var_H)
               
            end
            
            H_ests = real(H_ests(find(H_ests)));
            %mean( H_ests)
            var_H_ests = var_H_ests(find(var_H_ests));
            %mean(var_H_ests)
            
        
            weights = (1 ./ var_H_ests)/ sum((1 ./ var_H_ests));
            valid = isfinite(weights) & isfinite(H_ests);
    
            if  H_true < 1 
                H_median_all = weighted_median(H_ests(valid), weights(valid));
            else 
                 H_median_all = mean( H_ests);
                 %H_median_all = weighted_median(H_ests(valid), weights(valid));
            end
        
            H_sigma(s, n_rep) =  H_median_all;
            H_est(i,s,n_rep) =  H_median_all;
        
        end
    k = k+1;
    end
    
    % 
    subplot(2,5,i);
    groups = repelem(1:length(sigma), nrep)';
    data_column = reshape(H_sigma', [], 1);
    boxchart(groups, data_column); 
    hold on
    yline(H_true,'r--', 'LineWidth', 2.0); grid on
    ylim([0 1])
    xticks(1:5); % Positions at 1 through 5
    xticklabels(sigma); % Custom labels
    
    xlabel('Noise Level ($\sigma$)', 'Interpreter','latex');
    ylabel('Estimated Hurst Exponent($\hat{H}$)', 'Interpreter','latex');
    %yscale('log')
    title(sprintf('Scale Difference (|j_1 - j_2|) <= %d', D(i)))
    grid on
    % 
    
end
filename = sprintf('./NewFigs/Test_5c6.png');
%saveas(gcf, filename)


%%
h= figure('Renderer', 'painters', 'Position', [5 18 1800 600]);

trial1 = zeros(length(D),nrep); trial1(:) = H_est(:, 1, :); trial1 = trial1';
trial2 = zeros(length(D),nrep); trial2(:) = H_est(:, 2, :); trial2 = trial2';
trial3 = zeros(length(D),nrep); trial3(:) = H_est(:, 3, :); trial3 = trial3';
trial4 = zeros(length(D),nrep); trial4(:) = H_est(:, 4, :); trial4 = trial4';
trial5 = zeros(length(D),nrep); trial5(:) = H_est(:, 5, :); trial5 = trial5';

% These grouping matrices label the columns:
grp1 = repmat(1:length(D),size(trial1,1),1);
grp2 = repmat(1:length(D),size(trial2,1),1);
grp3 = repmat(1:length(D),size(trial3,1),1);
grp4 = repmat(1:length(D),size(trial4,1),1);
grp5 = repmat(1:length(D),size(trial5,1),1);

% These color matrices label the matrix id:
clr1 = repmat(1,size(trial1));
clr2 = repmat(2,size(trial2));
clr3 = repmat(3,size(trial3));
clr4 = repmat(4,size(trial4));
clr5 = repmat(5,size(trial5));

% Combine the above matrices into one for x, y, and c:
x = [grp1;grp2;grp3;grp4;grp5];
y = [trial1;trial2;trial3;trial4;trial5];
c = [clr1;clr2;clr3;clr4;clr5];

% Convert those matrices to vectors:
x = x(:);
y = y(:);
c = c(:);

% Multiply x by 2 so that they're spread out:
x = x*2;

% Make the boxchart, 
boxchart(x(:),y(:),'GroupByColor',c(:));
hold on
yline(H_true,'r--', 'LineWidth', 2.0); grid on
xticks(2:2:50);
xticklabels(D)
legend('\sigma_{\epsilon}= 0','\sigma_{\epsilon} = 0.25', '\sigma_{\epsilon} = 0.50', '\sigma_{\epsilon} = 0.75', '\sigma_{\epsilon} = 1.00', 'NumColumns', 1,'Location','best')
xlabel('Scale Difference ($|j_1 - j_2|$)', 'Interpreter','latex'); 
ylabel('Estimated Hurst Exponent($\hat{H}$)', 'Interpreter','latex');
%yscale('log')
%title('Estimated Hurst Exponent(H)');
%ylim([0.0 .0015]);
grid on 
filename = sprintf('./NewFigs/Test_5c6_all.png');
%saveas(gcf, filename)

