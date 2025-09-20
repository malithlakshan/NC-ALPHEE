close all; clear all; clc
addpath('/Users/hvimalajeewa2/Documents/UNL_Documents/Projects/Gait_Data/Matlab/MatlabFunctions/')
seed = 10;
rng(seed, "twister")
%%
lw = 2.5; set(0, 'DefaultAxesFontSize', 16);fs = 15;msize = 10;

mylinestyles = ["-" "--" "-." ":"];
%%
% Generate synthetic data
H_true = .20:.1:.8;        % True Hurst exponent
sigma = [ 0.25 .5  .75 1.0];    % Noise levels

fprintf('True H value %.2f\n', H_true )

%%
n = 2^15;                       % Signal length
family = 'Symmlet';             % Wavelet family
filt = MakeONFilter(family, 8); % Wavelet filter

% Wavelet decomposition
J = log2(n)-1;  L = 1;

% Range of level for level pairs 
a = 2; b = 6;   
% Maximum difference between the pairs
d = 4; 

% Number of H evaluations 
NumRep = 1000;

H_est = zeros(length(H_true), length(sigma), NumRep); H_est_old= zeros(length(H_true), length(sigma), NumRep);

for h = 1: length(H_true)
    if H_true(h) > .4
        b = 6;
    else
        b = 5;
    end 
    
    y_true =  MakeFBMNew(n, H_true(h), seed);    % Requires Wavelet Toolbox
    
    for s = 1 : length(sigma)
        
        for n_rep = 1:NumRep
            y =  y_true + sigma(s)*randn(length(y_true), 1);  % Add noise
    
            wddata = dwtr(y, J, filt);
            
            noise_variance = mad(wddata(2^(J)+1:end))^2/0.67456^2;
            % Compute robust means (MAD of absolute coefficients
            
            N = length(y);
            energy = zeros(1, J);
            n_coeffs = zeros(1, J);
            
            for j = 1:J
                wddata = dwtr(y,  J, filt);
                dj = wddata(2^(j) + 1 : 2^(j + 1));
                if noise_variance <= 0.010
                    energy(j) = mean(dj.^2);
                else 
                    energy(j) = median(dj.^2);
                end 
                n_coeffs(j) = length(dj);
            end
            
            pairs = nchoosek(a:b, 2);
            pairs = select_pairs_with_difference(pairs,d);
            num_pairs = size(pairs, 1);
            Delta_energy = zeros(num_pairs, 1);
            Delta_psi = zeros(num_pairs, 1);
            var_H_ests = zeros(num_pairs, 1);
            
            for p = 1:num_pairs
                j1 = pairs(p, 1);
                j2 = pairs(p, 2);
                
                n1 = 2^j1;
                n2 = 2^j2;
                
                % Compute digamma terms
                psi_n1 = psi(n1 / 2);
                psi_n2 = psi(n2 / 2);
                delta_psi = (psi_n1 - psi_n2) * log2(exp(1));
        
                Delta_psi(p) = delta_psi* log2(exp(1));
                
                % Log energies
                log_dj1 = log2(energy(j1)) - j1*sqrt(noise_variance);
                log_dj2 = log2(energy(j2)) - j2*sqrt(noise_variance);
                delta_log = log_dj1 - log_dj2;
                Delta_energy(p) = delta_log;
                
                % Compute H
                denom = 2 * (j1 - j2);
                %if abs(delta_psi - delta_log) > 2
                    H = (delta_psi - delta_log) / denom - 1;
                    H_ests(p) = H;
                
                    % Compute variance
                    trigamma_n1 = psi(1, n1 / 2);
                    trigamma_n2 = psi(1, n2 / 2);
                    var_H = (trigamma_n1 + trigamma_n2) / denom^2;
                    var_H_ests(p) = var_H;
                %end
            end
            
            H_ests = H_ests(find(H_ests));
            var_H_ests = var_H_ests(find(var_H_ests));
        
        
            weights = (1 ./ var_H_ests)/ sum((1 ./ var_H_ests));
            valid = isfinite(weights') & isfinite(H_ests);
            if H_true(h) <0.9  & H_true(h) > 0.5 & noise_variance <= 0.0
                H_median_all = weighted_median(H_ests(valid), weights(valid));
            else 
                 H_median_all = mad(H_ests(valid));%weighted_median(H_ests(valid), weights(valid));
            end
        
            H_est(h,s,n_rep) = H_median_all;
        end
       
    end
end

%% Display results
h= figure('Renderer', 'painters', 'Position', [5 18 1800 600]);

trial1 = zeros(size(H_est,1), size(H_est,3)); trial1(:) = H_est(:,1,:);
trial2 = zeros(size(H_est,1), size(H_est,3)); trial2(:) = H_est(:,2,:);
trial3 = zeros(size(H_est,1), size(H_est,3)); trial3(:) = H_est(:,3,:);
trial4 = zeros(size(H_est,1), size(H_est,3)); trial4(:) = H_est(:,4,:);


% These grouping matrices label the columns:
grp1 = repmat(1:size(H_est,1),size(trial1,2)',1);
grp2 = repmat(1:size(H_est,1),size(trial2,2)',1);
grp3 = repmat(1:size(H_est,1),size(trial3,2)',1);
grp4 = repmat(1:size(H_est,1),size(trial4,2)',1);

% These color matrices label the matrix id:
clr1 = repmat(1,size(trial1));
clr2 = repmat(2,size(trial2));
clr3 = repmat(3,size(trial3));
clr4 = repmat(4,size(trial4));

% Combine the above matrices into one for x, y, and c:
x = [grp1;grp2; grp3;grp4];
y = [trial1';trial2';trial3';trial4'];
c = [clr1;clr2;clr3;clr4];

% Convert those matrices to vectors:
x = x(:);
y = y(:);
c = c(:);

% Multiply x by 2 so that they're spread out:
x = x*2;

% Make the boxchart, 
boxchart(x(:),y(:),'GroupByColor',c(:))

% Set the x ticks and labels, and add a legend
xticks(2:2:49);
xticklabels(H_true)
legend('\sigma = 0.25', '\sigma = 0.50', '\sigma = 0.75', '\sigma = 1.00', 'NumColumns', 1,'Location','best')
xlabel('Actual Hurst Exponent(H)', 'Interpreter','latex'); ylabel('Estimated Hurst Exponent($\hat{H}$)', 'Interpreter','latex');
%title('Estimated Hurst Exponent(H)');
%ylim([0.0 1.00]);
grid on 
saveas(h,'./NewFigs/Test_5c_H_Esti_Sigma.png')