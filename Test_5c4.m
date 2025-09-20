close all; clear all; clc
addpath('/Users/hvimalajeewa2/Documents/UNL_Documents/Projects/Gait_Data/Matlab/MatlabFunctions/')
seed = 1;
rng(seed, "twister")

%% plotting parameters
lw = 2.5; set(0, 'DefaultAxesFontSize', 16);fs = 15;msize = 10;

mylinestyles = ["-" "--" "-." ":"];
%%
% Generate synthetic data
Hval = [.1 .2 .3 .4 .5 .6 .7 .8];  % True Hurst exponent
%Hval = [.6];        
sigma = [ 0.0 0.25 .5 .75 1];      % Noise Levels

n = 2^16;               % Signal length
a = 2; b = 10;          % Lower and Upper bounds of levels j1 & j2
d = 10;                 % Selected difference between the level pair  j1 - j2
ismean_pair = 1;        % method to compute wavelet energy
ismean_Hhat = 0;        % Method to aggregate candidate estimators  
J = log2(n)-1;  L = 1;  % Maximum number of decompositions J
        
family = 'Symmlet';               % Wavelet filter family
filt = MakeONFilter(family, 8);   % Wavelet filter

nrep = 20;  % Number of repeted estimations

%% ALPHEE method 
H_est = zeros(length(Hval), length(sigma), nrep); % Storing H estimations at differet noise levels

h= figure('Renderer', 'painters', 'Position', [5 18 1800 600]);

for i = 1:length(Hval)

    H_true = Hval(i);           % True Hurst exponent
    
    % if H_true > .6 
    %     a=1; b = 3;
    % else
    %    a= 1;b = 10;
    % end 
    
    %y_true = MakeFBMNew(n, H_true,seed); % Generate a fBm with H_true
    fprintf('True H value %.2f\n', H_true )
    
    H_sigma =zeros(length(sigma), nrep);
    k = 1;
    for  s = 1:length(sigma)
        
        for n_rep = 1:nrep
            y_true = MakeFBMNew(n, H_true); % Generate a fBm with H_true
            y = y_true + sigma(s)*randn(n, 1);  % Noisy signal: Add noise to the fBm, y_true 
            
             wddata = dwtr(y, J, filt);         % Perform DWT on the noisy signal, y 
            noise_variance = var(wddata(2^(J)+1:end)); % estimated noise variance using finest level of detalis
    
            if noise_variance >= .6 

                if H_true <= .6 
                    b = 9;
                else
                    b = 3;
                end
            end
            % Estimate H using ALPHEE method. 
            [H_hat, var_H, Energ_pair, Var_pair] = ALPHEE(y, L, filt, a, b,ismean_pair , ismean_Hhat,d);
        
            H_sigma(s, n_rep) =  H_hat;
            H_est(i,s,n_rep) = H_hat;
        
        end
    k = k+1;
    end
    
    % box plot per each H true H value at different noise levels
    
    subplot(2,4,i)
    groups = repelem(1:length(sigma), nrep)';
    data_column = reshape(H_sigma', [], 1);
    boxchart(groups, data_column); 
    hold on
    yline(H_true,'r--', 'LineWidth', 2.0); grid on
    ylim([-1 1])
    xticks(1:5); % Positions at 1 through 5
    xticklabels(sigma); % Custom labels
    
    xlabel('Noise Level ($\sigma_{\epsilon}$)', 'Interpreter','latex');
    ylabel('Estimated Hurst Exponent($\hat{H}$)', 'Interpreter','latex');
    
    title(sprintf('True H value (H) = %.2f', H_true))
    % 
    
end
filename = sprintf('./NewFigs/Test_5c4_ALPHEE.png');
%saveas(gcf, filename)

%% 
h= figure('Renderer', 'painters', 'Position', [5 18 1800 600]);

% Define parameters
% Sample Data:

trial1 = zeros(nrep,length(Hval));


trial1 = zeros(length(Hval),nrep); trial1(:) = H_est(:, 1, :); trial1 = trial1';
trial2 = zeros(length(Hval),nrep); trial2(:) = H_est(:, 2, :); trial2 = trial2';
trial3 = zeros(length(Hval),nrep); trial3(:) = H_est(:, 3, :); trial3 = trial3';
trial4 = zeros(length(Hval),nrep); trial4(:) = H_est(:, 4, :); trial4 = trial4';
trial5 = zeros(length(Hval),nrep); trial5(:) = H_est(:, 5, :); trial5 = trial5';

% These grouping matrices label the columns:
grp1 = repmat(1:length(Hval),size(trial1,1),1);
grp2 = repmat(1:length(Hval),size(trial2,1),1);
grp3 = repmat(1:length(Hval),size(trial3,1),1);
grp4 = repmat(1:length(Hval),size(trial4,1),1);
grp5 = repmat(1:length(Hval),size(trial5,1),1);

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
boxchart(x(:),y(:),'GroupByColor',c(:))

% Set the x ticks and labels, and add a legend
xticks(2:2:50);
xticklabels(Hval)
legend('\sigma_{\epsilon}= 0.00','\sigma_{\epsilon} = 0.25', '\sigma_{\epsilon} = 0.50', '\sigma_{\epsilon} = 0.75', '\sigma_{\epsilon} = 1.00', 'NumColumns', 1,'Location','best')
xlabel('Actual Hurst Exponent(H)', 'Interpreter','latex'); ylabel('Estimated Hurst Exponent($\hat{H}$)', 'Interpreter','latex');
%title('Estimated Hurst Exponent(H)');
%ylim([0.0 1.00]);
grid on 

filename = sprintf('./NewFigs/Test_5c4_ALPHEE_all.png');
%saveas(gcf, filename)

