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
pairs_all = nchoosek(1 :J-1, 2); 

% Selected level pairs for the H estimation
a = 3; b = 13;
pairs = pairs_all(find( pairs_all(:,1) >= a & pairs_all(:,2 ) <=b ),:);

% Maximum level pair difference used for the H estimation
d = b-a;

% Level-wise energy calculation: 0 - median, 1 - mean
ismean_pair = 1;

% Candidate estimators aggregation method: 0 - weigted median, 1 = weighted mean, 2 - Median, 3 - Mean 
ismean_Hhat = 0; 

% True H values
H     = .10:.1:.8;

% Number of repetited estimations
nrep = 1000;

H_Standard = zeros(length(H), nrep);
H_ALPHEE   = zeros(length(H), nrep); 
H_NCALPHEE = zeros(length(H), nrep);

for j = 1: length(H)
    %data_orig = MakeFBMNew(2^J, H(j),seed);
    for i = 1: nrep
        % Generate a fBM with known Hurst exponent H(j)
        % fBm is generated without fixing the seed to allow radnomness as
        % there is no random noise

        data = MakeFBMNew(2^J, H(j)) + 0.0*randn(n, 1); 
        %data = data_orig + 0.0*randn(n, 1);

        % Perform DWT with 14 decompostions; J = 15, L = 1
        wddata = dwtr(data, J - L, filt);
        
        %Standard wavelet spectrum based H estimation
        ismean_pair = 0; % mean - 0, median - 1, distance cov - 2
        [slope, levels, log2spec] = waveletspectra_new(data, L, filt, a, b, ismean_pair, 0);
        H_Standard(j, i) = (slope * (-1) - 1)/2;
        
        %ALPHEE based H estimation
        h_hat = MomentMatchHurst_new(wddata, pairs, L, ismean_Hhat);
        H_ALPHEE(j, i) = h_hat;
        
        % NC-ALPHEE based H estimation when noise variance is zero
        [H_hat, var_H, Var_pair] = NC_ALPHEE(data, L, filt, a, b,d, ismean_pair, ismean_Hhat);
        H_NCALPHEE(j,i) = H_hat;
    
    end
end 


%%
h= figure('Renderer', 'painters', 'Position', [5 12 1200 600]);

trial1 = H_Standard';
trial2 = H_ALPHEE';
trial3 = H_NCALPHEE';

% These grouping matrices label the columns:
grp1 = repmat(1:size(H_Standard,1),size(trial1,1),1);
grp2 = repmat(1:size(H_ALPHEE,1),size(trial2,1),1);
grp3 = repmat(1:size(H_NCALPHEE,1),size(trial3,1),1);

% These color matrices label the matrix id:
clr1 = repmat(1,size(trial1));
clr2 = repmat(2,size(trial2));
clr3 = repmat(3,size(trial3));

% Combine the above matrices into one for x, y, and c:
x = [grp1;grp2;grp3];
y = [trial1;trial2;trial3];
c = [clr1;clr2;clr3];

% Convert those matrices to vectors:
x = x(:);
y = y(:);
c = c(:);

% Multiply x by 2 so that they're spread out:
x = x*2;

% Make the boxchart, 
boxchart(x(:),y(:),'GroupByColor',c(:))

% Set the x ticks and labels, and add a legend
xticks(2:2:16);
xticklabels(H)
legend('Standard', 'ALPHEE','NC-ALPHEE(\sigma_{\epsilon} = 0)', 'NumColumns', 1,'Location','best')
xlabel('Actual Hurst Exponent(H)', 'Interpreter','latex'); 
ylabel('Estimated Hurst Exponent($\hat{H}$)', 'Interpreter','latex');
ylim([0.0 1.00]);
grid on 
%saveas(h,'./NewFigs/Test_10b_Standard_ALPHEE_NCALPHE_NoNoise.png')



