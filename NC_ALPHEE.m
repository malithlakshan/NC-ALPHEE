function [H_hat, var_H, Var_pair] = NC_ALPHEE(signal, L, filt, k1, k2,d, ismean_pair, ismean_Hhat)
    
    
    J = log2(length(signal));
    
    %distance between the level pair
    %d = k2 - k1;

    % Create paris of wavelet decomposition levels
    pairs = nchoosek(k1 :k2, 2);
    pairs = select_pairs_with_difference_new(pairs,d);

    % Perform wavelet transforms on data
    wddata = dwtr(signal, J - L, filt);
    
    W_finest = wddata(2^(J-1)+1:end);
    noise_variance = var( W_finest); % Robust estimate of noise std
    %noise_variance = (mad(W_finest)^2/0.67456^2)^2;

    E = zeros(size(pairs, 1), 2); Var_pair= zeros(size(pairs, 1), 1);
    
    for p = 1: size(pairs,1)
                    
        % select two arbitary levels 
        j1 = pairs(p,1); j2 = pairs(p,2); 
            
        % indexes to select wavelet coefficients 
        j1_indx =  2^(j1) + 1 : 2^(j1 + 1);  n1 = length(j1_indx);
        j2_indx =  2^(j2) + 1 : 2^(j2 + 1);  n2 = length(j2_indx);
                    
       % Log2 of median/ mean of squared wavelet coeffcients at level l1 and l2
        if ismean_pair == 0
           energy_j1 = median( wddata( j1_indx ).^2); 
           energy_j2 = median( wddata( j2_indx ).^2);

        elseif ismean_pair == 1
           energy_j1 = mean( wddata( j1_indx ).^2);  
           energy_j2 = mean( wddata( j2_indx ).^2);

        end 
        
        % Compute digamma terms
        psi_n1 = psi(n1 / 2); % Digamma function for n1/2
        psi_n2 = psi(n2 / 2); % Digamma function for n2/2

        delta_psi = psi_n1 - psi_n2;
        delta_psi = delta_psi / log(2);

        % Compute logarithmic terms
        numerator_j1   = n1*energy_j1 - 2*noise_variance*exp(psi_n1); 
        denominator_j2 = n2*energy_j2 - 2*noise_variance*exp(psi_n2); 
        
        delta_log_dj = log2(numerator_j1/denominator_j2);
        
        denom = 2 * (j1 - j2);
        % Compute H using the original equation
        H = (1 / denom ) * (  delta_psi - delta_log_dj ) - 1/2;
        
        % Compute variance
        % term1 
        trigamma_n1 = psi(1, n1 / 2);
        trigamma_n2 = psi(1, n2 / 2);

        term1 =  trigamma_n1 + trigamma_n2;
        
        % term2
        term2_1 = ( 8*noise_variance^4*exp(2*psi_n1) / (energy_j1^4*(n1 -2)^2*(n1-4)) );
        term2_2 = ( 8*noise_variance^4*exp(2*psi_n2) / (energy_j2^4*(n2 -2)^2*(n2-4)) );
        
        term2 = term2_1 + term2_2;

        % term3
        term3_1 = ( 8*noise_variance^2*exp(psi_n1) / (energy_j1^2*(n1 -2)*n1) ); 
        term3_2 = ( 8*noise_variance^2*exp(psi_n2) / (energy_j2^2*(n2 -2)*n2) ); 
        
        term3 = term3_1 + term3_2;

        v_denom = denom*log(2);

        var_H = (1/v_denom)^2 * (term1 + term2 - term3);
        
        E(p,:) = [real(H) var_H];
        Var_pair(p,:) =   var_H;
   end
    
    weights = ( 1./E(:,2)) ./ sum(1./E(:,2));
    
    % Estimate final Hurst exponent using pairs of levels
    if ismean_Hhat == 3
       h_hat = mean(E(:,1)); % arithmetic median

    elseif ismean_Hhat == 2 
       h_hat = median(E(:,1)); % arithmetic mean 

    elseif  ismean_Hhat == 1                
       h_hat =  [sum(E(:,1).* weights)]; % weighted mean

    elseif ismean_Hhat == 0    
       h_hat = [weighted_median(E(:,1), weights)]; % weighted median

    end
    H_hat = h_hat; var_H = mean(E(:,2)); 
end