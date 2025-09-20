# Neural Network-based Noise Resilient  Approach for Robust Hurst Exponent Estimation
We propose **Noise-Controlled ALPHEE (NC-ALPHEE)**, an enhancement of the **Average Level-Pairwise Hurst Exponent Estimator (ALPHEE)**, incorporating noise mitigation and generating multiple level-pairwise estimates from signal energy pairs. A neural network (NN) combines these estimates, replacing traditional averaging. This adaptive learning maintains ALPHEE’s behavior in noise-free cases while improving performance in noisy conditions. The steps below explain how to implement the **NC-ALPHEE** using **Matlab** and continue the Aggregation using **Python**.

### Matlab Codes
The repository includes Matlab files that are used to implement
  1. The estimation of Hurst exponents from fractional Brownian motion using wavelet-based methods and compares traditional aggregation with a neural network approach.
  2. The simulated fBM across H and noise levels, performs wavelet decompositions to compute pairwise H estimators with and without noise correction (m=1 vs m=2), derives variance-based weights, and exports CSVs for downstream neural-network aggregation and evaluation.
  3. The Hurst exponent via DWT-based, noise-corrected level-pair energy (with digamma/trigamma adjustments) and aggregates pairwise estimates (weighted median/mean or median/mean), returning estimated H its and variance.
  4. the respective figures included in the paper

The MatlabFunctions folder contains a set of functions used in the Matlab files.
In the following, a brief introduction, for each code, is provided to explain its functionality.

1. **Test_12_e_NN.m** : This MATLAB script simulates fractional Brownian motion signals with varying Hurst exponents and noise levels, applies wavelet-based decomposition, and computes multiple candidate H estimators from scale-pair energy statistics. It then aggregates these estimators using both traditional methods (weighted mean/median) and a neural network model to improve robustness. The script evaluates prediction accuracy across methods, visualizes actual vs. estimated Hurst exponents, and compares performance in noise-free and noisy conditions.
2. **Test_13_b.m** : This MATLAB script generates fractional Brownian motion signals across a range of Hurst exponents and noise levels, applies wavelet decomposition, and computes candidate H estimators using both noise-corrected and uncorrected formulations. For each estimator, it calculates variance-based weights to quantify reliability and stores both the estimates and weights in CSV files for later analysis. The framework is designed to support advanced aggregation methods, such as neural network models, to combine multiple estimators and improve accuracy in both noise-free and noisy conditions.
3. **NC_ALPHEE.m** : This function implements the NC-ALPHEE method to estimate the Hurst exponent of a signal using wavelet decomposition. After decomposing the signal into detail coefficients across multiple levels, it forms level pairs within a specified range and fixed gap. For each pair, it computes scale-wise energies (using either the mean or median of squared coefficients), applies noise correction based on the finest scale variance, and incorporates digamma/trigamma terms to adjust for distributional effects. Each pair produces a candidate Hurst estimate and its variance, which are then combined across all pairs using either weighted or unweighted aggregation methods (weighted median, weighted mean, arithmetic mean, or median). The function outputs the final Hurst estimate, the mean variance across pairs, and the individual variances for each pair.


### Python Codes
The repository includes **.ipynb** files that are used to implement
  1. Preprocess the dataset (scaling, splitting, tensor conversion).
  2. Define and train a tunable feedforward neural network (PyTorch).
  3. Optimize hyperparameters with Optuna (K-Fold CV + early stopping).
  4. Plot results/metrics to visualize model performance.

