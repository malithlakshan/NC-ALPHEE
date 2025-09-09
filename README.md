# Neural Network-based Noise Resilient  Approach for Robust Hurst Exponent Estimation
We propose **Noise-Controlled ALPHEE (NC-ALPHEE)**, an enhancement of the **Average Level-Pairwise Hurst Exponent Estimator (ALPHEE)**, incorporating noise mitigation and generating multiple level-pairwise estimates from signal energy pairs. A neural network (NN) combines these estimates, replacing traditional averaging. This adaptive learning maintains ALPHEE’s behavior in noise-free cases while improving performance in noisy conditions. The steps below explain how to implement the **NC-ALPHEE** using **Matlab** and continue the Aggregation using **Python**.

### Matlab Codes
  1. first
  2. second
  3. Third

The MatlabFunctions folder contains a set of functions used in the Matlab files.
In the following, a brief introduction, for each code, is provided to explain its functionality.

1. Implementation of the new method
2. Simulation Study


### Python Codes
The repository includes **.ipynb** files that are used to implement
  1. Preprocess the dataset (scaling, splitting, tensor conversion).
  2. Define and train a tunable feedforward neural network (PyTorch).
  3. Optimize hyperparameters with Optuna (K-Fold CV + early stopping).
  4. Plot results/metrics to visualize model performance.

