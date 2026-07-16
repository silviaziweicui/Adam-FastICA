%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adam-FastICA Online Decomposition Demo
% 
% This script demonstrates the offline and online decomposition of 
% high-density surface EMG (HDsEMG) signals using the Adam-FastICA algorithm.
% 
% Workflow:
%   1. Load and preprocess EMG data
%   2. Perform offline Adam-FastICA decomposition to obtain initial parameters
%   3. Perform online decomposition using the pre-trained parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
FS = 2048;
win_size_ms = 200;

%Please load data first

%% adam+FastICA offline (requires train_data)
[signal_a, parameters_a, signalprocess_a] = adam_fastica(200, train_data, FS);

%% adam+FastICA online (requires test_data and signalprocess_a)
online3 = struct();
[online3.ipt, online3.spikes] = online_decomp(test_data, win_size_ms, signalprocess_a, FS); 


                    







