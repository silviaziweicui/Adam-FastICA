%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adam-FastICA 
% 
% This function performs offline decomposition of high-density surface EMG 
% (HDsEMG) signals using the Adam-FastICA algorithm. 
%
% Key Steps:
%   1. Preprocessing (Notch filter, Bandpass filter, Signal extension)
%   2. Whitening (PCA-based with regularization)
%   3. FastICA with Adam optimization
%   4. Coefficient of Variation (CoV) minimization
%   5. SIL-based filtering and duplicate removal
%   6. Batch processing to extract pulse trains and centroids
%
% Input: 
%   niter: number of iterations (number of motor units to extract)
%   data: input sEMG data (channels x samples)
%   FS: sampling frequency (Hz)
%
% Output:
%   signal: struct containing pulse trains and discharge times
%       - signal.Pulsetrain: pulse trains of each MU
%       - signal.Dischargetimes: discharge times of each MU
%   parameters: struct containing decomposition parameters
%   signalprocess: struct containing intermediate processing results
%       - signalprocess.MUFilters: motor unit filters
%       - signalprocess.wSIG: whitened signal
%       - signalprocess.centroids: cluster centroids
%       - signalprocess.whiteningMatrix: whitening matrix
%       - signalprocess.CorrSig_all: accumulated correlation matrix
%       - signalprocess.eSigLen: length of extended signal
%       - signalprocess.Ctx: Ctx matrix for online update
%       - signalprocess.exFactor: extension factor
%   K: empty placeholder (for compatibility)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [signal, parameters, signalprocess] = adam_fastica(niter, data, FS)
% DECOMPOSITION PARAMETERS
parameters.NITER = niter; % total iterations
parameters.differentialmode = 0; % 0 = no; 1 = yes (filter out the smallest MU, can improve decomposition at the highest intensities
parameters.initialization = 0; % 0 = max EMG; 1 = random weights
parameters.covfilter = 1; % 0 = no; 1 = yes (filter out the motor units with a coefficient of variation of their ISI > than parameters.covthr)

% SPECIFIC VALUES
parameters.nbextchan = 256; % nb of extended channels 
parameters.contrastfunc = 'logcosh'; % contrast functions: 'skew', 'kurtosis', 'logcosh'
parameters.silthr = 0.6; 
parameters.covthr = 0.6; 
parameters.duplicatesthresh = 0.3; % threshold that define the minimal percentage of common discharge times between duplicated motor units

%% Step 0: Load the HDsEMG data
signal.data = data;
signal.fsamp = FS; 
signal.emgtype = 1;  
signal.EMGmask = zeros(size(signal.data,1),1);

signalprocess.data = signal.data;

signalprocess.coordinatesplateau = [1 size(signalprocess.data,2)];

%% 
% Step 1: Preprocessing
%       1a: Removing line interference (Notch filter)
signalprocess.data = notchsignals(signalprocess.data,signal.fsamp);
%       1b: Bandpass filtering
signalprocess.data = bandpassingals(signalprocess.data, signal.fsamp, signal.emgtype);
%       1c: Differentiation (perform only if there is many motor units,
%      filter out the smallest motor units) useful for high intensities
if parameters.differentialmode == 1
        signalprocess.data = diff(signalprocess.data,1,2);
end
%       1d: Signal extension (extension factor calculated to reach 1000
%       channels)
signalprocess.exFactor = round(parameters.nbextchan/size(signalprocess.data,1));
signalprocess.ReSIG = zeros(signalprocess.exFactor * size(signalprocess.data,1));
signalprocess.iReSIG = zeros(signalprocess.exFactor * size(signalprocess.data,1));

signalprocess.eSIG = extend(signalprocess.data,signalprocess.exFactor);
signalprocess.CorrSig_all = signalprocess.eSIG * signalprocess.eSIG'; 
signalprocess.ReSIG = signalprocess.eSIG * signalprocess.eSIG' / size(signalprocess.eSIG,2);
signalprocess.iReSIG = pinv(signalprocess.ReSIG);
signalprocess.eSigLen = size(signalprocess.eSIG,2);

%       1e: Removing the mean
signalprocess.eSIG = demean(signalprocess.eSIG);

%% -----------------------------------------------------------------%
% Step 2: Whitening
%Whitening with a regularization factor (average of the smallest half of 
%the eigenvalues of the covariance matrix from the extended signals)

%       2a: Get eigenvalues and eigenvectors (regularization factor =>
%       average smallest half of eigenvalues)
[E, D] = pcaesig(signalprocess.eSIG); %Returns the eigenvector (E) and diagonal eigenvalue (D) matrices

%       2b: Zero-phase component analysis
[signalprocess.wSIG, signalprocess.whiteningMatrix, ~] = whiteesig(signalprocess.eSIG, E, D);
clearvars E D

%% -----------------------------------------------------------------%
% Step 3: FastICA method

% Initialize matrix B (n x m) n: separation vectors, m: iterations 
% Initialize matrix MUFilters to only save the reliable filters
% Intialize SIL and PNR

signalprocess.B = zeros(size(signalprocess.wSIG,1), parameters.NITER); % all separation vectors
signalprocess.MUFilters = zeros(size(signalprocess.wSIG,1), parameters.NITER); % only reliable vectors
signalprocess.w = zeros(size(signalprocess.wSIG,1), 1);
signalprocess.icasig = zeros(parameters.NITER, size(signalprocess.wSIG,2));
signalprocess.SIL = zeros(1, parameters.NITER);
signalprocess.CoV = zeros(1, parameters.NITER);
idx1 = zeros(1, parameters.NITER);
% Find the index where the square of the summed whitened vectors is
% maximized and initialize W with the whitened observations at this time

% Supplement variables
signalprocess.MUFilters_nonnull = zeros(size(signalprocess.wSIG,1), parameters.NITER);

for j = 1:parameters.NITER
    %
    if j == 1 
        signalprocess.X = signalprocess.wSIG; % Initialize X (whitened signal), then X: residual

        if parameters.initialization == 0         
            actind = sum(signalprocess.X,1).^2;
            signalprocess.actind = actind;
            [~, idx1(j)] = max(actind);
            signalprocess.w = signalprocess.X(:, idx1(j)); % Initialize w
        elseif parameters.initialization == 1
            signalprocess.w = randn(size(signalprocess.X,1),1);
        end


    else
        if parameters.initialization == 0
            actind(idx1(j-1)) = 0; % remove the previous vector
            [~, idx1(j)] = max(actind);
            signalprocess.w = signalprocess.X(:, idx1(j)); % Initialize w
        elseif parameters.initialization == 1
            signalprocess.w = randn(size(signalprocess.X,1),1);  
        end
    end

    maxiter =200; % max number of iterations for the fixed point algorithm
    [signalprocess.w,k]  = adamalg(signalprocess.w, signalprocess.X, signalprocess.B , maxiter, parameters.contrastfunc);
  
    %% Step 4: Minimization of the CoV of discharge times (end when CoV is minimized)
    
    % Initialize CoV (variation of interspike intervals, %) Step 4a => 4e
    [signalprocess.icasig, signalprocess.spikes] = getspikes(signalprocess.w, signalprocess.X, signal.fsamp);
    
    if length(signalprocess.spikes) > 10
        ISI = diff(signalprocess.spikes/signal.fsamp); % Interspike interval
        signalprocess.CoV(j) = std(ISI)/mean(ISI); % Coefficient of variation
        Wini = sum(signalprocess.X(:,signalprocess.spikes),2); % update W by summing the spikes
        
        % Minimization of the CoV of discharge times (end when CoV is minimized)
        [signalprocess.MUFilters(:,j), signalprocess.spikes, signalprocess.CoV(j)] = minimizeCOVISI(Wini, signalprocess.X, signalprocess.CoV(j), signal.fsamp);
        signalprocess.B(:,j) = signalprocess.w;

        % Calculate SIL values
        [signalprocess.icasig, signalprocess.spikes, signalprocess.SIL(j)] = calcSIL(signalprocess.X, signalprocess.MUFilters(:,j), signal.fsamp);

    else
        signalprocess.B(:,j) = signalprocess.w;
    end
end
    
    
% Filter out MUfilters below the SIL threshold
signalprocess.MUFilters(:,signalprocess.SIL < parameters.silthr) = [];
if parameters.covfilter == 1
    signalprocess.CoV(signalprocess.SIL < parameters.silthr) = [];
    signalprocess.MUFilters(:,signalprocess.CoV > parameters.covthr) = [];
end


% Batch processing over each window
[PulseT, distime,signalprocess.centroids] = batchprocessfilters(signalprocess.MUFilters, signalprocess.wSIG, signalprocess.coordinatesplateau, signalprocess.exFactor, parameters.differentialmode, size(signal.data,2), signal.fsamp);

if size(PulseT,1) > 0
    % Remove duplicates
    [PulseT, distimenew] = remduplicates(PulseT, distime, distime, round((signal.fsamp/100)), 0.00025, parameters.duplicatesthresh, signal.fsamp);
    signal.Pulsetrain = PulseT;
    for j = 1:length(distimenew)
        signal.Dischargetimes{j} = distimenew{j};
    end

end

end