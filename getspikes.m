%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To identify the discharge times of the motor unit

% Input: 
%   w = weigths
%   X = whitened signal
%   fsamp = sampling frequency

% Output:
%   icasig = MU pulse train
%   spikes2 = discharge times of the motor unit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [icasig, spikes2] = getspikes(w, X, fsamp)

icasig = (w' * X).*abs(w' * X);
[~,spikes] = findpeaks(icasig, 'MinPeakDistance', round(fsamp*0.01));
icasig = icasig/mean(maxk(icasig(spikes),10));
if length(spikes)>1
    [L,C] = kmeans(icasig(spikes)',2);
    [~, idx2] = max(C);
    spikes2 = spikes(L==idx2);
    spikes2(icasig(spikes2)>mean(icasig(spikes2))+3*std(icasig(spikes2))) = [];
else
    spikes2 = spikes;
end
