%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To reapply the MU filters over each segment of decomposed data

% Input: 
%   MUFilters = matrix of MU filters
%   wSIG = whitened EMG signal
%   coordinates = onset and offset of each segment of data
%   exFactor = extension factor
%   differentialmode = differential mode on or off (1 or 0)
%   ltime = duration of the raw signal 
%   fsamp = sampling frequency


% Output:
%   PulseT = Pulse train of each MU
%   distime = discharge times of the MUs
%   Centroids = Cluster centroid of MU

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [PulseT, distime, Centroids] = batchprocessfilters(MUFilters, wSIG, coordinates, exFactor, differentialmode, ltime, fsamp)

f = waitbar(0,'Batch processing MUs');
MUn = 0;
MUn = MUn + size(MUFilters,2);
x = 1/MUn;

PulseT = zeros(MUn, ltime);
distime = cell(1,MUn);
Centroids = zeros(MUn, 2);
MUnb=1; 
for j = 1:size(MUFilters,2)
    PulseT(MUnb,coordinates(1):coordinates(2)+exFactor-1-differentialmode) = MUFilters(:,j)' * wSIG;

    t_new = PulseT(MUnb,:);
    IPT = abs(t_new).*t_new;
    if  -min(t_new) > max(t_new)
        IPT(find(t_new>0)) = 0;
    else
        IPT(find(t_new<0)) = 0;
    end
    IPT = abs(IPT);

    tT = normalize(IPT, 'range');

    PulseDistance = 20;  
    [L,rawC,~,D] = kmeans(tT',2); % Kmean ++ classification 
    Ind1 = find(L==1);
    Ind2 = find(L==2);
    if rawC(1)<rawC(2)
        spikeInd = remRepeatedInd(IPT,Ind2,PulseDistance);
    else
        spikeInd = remRepeatedInd(IPT,Ind1,PulseDistance);
    end  
    [~,spikes,~] = evaluateIPT(IPT,[],fsamp,0);
    distime{MUnb} = spikes;
    
    IPT = IPT/mean(maxk(IPT(1,spikeInd),10));
    Centroids(MUnb,1) = mean(IPT(spikeInd));
    Centroids(MUnb,2) = mean(IPT(setdiff([1:length(IPT)],spikeInd)));
    PulseT(MUnb,:) = IPT;

    MUnb = MUnb+1;
end   
PulseT = PulseT(:, 1:ltime);
close(f);

