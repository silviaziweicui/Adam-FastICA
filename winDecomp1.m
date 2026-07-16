%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% winDecomp1
% 
% This function decomposes EMG signals within a sliding window for real-time 
% FastICA. It performs online update of separation parameters and extracts 
% pulse trains for each motor unit.
%
% Input: 
%   winData: current window EMG data (channels x samples)
%   signalprocess: pre-trained signal processing structure
%   signalprocess_win: window-level signal processing structure (updated)
%   med: update mode (1 = single window, 2 = buffer-based update)
%   data_buffer: historical data buffer (used when med = 2)
%
% Output:
%   IPTs: independent potential traces for all MUs (MUnum x winLen)
%   Pulses: spike trains for all MUs (cell array)
%   signalprocess_win: updated window-level signal processing structure
%       - signalprocess_win.Ctx: updated Ctx matrix
%       - signalprocess_win.CorrSig_all: updated accumulated correlation matrix
%       - signalprocess_win.centroids: updated cluster centroids
%       - signalprocess_win.eSigLen: updated extended signal length
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [IPTs, Pulses, signalprocess_win] = winDecomp1(winData, signalprocess,  signalprocess_win, data_buffer) 


exFactor = signalprocess.exFactor;

eSigLen = signalprocess_win.eSigLen;
CorrSig_all = signalprocess_win.CorrSig_all;
Ctx = signalprocess_win.Ctx;
centroids = signalprocess_win.centroids;

winLen = size(winData,2);
tmpeSig = extend(winData,exFactor);

if nnz(tmpeSig) == 0
    W = signalprocess_win.MUFilters;
    munum = size(W,2);
    IPTs = zeros(munum,winLen);
    Pulses = cell(1,munum);
    return
end

Yq = tmpeSig(:,exFactor+1:end-exFactor); 

LR = 0.1;
if ~isempty(data_buffer)
    buffer_tmpeSig = extend(data_buffer, exFactor);  
    buffer_Yq = buffer_tmpeSig(:, exFactor+1:end-exFactor);

    buffer_Yq = sparse(buffer_Yq);
    CorrSig_all = CorrSig_all + buffer_Yq * buffer_Yq';
    eSigLen = eSigLen + size(data_buffer,2) + exFactor - 1;

else
    Yq = sparse(Yq);
    CorrSig_all = CorrSig_all + Yq * Yq';
    eSigLen = eSigLen + winLen + exFactor - 1;
end

CorrSig = CorrSig_all / eSigLen;
invCorrSig = pinv(CorrSig);
W = (Ctx' * invCorrSig)';  


munum = size(W,2);
Pulses = cell(1,munum);
IPTs = zeros(munum,winLen);
pNum = zeros(1,munum);

%% 
tmpIPTs = W' * tmpeSig; 

for mu = 1:munum
    tmpT = tmpIPTs(mu,:);
    tmpT([1:exFactor,end-exFactor+1:end]) = 0;
    tT = abs(tmpT).*tmpT;
    if  -min(tmpT) > max(tmpT)
        tT(find(tmpT>0)) = 0;     
    else
        tT(find(tmpT<0)) = 0;
    end
    tT = abs(tT);

    C1 = centroids(mu,1); 
    C2 = centroids(mu,2);

    % save IPT
    IPTs(mu,:) = tT(1:end-exFactor+1);

    aveC = (C1+C2)/2;
    tmpInd = find(tT(exFactor+1:end-exFactor)>aveC);         

    intervalLimit = 40;
    compInd = remRepeatedInd(tT(exFactor+1:end-exFactor),tmpInd,intervalLimit);
   
    compInd = compInd+exFactor; 
    Pulses{mu} = compInd;
    pNum(mu) = length(compInd);
    
    if ~isempty(compInd)          
        clen = length(compInd);

        tmpCtx = sum(tmpeSig(:,compInd),2);
        tmpCtx = tmpCtx/sqrt(sum(tmpCtx.^2));
        Ctx(:,mu) = Ctx(:,mu)+LR*tmpCtx;
        Ctx(:,mu) = Ctx(:,mu)/sqrt(sum(Ctx(:,mu).^2));

        w_C1 = 5;
        tmpC1 = (C1*w_C1+sum(tT(compInd)))/(w_C1+clen);            
        Ind2 = setdiff(1:length(tT),compInd);  
        len2 = length(Ind2);

        tmpC2 = (C2*w_C1+sum(tT(Ind2)))/(w_C1+len2);

        centroids(mu,1) = tmpC1;
        centroids(mu,2) = tmpC2;            
    end                      
end

if sum(pNum)<munum
    CorrSig_all = CorrSig_all-Yq*Yq';    
end
signalprocess_win.Ctx = Ctx;
signalprocess_win.CorrSig_all = CorrSig_all;
signalprocess_win.centroids = centroids;
signalprocess_win.eSigLen = eSigLen;

