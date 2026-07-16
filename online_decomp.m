%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% online FastICA 
% Input: 
% onlinedata: data for online decomposition
% win_size_ms: window size in milliseconds
% signalprocess: pre-trained signal processing structure
% FS: sampling frequency

% Output:
% PulseT_group: IPT from online process
% Distime_group: spike trains from online process

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PulseT_group, Distime_group] = online_decomp(onlinedata, win_size_ms, signalprocess,FS)

group_idx = 1;
PulseT_group = {[]};
Distime_group = {cell(1,size(signalprocess.MUFilters,2))};

signalprocess.Ctx = pinv(signalprocess.whiteningMatrix) * signalprocess.MUFilters;

win_size = round(win_size_ms / 1000 * FS);

buffer_size = 10;
data_buffer = [];
count = 0;
int_second = 3;
corr_interval = 1000 / win_size_ms * int_second;  

signalprocess_win = signalprocess;

for i = 1:win_size:length(onlinedata)-win_size+1
    current_window = onlinedata(:,i:i + win_size-1);    

    if size(data_buffer, 2) < buffer_size * win_size
        data_buffer = [data_buffer, current_window];
    else
        data_buffer = [data_buffer(:, win_size+1:end), current_window];
    end

    if count == corr_interval
        [ics_i, distime, signalprocess_win] = winDecomp1(current_window, signalprocess, signalprocess_win, data_buffer);
        count = 0;
    else
        [ics_i, distime, signalprocess_win] = winDecomp1(current_window, signalprocess, signalprocess_win, []);
        count = count+1;
    end

    PulseT_group{group_idx} = [PulseT_group{group_idx}, ics_i]; 
    for j = 1:size(distime, 2)
        Distime_group{group_idx}{j} = [Distime_group{group_idx}{j}, distime{j} + i - 1];
    end


end
