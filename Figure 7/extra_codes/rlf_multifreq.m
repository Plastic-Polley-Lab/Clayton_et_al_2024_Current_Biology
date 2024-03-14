%% here is a new code using not only the best frequency, but the nearby frequency
function [rlf,rlf_multifreqs, spont, fra, spls, freqs] = rlf_multifreq(data)

%% Step1: find the best frequency;

% i = 1;
analysis_window = [1:100]; % the windows used to get the fra
% data = temp(i);
tuning = data.tuning; % get the tuning curve
freqs = data.freqs;   % frequencies used
tuning_avg = mean(tuning,2); % average tuning curve
[~, best_indx] = max(tuning_avg); % get the index of the best frequency
best_freq = freqs(best_indx); % this is the best frequency

%% Step2: find k frequency away from the best frequency
k = 1; %  k different frequencies away from the best frequency
up_indx =[]; % frequency bigger than best frequency
low_indx =[]; % frequency smaller than best frequency


if best_indx == 1 % best frequency is the lowest frequency
    low_indx = [];
    up_indx = best_indx + (1:k);
end

if best_indx == length(freqs) % best frequency is the highest frequency
    low_indx = best_indx - (1:k);
    up_indx = [];
end

if best_indx > 1 && best_indx < length(freqs) % best frequency is in the mid range
    low_indx = best_indx - (1:k);
    up_indx  = best_indx + (1:k);
    
    if min(low_indx)<1 % in case where the lower frequencies could be outside the range
        low_indx = 1:max(low_indx);
    end
    
    if max(up_indx)>length(freqs) % in case where the higher frequencies could be outside the range
        up_indx = min(up_indx) : length(freqs);
    end
end
%% Step 3: use the frequency range low: high to get the rate level function
range = sort([low_indx, best_indx, up_indx]);
fra = data.fra;
rlf_multifreqs = squeeze(mean(fra(range, :, :), 1));
windows = analysis_window(end) -analysis_window(1) + 1;
rlf_multifreqs = mean(rlf_multifreqs/windows * 1000, 2);
rlf            = mean(data.rlf,2);
spont          = data.spont;
fra            = data.fra;
spls           = data.spls;
freqs          = data.freqs;
% figure;
% plot(mean(rlf_multifreqs,2))
% hold on
% plot(mean(data.rlf,2))


