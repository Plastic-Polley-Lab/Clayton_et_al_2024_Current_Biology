function result = cross_coeff(scmatrix, fanoparam, fig)
% by default there are 20 chords, and each chord last 50 ms
% let's get the index of each chord
% let's store the mean spike count and variane for each chord

timeOnset = fanoparam.timeOnset;
timeWindow = fanoparam.timeWindow;

psth    = scmatrix(:, timeOnset + 1 : end);
% w = gausswin(timeWindow);
% w = w/sum(w);
% psth_filter    = filter(w, 1, psth, [], 2); % introduce phase shift
psth_filter = smoothts(psth, 'g', timeWindow, 4);
% idx = [];
% for i = 1:size(psth_filter,1)
%     if sum(psth_filter(i,:)) ==0
%         idx = [idx,i];  % trials without spikes
%     end
% end
% psth_filter(idx,:) = []; % remove trials without spikes at all
R = corrcoef(psth_filter');
R_cross = R -  diag(diag(R));
R_avg   = sum(R_cross(:))/(size(R,1) *(size(R,1)-1));
if isnan(R_avg)
    idx = isnan(R_cross);
    R_cross(idx) = 0;
    R_avg   = sum(R_cross(:))/(size(R,1) *(size(R,1)-1));
end

result.corrR = R;
result.corrR_avg = R_avg;


