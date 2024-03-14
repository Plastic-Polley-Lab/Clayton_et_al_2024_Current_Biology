%% load sorted data
function spikes = loadSingleUnits(path)
if nargin < 1
else
    cd path
end
spike_times = readNPY('spike_times.npy');
spike_clusters = readNPY('spike_clusters.npy');
[~, ~, cluster_group] = tsvread('cluster_group.tsv');

k = 1;
for i = 1:size(cluster_group)
    if strcmp(cluster_group{i,2},'good')
        spike_good(k) = str2num(cluster_group{i,1}); %#ok<*ST2NM>
        k = k+1;
    end
end

for i = 1:length(spike_good)
    spikes(i).timestamps = spike_times(find(spike_clusters==spike_good(i)));
end

end