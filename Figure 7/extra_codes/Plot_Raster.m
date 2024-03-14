file = dir('*.mat');
for i = 1:10
    load(file(i).name)
    Generate_Raster(clusterData,2)
    xlabel('Time (ms)')
    ylabel('Avg Spike Count')
    Generate_Raster(clusterData,2)
    xlabel('Time (ms)')
    ylabel('Avg Spike Count')
end