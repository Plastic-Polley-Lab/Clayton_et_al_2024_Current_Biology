function kilosort_spiketimes_Tosca_clusterData(analdir)

% kilosort_spiketimes_to_scl.m

% Function to load in all of the spike times from a complete kilosort run,
% and save out all the relevant information (times, ticks etc) for each
% identified unit.
%
% Usage:
%       kilosort_spiketimes_to_scl(analdir)
%
% Options:
%       analdir:    analysis directory containing the kilosort output
%
% (c) Ross S Williamson, June 2016; ross.s.williamson@gmail.com;

%%Adapted by M.M.A, June 2018; 

%%%%%%
% %It needs to be forward slash since that's the delimiter used later
% 
% clearvars;
% analdir = 'F:/EPHYS/Passive/MMA012420/Sorted_Data/ACtx/Penetration_2'
% % chanmap_file = 'F:/EPHYS/Meenakshi_Code/Channel_Maps/chanMap64.mat'
% chanmap_file = 'F:/EPHYS/Meenakshi_Code/Channel_Maps/chanMap32x2.mat'



% p = setPaths();

Fs = 24414.0625; %sampling rate

%loads in the channel map (if needed)
%default: corresponds to 32 channel HZ package edge probe from Neuronexus
% load(fullfile(p,'pipeline-resources',chanmap_file));
% load(chanmap_file);

if exist(fullfile(analdir,  'spike_clusters.npy')) ~=0
    %phy has been opened for this file already
    clu = readNPY(fullfile(analdir,  'spike_clusters.npy'));
else
    %because spike_templates is the same as spike_clusters prior to manual
    %curation
    clu = readNPY(fullfile(analdir,  'spike_templates.npy'));
end

ss = readNPY(fullfile(analdir,  'spike_times.npy'));
st = double(ss)/Fs;

if  exist(fullfile(analdir,  'cluster_group.tsv')) ~=0
    %manual curation has already been carried out
    [cids, cgs] = readClusterGroupsCSV(fullfile(analdir,  'cluster_group.tsv'));
else
    %manual curation has not been carried out, so unique(clu) represents
    %the direct output of kilosort
    cids = unique(clu); %all clusters
    
    %label these as -1: unsorted
    cgs = zeros(1,length(cids)).*-1; 
end



[clusterID,isoDistances,contaminationRates,isiV] = assess_clustering_quality(analdir);
split_strings = strsplit(analdir,'\');

% %%%changed penDir (in Ross' code) to analdir
% temps = readNPY(fullfile(analdir, 'templates.npy'));
% spikeTemplates = readNPY(fullfile(analdir,  'spike_templates.npy')); % note: zero-indexed

%load in the timepoints (corresponding to breaks between parameter files)
cd(analdir)
% disp('loading timepoints...');
% % load(fullfile(analdir,  sprintf('%s-%d_complete_timepoints.mat',split_strings{5},str2double(split_strings{8}(end)))));
% % %load in the ticks
% % disp('loading ticks...');
% % load(fullfile(analdir,  sprintf('%s-%d_complete_ticks.mat',split_strings{5},str2double(split_strings{8}(end)))));
% % %load in the files
% % disp('loading files...');
% % load(fullfile(analdir,  sprintf('%s-%d_complete_files.mat',split_strings{5},str2double(split_strings{8}(end)))));
% file = [dir('*_complete_timepoints.mat'), dir('*_complete_ticks.mat'), dir('*_complete_files.mat')];
% for i = 1:length(file)
%     disp(file(i).name)
%     load(file(i).name)
% end
%st contains the spike times for all identified units
%clu contains the cluster identities
%cids tells us which clusters are good
goodClusters = cids(cgs==2);
allClusters = cids;

disp('creating SpikeData directory...');

%check for existence of SpikeData directory and wipe if so (because new
%clusters are going to be created)

spikedata_dir = [analdir '_SpikeData/'];
if exist(spikedata_dir) == 7 %7 means directory exists
    delete([spikedata_dir '*.mat']);
else
    mkdir(spikedata_dir);
end

for i=1:length(goodClusters)
    %for each goodCluster find their location in the clu vector, then the
    %st vector
    
    %good cluster i is cluster number "clus_id"
    clus_id = goodClusters(i);
    clus_rating = cgs(find(cids==clus_id));
    
    spikes = st(find(clu == clus_id));
    nSpikes = length(spikes);
%     idx = find(clusterID == clus_id+1);
    
%     %check isi is the same
%     spike_hist = histcounts(spikes.*1000,0:1:spikes(end).*1000);
%     test = xcorr(spike_hist,spike_hist,25);
%     test(26) = 0;
%     figure; bar(-25:1:25,test);

%     isoDist = isoDistances(idx); %clus_id is zero-indexed
%     isiViol = (isiV(idx) ./ nSpikes) .*100; %as a percentage
    
%     [templates, templateMinIndex] = get_mean_spatiotemporal_template(analdir, clus_id);
    [templates, templateMinIndex] = get_mean_spatiotemporal_template_ks2(analdir, clus_id);
    
    clusterData = [];
%     for j = 1:length(timepoints)
    j = 1;  % edit by Ke, as there is only 1 run for tosca data.
        
        %location of the SCL file corresponding to the particular run
%         clusterData.run(j).filename = save_files{j};
%         clusterData.run(j).chanmap_file = chanmap_file; %details of the electrode used for this recording

        %unit quality metrics
%         clusterData.run(j).isolationDistance = isoDist;
%         clusterData.run(j).numISIviols = isiV(idx);
        clusterData.run(j).nSpikes = nSpikes;
%         clusterData.run(j).isi_perc = isiViol;


        clusterData.run(j).clusterID = clus_id;
        clusterData.run(j).clusterRating = clus_rating;
        
        %if it is a merged cluster, this represents the average
        clusterData.run(j).spatiotemporalTemplate = templates;
        %channel on which the minimum template was found
        clusterData.run(j).templateMin = templateMinIndex;
        
        
        clusterData.run(j).spiketimes = spikes *1000;


     load('tosca_event.mat');      
    clusterData.run(j).stimData = Tosca; 
    disp(sprintf('...saving spiking data for cluster %d/%d',i,length(goodClusters)));    
    filename = [analdir '_SpikeData/'  split_strings{end-2} '-' split_strings{end-1} '-' sprintf('cluster_%d',clus_id) '.mat'];
    save(filename,'clusterData');    
end
end



