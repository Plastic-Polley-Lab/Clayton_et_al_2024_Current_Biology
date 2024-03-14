function [templates, templateMinIndex] = get_mean_spatiotemporal_template_ks2(penDir, clusterID)
    % returns the n-channel spatiotemporal template that is the mean of
    % all templates used for a given cluster.  Also returns which channel
    % of that template has the maximally minimum value.
    
    %Just to clarify the indexing here:
    %cluster id's are zero indexed, but a cluster id of 0 corresponds to a
    %template vector location of 1 (obviously). 
    %This code just adds 1 to cluster ID and grabs those templates from the
    %rez.Wraw matrix
    
     if exist(fullfile(penDir,  'spike_clusters.npy'), 'file') ~=0
        %phy has been opened for this file already
        clu = readNPY(fullfile(penDir,  'spike_clusters.npy'));
    else
        %because spike_templates is the same as spike_clusters prior to manual
        %curation
        clu = readNPY(fullfile(penDir,  'spike_templates.npy'));
        %clu = clu+1; %because this file is zero indexed
     end   
    
    spikeTemplates = readNPY(fullfile(penDir,  'spike_templates.npy')); % note: zero-indexed
%     temps = readNPY(fullfile(penDir,  'templates.npy'));
    
    templatesUsed = unique(spikeTemplates(clu==clusterID))+1;
    raw_templates = readNPY(fullfile(penDir,  'templates.npy'));
%     load(fullfile(penDir,'rez2.mat')); %the raw templates are in rez.Wraw
%     raw_templates = rez.Wraw;
%     clear rez;
    
 
%     templates = zeros(length(templatesUsed),size(raw_templates,1),size(raw_templates,2));
    templates = zeros(length(templatesUsed),size(raw_templates,2),size(raw_templates,3));
    for j = 1:length(templatesUsed)
        %get the 32-channel template with the given template index
        %idx = find(unique(spikeTemplates)==templatesUsed(j));
       idx = templatesUsed(j);
%        templates(j,:,:) = raw_templates(:,:,idx);
       templates(j,:,:) = raw_templates(idx,:,:);
    end
    templates = squeeze(mean(templates,1));

    %Which channel had the maximally minimum spike
    [~,templateMinIndex] = min(min(templates,[],1));
end


%compare temps to raw_templates



% figure;
% imagesc(squeeze(temps(7,:,:))');
% figure;
% imagesc(raw_templates(:,:,7));