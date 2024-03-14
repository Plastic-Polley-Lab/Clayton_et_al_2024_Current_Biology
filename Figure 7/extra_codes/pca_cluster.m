%% let's try to classify the responses
function [T, idx] = pca_cluster(data, bin, clusterN)
%INPUT:       data: a matrix where each row is a neuron, and each collumn
%                   represents a variable (e.g., sound intensity)
%             bin: just ignore this
%             clusterN: number of clusters you need,
%Example:

%    [T, idx] = pca_cluster(data, [], 6)


% resp = [data.resp];
% resp_indx = find(resp~=0);
% bin = 10;
% for i = 1:length(data)
%     psth_pop(i,:) = data(i).(['bin', num2str(bin)]).PSTH_auROC;
%     psth_pop_t(i,:)= data(i).(['bin', num2str(bin)]).time;
% end

% use pca to reduce the dimension of variable
[coeff,score,~,~,explained,mu] = pca(data);

sum_explained = 0;
idx = 0;
while sum_explained < 90 % using 90% variance explained
    idx = idx + 1;
    sum_explained = sum_explained + explained(idx);
end
idx 

% use tree to classify the responses
spike_score = score(:, 1:idx);
resp_ap= spike_score;

% Z = linkage(resp_ap,'average','euclidean');
Z = linkage(resp_ap,'ward','euclidean');
% Z = linkage(resp_ap,'average','correlation');

% clusterN = 4;
T = cluster(Z,'maxclust',clusterN); % edit by ke, increase the cluster from 9 to 16.
cutoff = median([Z(end-clusterN+2,3) Z(end-clusterN+2,3)]);
figure
[H,tt,outperm] = dendrogram(Z,0,'ColorThreshold',cutoff);

figure;
set(gcf, 'Position', [200, 200, 600, 400])

for i = 1:clusterN
    subplot(1,clusterN,i)
    ind_cluster = find(T == i);
    avg_pop = data(ind_cluster,:);
    plot(mean(avg_pop,1))
    hold on
    title(['Cluster ', num2str(i)])
    box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
xticks([1,9, 17])
xticklabels({'0', '40',  '80'})

end

    
figure;
set(gcf, 'Position', [200, 200, 800, 400])

for i = 1:clusterN
    subplot(1,clusterN,i)
    ind_cluster = find(T == i);
    avg_pop = data(ind_cluster,:);
    plot(avg_pop')
    hold on
    title(['Cluster ', num2str(i)])
    box off
    set(gca,'TickDir','out')
    set(gca,'fontsize',12)
    set(gca,'TickLengt', [0.015 0.015]);
    set(gca, 'LineWidth',1)
    xticks([1,9, 17])
    xticklabels({'0', '40',  '80'})
    xlabel('Intensity(dB SPL)')

end

% eva = evalclusters(resp_ap,'kmeans','CalinskiHarabasz','KList',[1:6])
% 
% opts = statset('Display','final');
% [idx,C] = kmeans(resp_ap,clusterN,'Distance','cityblock',...
%     'Replicates',5,'Options',opts);
% 
% figure;
% for i = 1:clusterN
%     subplot(1,clusterN,i)
%     ind_cluster = find(idx == i);
%     avg_pop = data(ind_cluster,:);
%     plot(mean(avg_pop,1))
%     hold on
% end
