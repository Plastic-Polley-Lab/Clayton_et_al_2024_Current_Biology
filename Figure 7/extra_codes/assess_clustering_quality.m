function [clusterIDs,unitQualities,contaminationRates, isiViolations] = assess_clustering_quality(penDir)

%unitQualities is a vector of isolation distances for the clusterIDs

[clusterIDs, unitQualities, contaminationRates, isiViolations] = sqKilosort.computeAllMeasures(penDir);

% %% plot the estimated false positive rate
% [counts,centers]=hist(isiViolations);
% h= figure;
% bar(centers, counts);
% title('Refractory Period Violations');
% xlabel('False Positive Rate');
% ylabel('Number of Clusters');


end