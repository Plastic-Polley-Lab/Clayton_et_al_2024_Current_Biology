% only choose the neurons activated for correlation analysis
% remove the units with isi violation bigger than 0.5%
function sig_noise_process(fra_summary, sign)
% INPUT:      fra_summary: preprocess data 
%             sign: 1, represent neurons exciated, -1 represent neurons
%             inhibited

fra_summary([fra_summary.refra_violation_ratio]>0.005) =[];
switch sign
    case 1
        label      = 'exc_indx';
        data_label = 'exc_data';
    case -1
        label      = 'inh_indx';
        data_label = 'inh_data';
end
sig_cor.(label) = find([fra_summary.resp]==sign);
n_exc    = length(sig_cor.(label));
% let's first calculate the signal correlation for 70 dB

sig_cor.(data_label) = fra_summary(sig_cor.(label));
temp_data = sig_cor.(data_label);
iter = 1;
sig_cor.cor70 = zeros( n_exc * (n_exc-1)/2, 1);
sig_cor.cor50 = zeros( n_exc * (n_exc-1)/2, 1);
sig_cor.corall = zeros( n_exc * (n_exc-1)/2, 1);


for i = 1:(length(temp_data)-1)  % calculate the pairwise correlation
    
    for j = (i+1): length(temp_data)
        fra_neuron1 = temp_data(i).fra;
        fra_neuron2 = temp_data(j).fra;
        
        fra_neuron1_tuning = mean(fra_neuron1,3); % the 1st column is 50 dB, 2nd column is 70 dB
        fra_neuron2_tuning = mean(fra_neuron2,3); % the 1st column is 50 dB, 2nd column is 70 dB
        % get the best freq for each neurons
        [~, sig_cor.cor70Best(iter,1)] = max(fra_neuron1_tuning(:,2));
        [~, sig_cor.cor70Best(iter,2)] = max(fra_neuron2_tuning(:,2));
        
        % add the location of the shank
        sig_cor.shank(iter,1) = temp_data(i).shank_x;
        sig_cor.shank(iter,2) = temp_data(j).shank_x;
        
        
        
        temp = corrcoef(fra_neuron1_tuning(:,2), fra_neuron2_tuning(:,2)); % 70 dB
        sig_cor.cor70(iter) = temp(1,2); 
        
        temp = corrcoef(fra_neuron1_tuning(:,1), fra_neuron2_tuning(:,1)); % 50 dB       
        sig_cor.cor50(iter) = temp(1,2); 

        temp = corrcoef(mean(fra_neuron1_tuning,2), mean(fra_neuron2_tuning,2)); % average tuning based on 50 and 70 dB
        sig_cor.corall(iter) = temp(1,2);
        iter = iter + 1; 
    end
    
end

%% check the signal correlation
figure;
h1 = scatter(sig_cor.cor70, sig_cor.cor50)
hold on
h2 = scatter(sig_cor.cor70, sig_cor.corall)
legend([h1, h2], {'70 vs 50 dB', '70 vs All'})
xlabel('Signal Correlation Coefficient at 70 dB')
ylabel('Signal Correlation Coefficent at 50 dB or All')
%% Let's calculate the noise correlation
% first noise correlation for each stimuli at 70 dB
iter = 1;
n_freq = 9; % number of center frequencies used
noise_cor.cor70all = zeros( n_exc * (n_exc-1)/2, 1);   % collapse all frequencies at 70 dB
noise_cor.cor70    = zeros( n_exc * (n_exc-1)/2, n_freq); % noise correlation for each frequence at 70 dB
noise_cor.cor50all = zeros( n_exc * (n_exc-1)/2, 1);   % collaps all frequencies at 50 dB
noise_cor.cor50    = zeros( n_exc * (n_exc-1)/2, n_freq); % noise correlation for each frequence at 50 dB
noise_cor.cor_all  = zeros(n_exc * (n_exc-1)/2, 1);   % collaps all frequencies at 50 and 70 dB


for i = 1:(length(temp_data)-1)  % calculate the pairwise correlation
    
    for j = (i+1): length(temp_data)
        fra_neuron1 = temp_data(i).fra;
        fra_neuron2 = temp_data(j).fra;
        
        % let's see 70 dB
        k = 2;
        fra_neuron1_70tuning = squeeze(fra_neuron1(:,k, :)); % the 1st column is 50 dB, 2nd column is 70 dB
        fra_neuron2_70tuning = squeeze(fra_neuron2(:,k, :)); % the 1st column is 50 dB, 2nd column is 70 dB
%         temp = corrcoef(fra_neuron1_70tuning(:),fra_neuron2_70tuning(:)); % here collapse all frequcies at 70 dB, it could also be z-scored based on each frequency
%         noise_cor.cor70all(iter) = temp(1,2);   

        fra_neuron1_tuning70_zscore = zscore(fra_neuron1_70tuning, 0, 2);
        fra_neuron2_tuning70_zscore = zscore(fra_neuron2_70tuning, 0, 2);
        temp1 = fra_neuron1_tuning70_zscore(:);
        temp2 = fra_neuron2_tuning70_zscore(:);
        temp = corrcoef(temp1, temp2); % here collapse all frequcies at 70 dB, it could also be z-scored based on each frequency
        noise_cor.cor70all(iter) = temp(1,2);
        
        % let's see 50 dB
        k = 1;
        fra_neuron1_50tuning = squeeze(fra_neuron1(:,k, :)); % the 1st column is 50 dB, 2nd column is 70 dB
        fra_neuron2_50tuning = squeeze(fra_neuron2(:,k, :)); % the 1st column is 50 dB, 2nd column is 70 dB   
%         temp = corrcoef(fra_neuron1_50tuning(:),fra_neuron2_50tuning(:));
%         noise_cor.cor50all(iter) = temp(1,2);
        fra_neuron1_tuning50_zscore = zscore(fra_neuron1_50tuning, 0, 2);
        fra_neuron2_tuning50_zscore = zscore(fra_neuron2_50tuning, 0, 2);
        temp3 = fra_neuron1_tuning50_zscore(:);
        temp4 = fra_neuron2_tuning50_zscore(:);
        temp = corrcoef(temp3, temp4); % here collapse all frequcies at 70 dB, it could also be z-scored based on each frequency
        noise_cor.cor50all(iter) = temp(1,2);

        % let's collapse across all
%         temp = corrcoef([fra_neuron1_70tuning(:); fra_neuron1_50tuning(:)],[fra_neuron2_70tuning(:); fra_neuron2_50tuning(:)]);
%         noise_cor.cor_all(iter) = temp(1,2);
        temp = corrcoef([temp1; temp3],[temp2;temp4]);
        noise_cor.cor_all(iter) = temp(1,2);

        % see noise correlation for all stimuli
        R_noise =zeros(size(fra_neuron1,1), size(fra_neuron1,2));
        for nn = 1: size(fra_neuron1,1)
            for mm = 1:size(fra_neuron1,2)
                temp = corrcoef(fra_neuron1(nn,mm,:), fra_neuron2(nn,mm,:));
                R_noise(nn,mm) = temp(1,2);
            end
        end

        noise_cor.cor70(iter,:) = R_noise(:,2)'; % the 1st column is 50 dB, 2nd column is 70 dB
        noise_cor.cor50(iter,:) = R_noise(:,1)'; % the 1st column is 50 dB, 2nd column is 70 dB
       iter = iter + 1; 
    end
    
end

%% check the noise correlation
figure;
scatter(noise_cor.cor70all, noise_cor.cor50all)
hold on
scatter(noise_cor.cor70all, noise_cor.cor_all)
legend([h1, h2], {'70 vs 50 dB', '70 vs All'})
xlabel('Collapse Noise Correlation Coefficient at 70 dB')
ylabel('Collapse Noise Correlation Coefficent at 50 dB or All')
%% check sig-noise correlation
figure;
scatter(sig_cor.cor70, noise_cor.cor70all)
xlabel('Signal Correlation Coefficient at 70 dB')
ylabel('Noise correlation Coefficient at 70 dB')
%% check the sig-noise correlation for each stimuli; it is very messy, unable to pull out structures
% X = sig_cor.cor70;
% Y = 1:n_freq;
% Z = (noise_cor.cor70)';
% figure;
% mesh(X, Y, Z)
% rotate3d on
% xlabel('Signal correlation')
% ylabel('Frequencies')
% zlabel('Noise correlation')

% figure;
% hold on
% for i = 1:length(X)
%    zz = repmat(X(i), n_freq, 1);
%     yy = Z(:, i);
%     xx = 1:n_freq;
%     if X(i)< 0
%         
%         plot3(xx, yy, zz, '-k')
%     elseif X(i)>0
%         plot3(xx, yy, zz, '-r')
%     end
%     
% end
% xlabel('frequence')
% ylabel('noise correlation')
% zlabel('signal correlation')
% rotate3d on
%% check the average; still messy
% sig_positive = find(sig_cor.cor70>0);
% sig_negative = find(sig_cor.cor70<0);
% figure; plot(mean(noise_cor.cor70(sig_positive,:),1), '-k')
% hold on; plot(mean(noise_cor.cor70(sig_negative,:),1), '-r')

% figure
% x = 1:n_freq;
% y = noise_cor.cor70(sig_positive,:);
% h1 = line_sem_plot(x, y, 'k');
% hold on
% h2 = line_sem_plot(x, noise_cor.cor70(sig_negative,:), 'r')
% legend([h1, h2], {'Sig Corr > 0', 'Sig Corr < 0'})
% xlabel('Frequency')
% ylabel('Noise Correlation Coefficient')
% lineplot_error(x,noise_cor.cor70(sig_negative,:), 'g')
%% check noise correlation between neurons sharing the same best frequencies
sameBest = find(sig_cor.cor70Best(:,1) == sig_cor.cor70Best(:,2));
sameBestFreq = sig_cor.cor70Best(sameBest,1);

sameBest_sig_corr = sig_cor.cor70(sameBest);

sameBest_noise_corr = [];
for i = 1:length(sameBest)
    sameBest_noise_corr(i) = noise_cor.cor70(sameBest(i), sameBestFreq(i));
end
figure;
scatter(sameBest_sig_corr, sameBest_noise_corr,'o')
xlabel('Signal Correlation Coefficient at 70 dB')
ylabel('Noise Correlation Coefficient at 70 dB')
title('Neuron pairs with the same BF')
sig_cor.sameBest = sameBest;
sig_cor.sameBestFreq = sameBestFreq;
sig_cor.sameBest_sig_corr = sameBest_sig_corr;
noise_cor.sameBest_noise_corr = sameBest_noise_corr;
% save(['sumamry_sig_noise_', data_label, '_v2.mat'], 'sig_cor', 'noise_cor')
save(['sumamry_sig_noise_', data_label, '_v3.mat'], 'sig_cor', 'noise_cor') % add the shank location of each units

