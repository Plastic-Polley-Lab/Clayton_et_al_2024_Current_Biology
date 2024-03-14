clear all
close all

%CHANGE THIS PATH 
files = {'C:\Users\kx776\Dropbox\codeBase\Videography\summary_data\speech_in_noise\ChAT_151_speech_in_noise.mat',...
    'C:\Users\kx776\Dropbox\codeBase\Videography\summary_data\speech_in_noise\ChAT_152_speech_in_noise.mat',...
    'C:\Users\kx776\Dropbox\codeBase\Videography\summary_data\speech_in_noise\ChAT_153_speech_in_noise.mat',...
    'C:\Users\kx776\Dropbox\codeBase\Videography\summary_data\speech_in_noise\ChAT_154_speech_in_noise.mat',...
    'C:\Users\kx776\Dropbox\codeBase\Videography\summary_data\speech_in_noise\ChAT_155_speech_in_noise.mat',...
    'C:\Users\kx776\Dropbox\codeBase\Videography\summary_data\speech_in_noise\ChAT_156_speech_in_noise.mat',...
    'C:\Users\kx776\Dropbox\codeBase\Videography\summary_data\speech_in_noise\ChAT_157_speech_in_noise.mat',...
    'C:\Users\kx776\Dropbox\codeBase\Videography\summary_data\speech_in_noise\ChAT_150_speech_in_noise.mat'}

files = unique(files); %ensure no duplicates

%% Load in data
sprate_mat = [ ];
for m = [1:8]% length(files)
    clear dat
    load(files{m})
    sprate_mat(m,:,:,:,:) = dat.sorted_face_stim;
end

%% Run Classifier
sprate_mat;
for m = 1:8
    for snr = 1:6
        clear temp_mat temp_mat_ds Mdl
        figure
        x = [squeeze(sprate_mat(m,1,snr,:,800:1699));squeeze(sprate_mat(m,2,snr,:,800:1699))];
        labels_ind = randsample(1:40,40,'false');
        labels = [zeros(20,1);ones(20,1)];
        labels_shuffle = labels(labels_ind);
        for k = 1:1000
            [xx svmAccuracy(m,snr,k)] = trainClassifier_speech_tokens_v3(x,labels);
            [xx svmAccuracy_shuffled(m,snr,k)] = trainClassifier_speech_tokens_v3(x,labels_shuffle);
        end
    end
end

%Plot out average classification 
% figure;
% box off
% scatter(ones(1,8)*1+.25*rand(1,8)-0.125,mean(svmAccuracy,2) *100,400,'r.')
% hold  on
% plot([1-0.125 1+0.125],ones(1,2)*mean(mean(svmAccuracy,2))*100,'k','LineWidth',2)
% scatter(ones(1,8)*2+.25*rand(1,8)-0.125,mean(svmAccuracy_shuffled,2)*100,400,'r.')
% plot([2-0.125 2+0.125],ones(1,2)*mean(mean(svmAccuracy_shuffled,2))*100,'k','LineWidth',2)
% xlim([-.5 3.5])
% yline([50],'--') 
% xticks([1 2])
% ylabel('Classification accuracy (%)')
% xticks([1 2]) 
% xticklabels({'True trial labels','Shuffled trial labels'})
% ylim([30 100]) 
% savedir = 'C:\Users\kx776\Dropbox\Documents\Presentations\MEE research day 2023\schematics'; 
% print(gcf,[savedir 'Classifier_performance.pdf'],'-dpdf')

% %Save model fits 
% savedir = 'C:\Users\kx776\Dropbox\codeBase\Videography\figure_making\summary_data\figure 2\'; 
% save([savedir 'SVM_data.mat'],'svmAccuracy','svmAccuracy_shuffled')


 

