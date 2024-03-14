%Plot Figure 3E 
clear all
close all

%CHANGE THIS PATH 
%Load decoder results used in manuscript
load('C:\Users\kx776\Dropbox\codeBase\Videography\figure_making\summary_data\figure 2\SVM_data_speech_masker.mat')

mSVM_accuracy = squeeze(mean(svmAccuracy,3))*100
mSVM_accuracy_Shuffled = squeeze(mean(svmAccuracy_shuffled,3))*100

%% Now generate line plots
figure;
add_ind = 0:2:12; 
for i = 1:6
   
   plot([1 2]+add_ind(i), [mSVM_accuracy(:,i) mSVM_accuracy_Shuffled(:,i)],'Color',[0.5 0.5 0.5]); 
   hold on
   plot([0.75 1]+add_ind(i), ones(1,2)*mean(mSVM_accuracy(:,i)),'k','LineWidth',2)
    plot([2.0 2.25]+add_ind(i), ones(1,2)*mean(mSVM_accuracy_Shuffled(:,i)),'k','LineWidth',2)
end
xlim([0 13]) 
box off 
hold on
yline(50,'--') 
ylabel('Classification accuracy (%)') 
xticks([1.5:2:12]) 
xticklabels({'-Inf','10','20','30','40','50'}) 
xlabel('Masker level (dB SPL)') 
set(gcf,'Position',[178         504        1161         420])

%% Paired ttests on shuffled vs. non-shuffled data 
for i = 1:6
    [h(i) p(i)] = ttest(mSVM_accuracy(:,i),mSVM_accuracy_Shuffled(:,i))
end
