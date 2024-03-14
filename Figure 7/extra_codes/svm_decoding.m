function [accuracy_decoding,Confusion_all]= svm_decoding(data,n_kfold, n_repeats)
% INPUT:
%       data: a table structure; inclues different features; the target is
%       named as classes
%       for instance; load fisheriris; data = table(meas); data.classes =
%       species.
C = [0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30];
repeats = n_repeats;
svm_model =[];
CV_svmModel =[];
accuracy =[];
try
    temp = fitcsvm(data,'classes','Standardize', true, 'KernelScale', 'auto', 'BoxConstraint', C(1));
    idx = 1
catch
    idx = 2
    display ('More than two classes for prediction')
end
    
for j = 1: repeats
    for i = 1:length(C)
        switch idx
            case 1
             svm_model{i,j} = fitcsvm(data,'classes','Standardize', true, 'KernelScale', 'auto', 'BoxConstraint', C(i)); % data can be a table
            case 2
             t = templateSVM('Standardize',1,'KernelScale', 'auto','BoxConstraint',C(i));
             svm_model{i,j} = fitcecoc(data,'classes', 'Learners', t); % data can be a table
        end
             CV_svmModel{i,j} = crossval(svm_model{i,j},'kfold', n_kfold);
        % Compute validation accuracy
        accuracy(i,j) = 1 - kfoldLoss(CV_svmModel{i,j}, 'LossFun', 'ClassifError');
        fprintf('Repeats %d, Training with BoxConstraint value %d, accuracy is %0.4f\n',j, C(i), accuracy(i))
    end
end
accuracy_avg = mean(accuracy,2);
[M,I] = max(accuracy_avg);
C_optimal = C(I)
accuracy_decoding = accuracy(I,:);
predictions =[];
validationScores =[];
order = [];
for i = 1: repeats
    [predictions{i}, validationScores{i}] = kfoldPredict(CV_svmModel{I,i});
    [Confusion_matrix(:,:,:,i),order{i}] = confusionmat(data.classes, predictions{i});
end
Confusion_all = sum(Confusion_matrix,4);
true_tot = sum(Confusion_all,2);
true_tot = repmat(true_tot, 1, size(Confusion_all,2));
figure; 
imagesc(Confusion_all./true_tot)
xticks([1:length(order{i})])
xticklabels(order{1})
yticks([1:length(order{i})])
yticklabels(order{1})
ylabel('True Classes')
xlabel('Predicted Classes')
colormap(jet)
box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,400,300])
hold off