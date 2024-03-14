%PSTH CLASSIFIER!
clearvars

mean_perf_sm_ensemble = [];
k= 1;
for reg = 1:3
    
% for reg = 1
    if (reg == 1)
        load('F:\EPHYS\Passive\MMA071219\Sorted_Data\IC\units.mat');
    elseif (reg == 2)
        load('F:\EPHYS\Passive\MMA071219\Sorted_Data\MGB\units.mat');
    else
        load('F:\EPHYS\Passive\MMA071219\Sorted_Data\ACtx\units.mat');
    end

s = 9;
t = 25;
sample_t = 25;%trial dropping

% figure();
% mat = [];
% master_raster=[];
% master_raster_4msbin=[];
% master_raster_10msbin=[];
% master_raster_20msbin=[];
% mtfcluster = units(find([units.mtfcluster]==2));
% mtfcluster = units(union(find([units.mtfcluster]==1),find([units.mtfcluster]==2)));
mtfcluster = units(union(union(find([units.mtfcluster]==1),find([units.mtfcluster]==2)),find([units.mtfcluster]==3)));
ensemble_size = [1,3,5,10,20];
n_iter = 100;
perf_sm_ensemble = [];
perf_sm_ensemble_all = [];
perf_sm_ensemble_high = [];
perf_sm_ensemble_low = [];
for enss = 1:length(ensemble_size)
    ens = ensemble_size(enss);
for iter = 1:n_iter
master_raster=[];
master_raster_4msbin=[];
master_raster_10msbin=[];
master_raster_20msbin=[];
y = datasample(1:length(mtfcluster),ens,'Replace',false);
% for clu=1:length(units)
% for clu=1:length(mtfcluster)
for i=1:length(y)
    clu = y(i);
% clu = 1;
    master_raster = horzcat(master_raster,mtfcluster(clu).raster(:,501:1500));
%     master_raster_4msbin = horzcat(master_raster_4msbin,mtfcluster(clu).raster_4msbin(:,126:375));
%     master_raster_10msbin = horzcat(master_raster_10msbin,mtfcluster(clu).raster_10msbin(:,51:150));
%     master_raster_20msbin = horzcat(master_raster_10msbin,mtfcluster(clu).raster_20msbin(:,26:75));

end
% master_raster = master_raster_4msbin;
% master_raster = master_raster_10msbin;
master_raster = master_raster_20msbin;

%     master_raster = mtfcluster(clu).raster;
%     master_raster = master_raster(:,501:1500);

%     master_raster_4msbin = mtfcluster(clu).raster_4msbin;
%     master_raster = master_raster_4msbin(:,126:375);

%     master_raster_10msbin = mtfcluster(clu).raster_10msbin;
%     master_raster = master_raster_10msbin(:,51:150);

template = zeros(s,size(master_raster,2));
for i = 1:s
    sub_master_raster = master_raster((i-1)*t+1:(i-1)*t+sample_t,:);
    template(i,:) = mean(sub_master_raster,1);
end

%testing stage
min_euc = zeros(s,sample_t);
for i = 1:s
    sub_master_raster = master_raster((i-1)*t+1:(i-1)*t+sample_t,:);
    for j = 1:size(sub_master_raster,1)
        test_vector = sub_master_raster(j,:);
        euc_dist = [];
        for temp = 1:size(template,1)
            if(temp == i)
                train_vector = ((template(temp,:)*sample_t)-test_vector)/(sample_t-1);
            else
                train_vector = template(temp,:);
            end
%         train_vector = template(temp,:)-test_vector);
        euc_dist(temp) = sqrt(sum((test_vector - train_vector).^2));
        end
        tmp = find(euc_dist == min(euc_dist));
        min_euc(i,j) = tmp(1);
    end
end
% figure();
% imagesc(min_euc)

actual_vs_classified = zeros(s,s);
for i = 1:size(min_euc,1)
    for j = 1:s
    actual_vs_classified(i,j) = length(find(min_euc(i,:)==j))/sample_t;
    end
end
actual_vs_classified = fliplr(actual_vs_classified);
% subplot(ceil(length(units)/4),4,clu);
% figure();
% subplot(3,1,3);
% colormap(hot);
% imagesc(template);
% subplot(2,1,2);
% colormap(hot);
% imagesc(actual_vs_classified);
% colorbar();
% caxis([0 1])
% subplot(2,1,1);
% L = logical(master_raster);
%     MarkerFormat.Color = 'b';
%     plotSpikeRaster(L,'PlotType','scatter','MarkerFormat',MarkerFormat);
%     yticks(1:25:225);
%     mod_rates = 2.^(0:1:8);
%     ylabeling = cell(1,length(yticks));
%     for i = 1:s
%         ylabeling{i} = num2str(mod_rates(s-i+1));
%     end
%     yticklabels(ylabeling);

% perf(clu) = mean(diag(flipud(actual_vs_classified)));
% mat(:,:,clu) = actual_vs_classified;
% clu
% end
% perf_su(reg) = mean(perf);
% perf_mean = mean(mat,3);
% figure(1);
% subplot(3,1,reg);
% colormap(hot);
% imagesc(perf_mean);
% colorbar();
% caxis([0 0.5])

% perf_complete(reg) = mean(diag(flipud(actual_vs_classified)));
% perf1(reg) = mean(diag(flipud(actual_vs_classified)));
% perf2(reg) = mean(diag(flipud(actual_vs_classified)));
% perf3(reg) = mean(diag(flipud(actual_vs_classified)));
perf_sm_ensemble_all(iter) = mean(diag(flipud(actual_vs_classified)));
all = diag(flipud(actual_vs_classified));
perf_sm_ensemble_high(iter) = mean(all(6:9));
perf_sm_ensemble_low(iter) = mean(all(1:5));

mat(:,:,iter) = actual_vs_classified;
end
mean_perf_sm_ensemble_all(enss,reg) = mean(perf_sm_ensemble_all);
mean_perf_sm_ensemble_high(enss,reg) = mean(perf_sm_ensemble_high);
mean_perf_sm_ensemble_low(enss,reg) = mean(perf_sm_ensemble_low);
figure(1);
subplot(3,length(ensemble_size),k);
colormap(hot);
imagesc(mean(mat,3));
colorbar();
k = k+1;
enss
end
end
figure();
colormap(hot);
imagesc(mean_perf_sm_ensemble_all);
colorbar();
figure();
colormap(hot);
imagesc(mean_perf_sm_ensemble_high);
colorbar();
figure();
colormap(hot);
imagesc(mean_perf_sm_ensemble_low);
colorbar();
%%
figure();
subplot(1,4,1);
colormap(hot);
load('F:\EPHYS\Passive\MMA071219\mean_perf_sm_ensemble_high_1msbin.mat');
imagesc(mean_perf_sm_ensemble_high);
colorbar();
% caxis([0 1])
subplot(1,4,2);
colormap(hot);
load('F:\EPHYS\Passive\MMA071219\mean_perf_sm_ensemble_high_4msbin.mat');
imagesc(mean_perf_sm_ensemble_high);
colorbar();
% caxis([0 1])
subplot(1,4,3);
colormap(hot);
load('F:\EPHYS\Passive\MMA071219\mean_perf_sm_ensemble_high_10msbin.mat');
imagesc(mean_perf_sm_ensemble_high);
colorbar();
% caxis([0 1])
subplot(1,4,4);
colormap(hot);
load('F:\EPHYS\Passive\MMA071219\mean_perf_sm_ensemble_high_20msbin.mat');
imagesc(mean_perf_sm_ensemble_high);
colorbar();
% caxis([0 1])
%%
temp = load('F:\EPHYS\Passive\MMA071219\mean_perf_sm_ensemble_high_1msbin.mat');
mean_perf_sm_ensemble_high_1msbin = temp.mean_perf_sm_ensemble_high;
temp = load('F:\EPHYS\Passive\MMA071219\mean_perf_sm_ensemble_high_4msbin.mat');
mean_perf_sm_ensemble_high_4msbin = temp.mean_perf_sm_ensemble_high;
temp = load('F:\EPHYS\Passive\MMA071219\mean_perf_sm_ensemble_high_10msbin.mat');
mean_perf_sm_ensemble_high_10msbin = temp.mean_perf_sm_ensemble_high;
temp = load('F:\EPHYS\Passive\MMA071219\mean_perf_sm_ensemble_high_20msbin.mat');
mean_perf_sm_ensemble_high_20msbin = temp.mean_perf_sm_ensemble_high;
%%
matrix_IC = [];
matrix_IC = horzcat(mean_perf_sm_ensemble_high_1msbin(:,1),mean_perf_sm_ensemble_high_4msbin(:,1),mean_perf_sm_ensemble_high_10msbin(:,1));

matrix_MGB = [];
matrix_MGB = horzcat(mean_perf_sm_ensemble_high_1msbin(:,2),mean_perf_sm_ensemble_high_4msbin(:,2),mean_perf_sm_ensemble_high_10msbin(:,2));

matrix_ACtx = [];
matrix_ACtx = horzcat(mean_perf_sm_ensemble_high_1msbin(:,3),mean_perf_sm_ensemble_high_4msbin(:,3),mean_perf_sm_ensemble_high_10msbin(:,3));

figure();
subplot(1,3,1);
colormap(hot);
imagesc(matrix_IC);
colorbar();
caxis([0.3 1])
subplot(1,3,2);
colormap(hot);
imagesc(matrix_MGB);
colorbar();
caxis([0.3 1])
subplot(1,3,3);
colormap(hot);
imagesc(matrix_ACtx);
colorbar();
caxis([0.3 1])
%%
figure();
plot(perf_complete,'k-o','linewidth',1.5);
hold on
plot(perf1,'r:o','linewidth',1.5);
hold on
plot(perf2,'r-o','linewidth',1.5);
hold on
% plot(perf3,'r:o','linewidth',1.5);
% hold on
plot(perf_su,'b-o','linewidth',1.5);
set(gca,'xtick',[1,2,3]);
set(gca,'xticklabel',{'IC','MGB','ACtx'});
xlim([0 4]);
legend('Complete ensemble','Ensemble of rate/mixed coding units','Ensemble of temporal coding units','Single units');
box off
title('PSTH Classifier performance');