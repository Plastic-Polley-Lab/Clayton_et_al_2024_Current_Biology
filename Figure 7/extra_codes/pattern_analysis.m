%% frequency pattern
% neuron_num = 54;
% pattern_idx = 1:2:length(spikedata(1).clusterData.psth.raster);
% figure
% psth_idx = pattern_idx;
% psth = spikedata(neuron_num).clusterData.psth;
% psth1.raster = psth.raster(psth_idx );
% psth1.scmatrix = psth.scmatrix(psth_idx ,:);
% psth1.stimulus = psth.stimulus;
% 
% pattern_idx = 2:2:length(spikedata(1).clusterData.psth.raster);
% figure
% psth_idx = pattern_idx;
% psth = spikedata(neuron_num).clusterData.psth;
% psth2.raster = psth.raster(psth_idx );
% psth2.scmatrix = psth.scmatrix(psth_idx ,:);
% psth2.stimulus = psth.stimulus;
% 
% % plot the data
% psth_plot_freq(psth1,psth2,1:3000)
% set(gcf,'position',[100,200,1000,400])
%% noise pattern
% clear
keep = find([spikedata.keep] ==1);

for i = 1:length(keep)
    try
        neuron_num = keep(i);
        pattern_result(neuron_num) = plot_pattern_raster(neuron_num, spikedata, 1)
    catch
        warning('some errors happend, discard that units')
        continue
    end
end
save('summary_pattern.mat', 'pattern_result')

%% pull data together
%% get noise burst response
clear
mydir  = pwd;
idcs   = strfind(mydir,'\');
file_details = split(mydir, '_');
set_num = char(file_details(end-1))
newdir = [mydir(1:idcs(end)-1), '\NoiseCSD']
load([newdir, '\summary_noise.mat'])
resp_ind = find([spikedata.resp] == 1);
%%
batch_analysis = 1
% load('E:\Ke_Chen\Processed Data\Rach Recording\KeC16\ACtx\080520\regrand_cyc4_set101\summary_pattern.mat')
load('E:\Ke_Chen\Processed Data\Rach Recording\KeC20\090420\regrand_cyc4_set101\summary_pattern.mat')

iter = 1;
reg = [];
rand = [];
for i = 1:length(resp_ind)
    if isempty(pattern_result(resp_ind(i)).model)
    else
        reg(iter) = pattern_result(resp_ind(i)).model(1).tau_final;
        rand(iter) = pattern_result(resp_ind(i)).model(2).tau_final;
        iter = iter + 1;
    end
end

load('E:\Ke_Chen\Processed Data\Rach Recording\KeC17\ACtx\081720\regrand_cyc4_set101\summary_pattern.mat')
iter = 1;
reg2 = [];
rand2 = [];
for i = 1:length(resp_ind)
    if isempty(pattern_result(resp_ind(i)).model)
    else
        reg2(iter) = pattern_result(resp_ind(i)).model(1).tau_final;
        rand2(iter) = pattern_result(resp_ind(i)).model(2).tau_final;
        iter = iter + 1;
    end
end


figure;
h1 = scatter(rand, reg, 'ko', 'filled')
hold on
h2 = scatter(rand2, reg2, 'o', 'filled', 'MarkerFaceColor', [0.7,0.7,0.7], 'MarkerEdgeColor', [0,0,0])
plot([0,60],[0,60],'--', 'color', [0.5, 0.5, 0.5])
legend([h1, h2], {'KeC16', 'KeC17'})
axis square
xlabel('RAND tau (ms)')
ylabel('REG tau (ms)')
set(gca,'TickDir','out')
% set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,600,400,400])


reg = [reg, reg2];
rand = [rand, rand2];
[h, p] = ttest(reg',rand');

%%
