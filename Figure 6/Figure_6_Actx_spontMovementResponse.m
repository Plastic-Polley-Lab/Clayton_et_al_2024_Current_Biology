%Clayton et al, Figure 6, panels C, D, E

%Input data: spontaneous facial motion bouts (-500 ms to +1500 ms of the
%facial motion trace) and corresponding single-unit neural recordings from
%all layers of ACtx.


function [] = Actx_spontMovementResponse()

%CHANGE THIS PATH!
file_dir = 'C:\Users\AY92\Documents\Data\Orofacial Movements ACtx\Fig6_A1_ephys';
cd(file_dir)

R = dir('arc*.mat');
for i = 1:length(R)
    load(R(i).name);
    movementResp(i).movement = recording.movement;
    movementResp(i).units = recording.units;
end

% %calculate the response per cell to each movement bout (concatenating
% %across recording blocks from the same recording)
cellDepth = [];
cellType = [];
t = 1;
for i = 1:length(movementResp)
    for j = 1:length(movementResp(i).units)
        cellDepth(t) = movementResp(i).units(j).depth;
        cellType(t) = movementResp(i).units(j).spikeWidth;
        t = t + 1;
    end
end

%extract stillness and movement spikes per cell

t = 1;
for f = 1:length(movementResp)
    units = movementResp(f).units;
    for i = 1:length(units)
        tempStill = sum(units(i).raster(:,1:250),2); %sum of spikes in 'stillness' window
        tempMove = sum(units(i).raster(:,501:750),2);%sum of spikes in 'movement' window
        %test for significant modulation
        if nansum(nansum(units(i).raster(:,1:750),2)>0) > 0.5*length(tempMove)
            [~,tempSig] = ttest(tempStill,tempMove);
            summary.sigMod(t) = tempSig<0.001;
        else
            summary.sigMod(t) = false;
        end
        summary.meanResponse(t,:) = mean(units(i).raster);
        summary.moveResp(t) = mean(tempMove);
        summary.stillResp(t) = mean(tempStill);
        summary.cellDepth(t) = cellDepth(t);
        summary.cellType(t) = cellType(t);
        t = t + 1;
        clear tempStill tempMove
    end
end

%%plot neurograms, separated by depth and z-scored

superficial_response = summary.meanResponse(summary.cellDepth>=0,:);
deep_response = summary.meanResponse(summary.cellDepth<0,:);

%100 ms bin
binWidth = 100;
binEdges = round(1:binWidth:size(superficial_response,2));
for i = 1:length(binEdges)-1
    binnedSuperficial(:,i) = nansum...
        (superficial_response(:,binEdges(i):binEdges(i)+binWidth-1),2);
    binnedDeep(:,i) = nansum...
        (deep_response(:,binEdges(i):binEdges(i)+binWidth-1),2);
end

%sort by max movement response
for i = 1:size(binnedSuperficial,1)
    binnedSuperficial(i,:) = zscore(binnedSuperficial(i,:));
    max_superficial(i,1) = max(binnedSuperficial(i,:));
end

max_superficial(:,2) = 1:size(binnedSuperficial,1);
max_superficial = sortrows(max_superficial,1);

for i = 1:size(binnedDeep,1)
    binnedDeep(i,:) = zscore(binnedDeep(i,:));
    max_deep(i,1) = max(binnedDeep(i,:));
end

max_deep(:,2) = 1:size(binnedDeep,1);
max_deep = sortrows(max_deep,1);

% %generate plots
figure;imagesc(linspace(-500,1500,19),[],binnedSuperficial(flipud(max_superficial(:,2)),:));
box off
caxis([-3 3])
colorbar
title('RS cells, superficial layers')
ylabel('Cell')
xlabel('Time re: start of facial movement bout')
figure;imagesc(linspace(-500,1500,19),[],binnedDeep(flipud(max_deep(:,2)),:));
box off
caxis([-3 3])
colorbar
title('RS cells, deep layers')
ylabel('Cell')
xlabel('Time re: start of facial movement bout')


%%plot horizontal error bar with response by depth
for i = 1:length(summary.moveResp)
    logResp_movement(i) = log(summary.moveResp(i))-log(summary.stillResp(i));
end
logResp_movement(logResp_movement == Inf | logResp_movement == -Inf) = NaN;

%restrict to RS cells only
logResp_RS = logResp_movement(cellType == 2);
cellDepth = cellDepth(cellType == 2);

figure;
depths = -600:50:400;
logResp_RS_abs = abs(logResp_RS);
for i = 1:length(depths)-1
x = logResp_RS_abs(cellDepth>=depths(i) & cellDepth<depths(i+1));
hold on
errorbar(nanmean(x),([depths(i+1)]),nansem(x'),'horizontal','ok');
end
box off
set(gca,'TickDir','out')
set(gca,'FontSize',12)
ylabel('depth re: L4 sink')
xlabel('abs(log(movement)-log(stillness))')
xlim([0 1.6])

%plot pie chart
%upper layer RS units
cellDepth_RS = summary.cellDepth(summary.cellType==2);
sigMod_RS = summary.sigMod(summary.cellType==2);

deep_layer_modulation = logResp_RS(sigMod_RS & cellDepth_RS <= 0);
superficial_layer_modulation = logResp_RS(sigMod_RS & cellDepth_RS > 0);

X_superficial = [sum(superficial_layer_modulation<0) sum(superficial_layer_modulation>0)...
    sum(cellDepth_RS > 0)-sum(superficial_layer_modulation<0)-sum(superficial_layer_modulation>0)];

figure;subplot(1,2,1),pie(X_superficial); colormap(pink)
title('Upper layer RS units')

X_deep = [sum(deep_layer_modulation<0) sum(deep_layer_modulation>0) ...
    sum(cellDepth_RS <= 0)-sum(deep_layer_modulation<0)-sum(deep_layer_modulation>0)];
subplot(1,2,2);pie(X_deep); colormap(pink)
title('Deep layer RS units')







