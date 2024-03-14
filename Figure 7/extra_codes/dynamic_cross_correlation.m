analyze_sig = 'cross_signal';
location = {'NDB', 'SI/GP'}
type = {'RAND1','REG', 'RAND'}
fiber = {'fiber1', 'fiber2'};
fiberid = 2
clear data_ana
for i = 1:6
    for kk = 1:3
        numFrame = size(data(1).REG.cross_signal.fiber2{1},1);
        data_ana(i).(type{kk}).(analyze_sig) = cell([size(data(i).(type{kk}).(analyze_sig).(fiber{fiberid}))]);
        for j = 1:size(data_ana(i).(type{kk}).(analyze_sig), 1)
            for z = 1:size(data_ana(i).(type{kk}).(analyze_sig),2)
                if strcmp(analyze_sig, 'avg_signal')
                data_ana(i).(type{kk}).(analyze_sig){j,z} = max(data(i).(type{kk}).(analyze_sig).(fiber{fiberid})(j, z));
                else
                 range_idx = 800:1200;
                [data_ana(i).(type{kk}).(analyze_sig){j,z}, I] = max(data(i).(type{kk}).(analyze_sig).(fiber{fiberid}){j, z}(:,range_idx),[], 2);
                I
                end
                
            end
        end
    end
end

%%
cyc = 4
switch cyc
    case 4
        idx = 1
    case 12
        idx = 5
end
reg = []
control =[];
for i = 1:length(data_ana)
    set1 = [data_ana(i).RAND1.cross_signal{idx,1}; data_ana(i).REG.cross_signal{idx,1}; data_ana(i).RAND.cross_signal{idx,1}];
    set2 = [data_ana(i).RAND1.cross_signal{idx,2}; data_ana(i).REG.cross_signal{idx,2}; data_ana(i).RAND.cross_signal{idx,2}];
    set3 = [data_ana(i).RAND1.cross_signal{idx,2}; data_ana(i).REG.cross_signal{idx,3};data_ana(i).RAND.cross_signal{idx,3}];
    reg(:,i)  = set1;
    control(:,i) = [data_ana(i).RAND1.cross_signal{idx,4}; data_ana(i).REG.cross_signal{1,4};data_ana(i).RAND.cross_signal{1,4}];
    
end
%%
% figure;
% figure
hold on
switch cyc
    case 4
        idx = 5;
    case 12
        idx = 15;

end

plot_value = control;
color =cbrewer('div', 'RdYlBu',4);
rectangle('Position',[0.5,-0.1, idx, 0.5], 'FaceColor',[color(3,:),0.5], 'EdgeColor',[color(3,:),0.5])
hold on
rectangle('Position',[idx+0.5,-0.1, idx, 0.5], 'FaceColor',[0.5,0.5,0.5,0.5], 'EdgeColor',[0.5,0.5,0.5,0.5])
% for i = 1:size(reg,2)
%     h1 = plot(reg(:,i),'-','Color','k')
%     h1.Color(4) = 0.5
%     hold on
%     h2 = scatter(1:size(reg,1),reg(:,i), 'MarkerFaceColor', 'k','MarkerEdgeColor', 'k')
%     h2.MarkerFaceAlpha =0.5
%     h2.MarkerEdgeAlpha = 0.5
% end
xlim([-2,idx*2+1])
ylim([-0.1, 0.5])
hold on
mean_value = mean(plot_value,2);
sem_value  = std(plot_value,0,2)/sqrt(size(plot_value,6));
h2 = errorbar((1:length(plot_value))-2,mean_value,sem_value, '-d','color',color(4,:), 'MarkerFaceColor', color(4,:), 'LineWidth', 1, 'CapSize',12)
xlabel('Frames (20 noise bursts)')
ylabel('Peak Cross-correlation')
% text(2,0.42, 'RAND', 'FontSize', 12)
% text(7,0.42, 'RAND', 'FontSize', 12)
title([location{fiberid} 'Cycle Size ', num2str(cyc)])
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,1200,600])
plot([5.5,5.5], [-0.1, 0.4], '--','color',[0.5,0.5,0.5])
%%
legend([h1,h2], {'REG-RAND','RAND-RAND'})