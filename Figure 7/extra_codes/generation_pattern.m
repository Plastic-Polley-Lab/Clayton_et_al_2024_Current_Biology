function [t_jitter, jitters] = generation_pattern(reg_interval, intervals, jitters, Variable, outputFolder, jitter_type)
% intervals = 100:10:300;
% reg_interval = [260,120,190,230];
% outputFolder = 'Z:\File Transfer\KeChen\Jitter_Sound_pattern_Blaise\'
% set_num = 1;
% Variable.CycleLength = 4;
% Variable.NumCycles   = 25;
% jitter_type = 'jitter-scale';

% jitter_type = 'jitter-rand';
% jitters = [5, 10, 15, 20]/100;
t_jitter = [];
for i = 1:length(jitters)
    t_jitter(i,:) = interval_jitter(intervals, reg_interval, jitters(i), Variable, jitter_type);
end
t_reg = repmat(reg_interval, 1, Variable.NumCycles);
% % t_jitter = [t_reg; t_jitter];
% jitters =[0, jitters];
% plot_patternSeq_gui(t_jitter, 4, jitters, jitter_type)

% k = [];
% for i = 1:length(jitters)
%     t_pre_rand = reg_interval(randi(length(reg_interval), [1 50]));
%     k(i,:) =t_pre_rand;
% end
% 
% t = [k, t_jitter];
% plot_patternSeq_gui(t, 4, jitters, jitter_type)
% tt = t;
% t = [];
% for i = 1:size(tt,1)
%     t = tt(i,:);
%     save(fullfile(outputFolder, sprintf('final_cyc%d_reg_set%d_%s%d.mat',Variable.CycleLength, set_num, jitter_type, jitters(i)*100)),'t');
% end

% t_rand = [];
% for k = 1:Variable.NumCycles
%     t_rand = [t_rand reg_interval(randperm(length(reg_interval)))];
% end
% t_pre_rand = reg_interval(randi(length(reg_interval), [1 50]));
% t_rand = [t_pre_rand, t_rand];
% t = t_rand;
% save(fullfile(outputFolder, sprintf('final_cyc%d_reg_set%d_%s%s.mat',Variable.CycleLength, set_num, jitter_type, '-rand')),'t')

%%
% file = dir('*.mat');
% for i = 1:length(file)
%     load(file(i).name)
%     temp(i,:) = t;
% end
% plot_patternSeq_gui(temp(6:10,:), 4, [0, 0.05, 0.1, 0.15, 0.2], jitter_type)