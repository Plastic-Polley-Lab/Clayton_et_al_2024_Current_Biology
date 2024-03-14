function [t_sequencies, t_sequence, t] = complex_sound(intervals, reg_interval, jitters, set_num, CycleLength, outputFolder)
% intervals = 100:10:300;
% reg_interval = [260,120,190,230];
% outputFolder = 'Z:\File Transfer\KeChen\Jitter_Sound_pattern_Blaise\'
% set_num = 1;
Variable.CycleLength = CycleLength;
Variable.NumCycles   = 25;
% jitter_type = 'jitter-rand';
[t_randJitter, jitters] = generation_pattern(reg_interval, intervals, jitters, Variable, outputFolder, 'jitter-rand');
% jitter_type = 'jitter-scale';
[t_scale, jitters] = generation_pattern(reg_interval, intervals,jitters, Variable, outputFolder, 'jitter-scale');

% t_pre_rand = reg_interval(randi(length(reg_interval), [1 50]));
% t_reg = [t_pre_rand, repmat(reg_interval, 1, Variable.NumCycles)];
t_reg = repmat(reg_interval, 1, Variable.NumCycles);

t_rand = [];
for k = 1:Variable.NumCycles
    t_rand = [t_rand reg_interval(randperm(length(reg_interval)))];
end

t = [t_randJitter; t_scale; t_reg; t_rand];
n_pattern = size(t, 1);
t_sequence = randperm(n_pattern);
t_sequencies = []
for i = 1: length(t_sequence)
    t_sequencies = [t_sequencies, t(t_sequence(i), :)];
end


% t_pre_rand = reg_interval(randi(length(reg_interval), [1 50]));
% t_rand = [t_pre_rand, t_rand];

plot_patternSeq_gui([t_reg; t_randJitter], 4, [0,jitters], 'Jitter')
plot_patternSeq_gui([t_scale; t_rand], 4, [jitters,1], 'Scale')
% set(gcf,'position',[900,200,800,800])
