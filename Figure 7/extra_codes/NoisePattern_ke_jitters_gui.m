function [t_jitter, t_reg, t_rand, jitters, reg_interval] = NoisePattern_ke_jitters_gui(set_num, CycleLength, jitter_type,output)
%%
% NOISEPATTERN -- creates random or regular pattern for Tosca.
% the transition of regular to random starts from swap of two/three
% interval
% Determine where sound generation folder is.
jitters = [5, 10, 15, 20];
jitters = jitters/100;
folder = 'E:\Ke_Chen\MATLAB\MATLAB-Noise';
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));
% clearvars -except Fs args
Fs = 100000;
args = '';
% Begin parameters
Fixed.StimDuration_ms = 20;
Fixed.StimRamp_ms = 5;
Fixed.Intervals_ms = '150:20:350';

% Variable.IsRegular = 1;
Variable.IsRegular = 0;
Variable.CycleLength = CycleLength; % define how many elements in each cycle
Variable.NumCycles = 35;
Variable.Seed = -1;
% Variable.swap = 2;

% Variable.Seed = 1;
% End parameters

eval(args);
% set_num = 3;
% Private parameters
frameSize_s = 0.5;
% outputFolder = 'F:\Behavior\Patterns\test';
% calFile = fullfile('F:\Behavior\Patterns', 'Freefield.calib');
outputFolder = output;
% outputFolder = 'F:\KeChen\Behavior\Patterns\test';
calFile = fullfile('E:\Ke_Chen\Behavior\Patterns', 'Freefield.calib');

% Do it
result = '';
output = '';
load('E:\Ke_Chen\Behavior\Patterns\frozen_noise.mat');
try
    % Create gate
    dac.TotalDuration = frameSize_s * 1000;
    dac.SampleRate = Fs;
    g = Gate();
    g.IsActive = true;
    g.RiseFallTime = Fixed.StimRamp_ms;
    g.Width = Fixed.StimDuration_ms;
    win = create(g, dac);
    
    %    frozen_noise = normrnd(0, 1, length(win), 1); %save this and keep and
    % try to generate pure tone
    % tt = 0:1/Fs: 0.020;
    % tt = tt(1:end-1);
    % f = 2000;
    %    puretone = sin(2*pi*f*tt); % pure tone for 20 ms
    % frozen_noise = puretone';
    
    %reuse it
    
    %%
    % Generate intervals
    intervals = eval(Fixed.Intervals_ms);
    reg_interval = randomSeq_noRepeats(intervals, Variable)
    t_reg = repmat(reg_interval, 1, Variable.NumCycles);
    %pre_rand from the reg_interval
    %now append randomness with cycleduration constraint
    t_rand = [];
    for k = 1:Variable.NumCycles
        t_rand = [t_rand reg_interval(randperm(length(reg_interval)))];
    end
    t_jitter = zeros(CycleLength-2, size(t_reg, 2));
    if strcmp(jitter_type, 'jitter-swap')
        t_jitter = interval_jitter_swap(t_jitter,reg_interval, Variable);
    else
        t_jitter = [];
        for i = 1:length(jitters)
            t_jitter(i,:) = interval_jitter(intervals, reg_interval, jitters(i), Variable, jitter_type);
        end
    end
    
    %%
    % Create signal
%     jitters = [0,jitters];
%     for i = 1:size(t_jitter,1)
%         t = t_jitter(i,:);
%         duration = ceil(0.001*(sum(t) + Fixed.StimDuration_ms) / frameSize_s) * frameSize_s;
%         npts = round(duration * Fs);
%         
%         signal = zeros(npts, 1);
%         nw = length(win);
%         starts = [0 cumsum(round(0.001*t*Fs))];
%         for k = 1:length(starts)-1,
%             signal(starts(k) + (1:nw)) = win .* frozen_noise;
%         end
%         
%         % Compute SPL reference
%         sf = max(abs(signal));
%         refSPL = compute_noise_maxSPL(calFile, Fs);
%         refSPL = refSPL - 20*log10(sf);
%         signal = signal / sf;
% 
%         % Write wav files
%         EPLwavwrite(signal, Fs, 16, fullfile(outputFolder, sprintf('final_cyc%d_reg_set%d_jitter%d%s.wav',Variable.CycleLength, set_num, jitters(i)*100,jitter_type)),  'refSPL', refSPL);
%         %save output
%         save(fullfile(outputFolder, sprintf('final_cyc%d_reg_set%d_jitter%d%s.mat',Variable.CycleLength, set_num, jitters(i)*100,jitter_type)),'t');
%         % save(sprintf('final_cyc%d_regrand_set_para%d.mat',Variable.CycleLength,set_num),'reg_interval','swap_interval');
%         % Plot results
%         %    ti = (0:npts-1) / Fs;
%         %    figure(1);
%         %    clf;
%         %    hold on;
%         %    plot(ti, signal);
%         %    plot(ti, reward + 1, 'r');
%     end
catch ex
    result = sprintf('%s Line %d: %s', ex.stack(1).name, ex.stack(1).line, ex.message);
end

%--------------------------------------------------------------------------
function SPL = compute_noise_maxSPL(calFile, Fs)

S = load_probetube_sens(calFile);
spec = S.sweep.spec;

signal = normrnd(0, 1, 1e5, 1);

sigspec = ffa(signal, Fs, [], 1, '', [0 Fs/2]);
calspec = interp1(spec, get(sigspec, 'Frequency'));

s = double(sigspec)/(length(signal)) .* double(calspec)';

SPL = 20*log10(sqrt(nansum(2*abs(double(s)).^2)));
