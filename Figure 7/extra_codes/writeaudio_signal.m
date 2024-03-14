function signal = writeaudio_signal(t, Fs)
%%
% NOISEPATTERN -- creates random or regular pattern for Tosca.
% the transition of regular to random starts from swap of twFso/three
% interval
% Determine where sound generation folder is.
folder = 'F:\KeChen\MATLAB\MATLAB-Noise'; 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));
% clearvars -except Fs args

args = '';
% Begin parameters
Fixed.StimDuration_ms = 20;
Fixed.StimRamp_ms = 5;

% Variable.IsRegular = 1;
Variable.IsRegular = 0;
% Variable.CycleLength = CycleLength; % define how many elements in each cycle
Variable.NumCycles = 25;
Variable.Seed = -1;
% Variable.swap = 2;

% Variable.Seed = 1;
% End parameters

% set_num = 3;
% Private parameters
frameSize_s = 0.5;
% outputFolder = 'F:\Behavior\Patterns\test';
% calFile = fullfile('F:\Behavior\Patterns', 'Freefield.calib');
% outputFolder = output;
% outputFolder = 'F:\KeChen\Behavior\Patterns\test';
calFile = fullfile('F:\KeChen\Behavior\Patterns', 'Freefield.calib');

% Do it
result = '';
% output = '';
% load('F:\KeChen\Behavior\Patterns\frozen_noise.mat');

try
   % Create gate
   dac.TotalDuration = frameSize_s * 1000;
   dac.SampleRate = Fs;
   g = Gate();
   g.IsActive = true;
   g.RiseFallTime = Fixed.StimRamp_ms;
   g.Width = Fixed.StimDuration_ms;
   win = create(g, dac);
   load('F:\KeChen\Behavior\Patterns\frozen_noise.mat');
%    frozen_noise = ones(1, length(win))';
%    frozen_noise = normrnd(0, 1, length(win), 1); %save this and keep and
% try to generate pure tone
% tt = 0:1/Fs: 0.020;
% tt = tt(1:end-1);
% f = 2000;
%    puretone = sin(2*pi*f*tt); % pure tone for 20 ms
% frozen_noise = puretone';

   %reuse it
   
%%

    
  %%    
   % Create signal
   duration = t(end)-t(1);
   npts = round(duration * Fs);
   
   signal = zeros(npts, 1);
   nw = length(win);
   starts = round((t-t(1))*Fs);
   for k = 1:length(starts)-1
      signal(starts(k) + (1:nw)) = win .* frozen_noise (1:length(win));
   end

   % Compute SPL reference
   sf = max(abs(signal));
   refSPL = compute_noise_maxSPL(calFile, Fs);
   refSPL = refSPL - 20*log10(sf);
   signal = signal / sf;

   % Write wav files
%       EPLwavwrite(signal, Fs, 16, outputFile,  'refSPL', refSPL);

%    EPLwavwrite(signal, Fs, 16, fullfile(outputFolder, sprintf('final_cyc%d_reg_set%d_%s%d.wav',Variable.CycleLength, set_num, jitter_type,jitters*100)),  'refSPL', refSPL);
%save output
%    save(fullfile(outputFolder, sprintf('final_cyc%d_reg_set%d_%s%d.mat',Variable.CycleLength, set_num, jitter_type,jitters*100)),'t');
 % save(sprintf('final_cyc%d_regrand_set_para%d.mat',Variable.CycleLength,set_num),'reg_interval','swap_interval');
   % Plot results
%    ti = (0:npts-1) / Fs;
%    figure;
%    clf;
%    hold on;
%    plot(ti, signal);
%    plot(ti, reward + 1, 'r');
   
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
