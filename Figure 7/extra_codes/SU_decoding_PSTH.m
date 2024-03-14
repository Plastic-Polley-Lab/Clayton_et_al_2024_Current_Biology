% Decode single units, find selectivity and decoding latencies
% Plot results
%
% Luca Mazzucato August 2017
function [results, options] = SU_decoding_PSTH(data)
% addpath('E:\Ke_Chen\MATLAB\Ephys-analysis\NeuralDecoder')
%%
% data = summary(9);
trig=fieldnames(data);
for i = 1: length(trig)
    spikes.(trig{i}) = data.(trig{i}).raster'; 
end
nev=numel(trig);
%%
% decoding parameters
kind='units';
xval=1; % number of leave-out trials for each x-validation run
binsize=1; % size of decoding window; 1 ms
% windows=[0 2.5]; % total trial interval to be decoded aligned at t=0.
windows=[251 1000]; % total trial interval to be decoded aligned at t=250.
nboot=1000; % number of bootstrap runs for each shuffled stimulus
filesave_dec='decoding'; % prefix of files to be saved
%%
[results, options] = fun_decode_units_selectivity(spikes,binsize,windows,xval,nboot,trig);

%%
%--------
% PLOT
%--------
% numfig=1;
% fun_decode_units_selectivity_plot(results,trig);