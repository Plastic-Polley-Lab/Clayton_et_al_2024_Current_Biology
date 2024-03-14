clearvars;
%% Server
addpath(genpath('F:\KeChen\MATLAB\MATLAB-Noise'));
if ~exist('TT','var')
    clear;

    TT = actxserver('TTank.X'); % Tank50
    if TT.ConnectServer('Local', 'Impale') == 0,
        error('Cannot connect to server.');
    end
end

%% Parameters
path = uigetdir('F:\KeChen\RawData\Rach_Ephys\','Select Noise CSD recording');
path_str = split(path,'\');
blockName = char(path_str(end));
tankName = char(join(path_str(1:end-1),'\'));

Fs = 24414;
trialNum = 50;%how many trials


% close all;

% chans = [9 8 10 7 13 4 12 5 15 2 16 1 14 3 11 6]+16;
chans = [43 41 39 40 42 6 34 45 33 47 36 37 35 44 48 46 4 15 11 9 7 5 3 1 38 8 10 12 16 14 13 2 30 ...
      17 21 23 25 27 29 31 60 26 24 22 18 20 19 32 53 55 57 58 56 28 64 51 63 49 62 59 61 54 52 50];
% chans = 1:64;
% chans = chans([1:end-8 end-6:end]);%removing the one channel with high impedance
% chans(53) = chans(52);
% chans(57) = chans(56);
chans = chans([1:32]);%shankA
% chans = chans([1:32]+32);%shankB
% chans = chans([1:7,9:16,18:32]);%shankB

%% Read
ev = read_ticks_for_rach(blockName,tankName, TT);
ev = ev - ev(1);
trialLength = round((ev(2)-ev(1))*1000)/1000;

blockSize = round(50/trialLength);%????

trace = [];

disp('  Reading traces...');
for i = 1:length(chans)
    disp(['    Chan ',num2str(chans(i))]);
    temp = [];
    for j = 1:ceil(trialNum/blockSize)
        temp = [temp;read_waves_for_rach(blockName,chans(i),[(j-1)*blockSize+1 min(j*blockSize,trialNum)], ... 
            tankName, TT)];
    end
%     trace(1:size(temp,1),i) = double(temp)/4; % if it's int
    trace(1:size(temp,1),i) = double(temp) * 1e6; %wei change it to match the setting 16 float
    %trace(1:size(temp,1),i) = double(temp);% * 1e6;
%     trace(:,i) = trace(:,i)-mean(trace(1:Fs*10,i));
%     trace(:,i) = trace(:,i)-mean(trace(1:end,i));
end

%bad channel
% trace(:,7) = mean(trace(:,[6,9]),2);
trace(:,8) = mean(trace(:,[7,9]),2);

trace(:,17) = mean(trace(:,[16,18]),2);
trace(:,32) = mean(trace(:,[30,31]),2);
% trace(:,40) = mean(trace(:,[39,41]),2);
%shankB
% trace(:,40-32) = mean(trace(:,[39,41]-32),2);
% trace(:,48-32) = mean(trace(:,[47,49]-32),2);
% trace(:,46-32) = mean(trace(:,[45,47]-32),2);
% trace(:,44-32) = mean(trace(:,[43,45]-32),2);

disp('  60Hz filtering...');
trace2 = zeros(size(trace));
for j = 1:floor(size(trace,1)/Fs)
    for k = 1:size(trace,2)
        ft = fft(trace((1:Fs)+(j-1)*Fs,k));
        ft(61) = 0;
        ft(121) = 0;
        ft(181) = 0;
        ft(241) = 0;
        ft(24355) = 0;
        ft(24295) = 0;
        ft(24235) = 0;
        ft(24175) = 0;
        trace2((1:Fs)+(j-1)*Fs,k) = ifft(ft);
    end;
end
trace = trace2;
clear trace2;

% % dTrace = trace-mean(trace(:,1:64),2)*ones(1,length(chans));
% dTrace = trace-mean(trace(:,1:length(chans)),2)*ones(1,length(chans));
% 
% % dTrace = trace;
% 
% disp(['  Filtering data...']);
% Fs = 24414;
% LowCutoff = 300;
% [B,A] = butter(2,LowCutoff/(Fs/2),'high');
% HighCutoff = 3000;
% [D,C] = butter(2,HighCutoff/(Fs/2),'low');
%   
% MUtrace = filter(B,A,dTrace);
% MUtrace = filter(D,C,MUtrace);


%% #ATL; Below is for CSD plotting
disp('  Downsampling LFP...');
newFs = 1000;
LStrace = resample(trace, newFs, Fs);


disp(['  Filtering LFP...']);
LowCutoff = 0.1;
[B,A] = butter(4,LowCutoff/(newFs/2),'high');

LStrace = filtfilt(B,A,LStrace);


LStrace = [LStrace;zeros(1000,length(chans))];%%???????

% LStrace = LStrace*0.5+LStrace(:,[1,1:end-1])*0.25+LStrace(:,[2:end,end])*0.25;
% LStrace = LStrace*0.5+LStrace(:,[1,1:end-1])*0.25+LStrace(:,[2:end,end])*0.25;

LSTrace = filtfilt(hann(5),sum(hann(5)),LStrace')';


meanTrace = [];

traceRaster = [];
diffThr = 1e3;
rejThr = 5e3;
for i = 1:size(trace,2);
    goodTrials = [];
    for j = 1:trialNum
        traceRaster(i,j,:) = LStrace(round((ev(j)*newFs))+(1:newFs*trialLength),i);
        if (max(abs(diff(LStrace(round((ev(j)*newFs))+(1:newFs*trialLength),i)))) < diffThr)...
            && (max(abs(LStrace(round((ev(j)*newFs))+(1:newFs*trialLength),i))) < rejThr)
            goodTrials = [goodTrials,j];
        end
    end
    meanTrace(:,i) = mean(traceRaster(i,goodTrials,:),2);
end

disp(size(goodTrials));

%meanTrace = squeeze(mean(traceRaster(:,1:8:end,:),2))';

meanCSD = [];
CSDRaster = diff(traceRaster,2,1);%second derivative along the 1st dimension
for i = 1:size(CSDRaster,1)
meanCSD(:,i) = mean(CSDRaster(i,goodTrials,:),2);
end

disp('  Calculating CSD...');
intERP = [];
for j = 1:size(meanTrace,1);
    intERP(j,:) = interp1(1:size(meanTrace,2),meanTrace(j,:),1:0.1:size(meanTrace,2),'spline');
end;

% ERPCSD = [];
% for i = 1:size(intERP,2)-21;
%     ERPCSD(:,i) = intERP(:,i)+intERP(:,i+20)-2*intERP(:,i+10);
% end
% ERPCSD = ERPCSD/(1000*(0.1)^2);

ERPCSD = [];
for i = 1:size(intERP,2)-81; % change by ke, orginal it is 100; as the distance between site is 25 um not 20 um
    ERPCSD(:,i) = intERP(:,i)+intERP(:,i+80)-2*intERP(:,i+40);
end
ERPCSD = ERPCSD/(1000*(0.1)^2);

% load YBMap;
load('YBMap.mat')
%% Plotting

figure('units','normalized','outerposition',[0 0.05 1 0.95]);
set(gcf,'PaperPositionMode','auto');

subplot(211);
imagesc(-ERPCSD');
maxval = max(abs(ERPCSD(:)));
caxis([-maxval maxval]);
colormap(YBMap);
h = colorbar;
% oldticks = get(h,'yticklabel');
% set(h,'yticklabel',oldticks(end:-1:1,:));
ylabel(h,'mV/mm^2','fontsize',16);
%set(gca,'FontSize',16);
xlabel('Time (ms)');
ylabel('Channels');
% set(gca,'ytick',1:100:430);
% set(gca,'yticklabel',11:10:64);
title('CSD');

subplot(212);
imagesc(-intERP');
maxval = max(abs(intERP(:)));
caxis([-maxval maxval]);
colormap(YBMap);
h = colorbar;
oldticks = get(h,'yticklabel');
set(h,'yticklabel',oldticks(end:-1:1,:));
ylabel(h,'\muV','fontsize',16);
%set(gca,'FontSize',16);
xlabel('Time (ms)');
ylabel('Channels');
set(gca,'ytick',1:20:(length(chans)-1)*10);
set(gca,'yticklabel',1:2:length(chans)-1);
title('LFP');
% %%
% %LFP
% figure(2);
% for i = 1:size(meanTrace,2)
% plot(meanTrace(:,i)-(i-1)*100);
% hold on
% end
% %CSD
% figure(3);
% for i = 1:size(meanCSD,2)
% plot(meanCSD(:,i)-i*40);
% hold on
% end
% %%
% %clear;
% 
% figure('units','pixels','outerposition',[100 100 700 600]);
% set(gcf,'PaperPositionMode','auto');
% winLeng = 20;
% hold on;
% for i = 1:length(chans)
% plot(1/Fs:1/Fs:winLeng,dTrace(1:Fs*winLeng,i) - i*10000);
% end
% hold off;  set(gca,'ytick',[]);
% set(gca,'FontSize',16);
% xlabel('Time (s)');
% 
% figure('units','pixels','outerposition',[100 100 700 600]);
% set(gcf,'PaperPositionMode','auto');
% winLeng = 10;
% hold on;
% for i = 1:16
%    subplot(211);
%    plot(1/Fs:1/Fs:winLeng,MUtrace((1:Fs*winLeng)+30*Fs,i) - i*200);
%    set(gca,'FontSize',16);
%    title('Channel 16');
%    set(gca,'ytick',[]);
%    subplot(212);
%    plot(1/Fs:1/Fs:0.2,MUtrace((1:Fs*0.2)+30*Fs,i) - i*150);
% end
% hold off;
% set(gca,'ytick',[]);
% set(gca,'FontSize',16);
% xlabel('Time (s)');


