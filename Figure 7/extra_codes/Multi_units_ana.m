clearvars;
%% Server
%% Parameters
blockName = 'KeC030920-8-4';
tankName = 'F:\KeChen\RawData\Rach_Ephys\KeC022820\Raw_Waveform';


%% Read
trace = double(results.waveform)*1e6;
trace2 = double(results.waveform)*1e6;

Fs = 24414;
d1 = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',Fs);

d2 = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',119,'HalfPowerFrequency2',121, ...
               'DesignMethod','butter','SampleRate',Fs);

trace2 = filtfilt(d,trace); % remove 60Hz noise
trace3 = filtfilt(d2,trace2); % remove 120Hz noise
          
d4 = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',0.1,'HalfPowerFrequency2',500, ...
    'SampleRate',Fs,'DesignMethod','butter');    
lfp = filtfilt(d4,trace3);

[pxx,f] = pwelch(trace(1:240000));
[pxx2,f2] = pwelch(trace3(1:240000));
[pxx3,f3] = pwelch(lfp(1:240000));
figure;
plot(f/(2*pi)*Fs,10*log10(pxx))
xlim([0,5000])
hold on
plot(f2/(2*pi)*Fs,10*log10(pxx2))
plot(f3/(2*pi)*Fs,10*log10(pxx3))

figure;
plot(lfp(2:240000))

d = designfilt('bandpassiir','FilterOrder',4, ...
               'HalfPowerFrequency1',300,'HalfPowerFrequency2',3000, ...
               'DesignMethod','butter','SampleRate',Fs);
       
trace3 = filtfilt(d,trace3); % high pass [300,3000]     

figure;plot(trace(1:24000))
electrode3 = trace3;
eletrode2  = trace3;
electrode4 = trace3;

ref = 1/3*(electrode3+eletrode2+electrode4);
trace2 = zeros(size(trace));
ft = fft(trace);
ft(61) = 0;
ft(121) = 0;
ft(181) = 0;
ft(241) = 0;
ft(24355) = 0;
ft(24295) = 0;
ft(24235) = 0;
ft(24175) = 0;
trace2 = ifft(ft);

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

load YBMap;

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


