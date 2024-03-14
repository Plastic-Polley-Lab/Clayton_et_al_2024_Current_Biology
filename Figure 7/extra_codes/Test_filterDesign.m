%% filter design practise
% Please check signal processing toolbox tutorial in matlab
%FIR filters
Fstop1 = 10;
Fpass1 = 15;
Fpass2 = 20;
Fstop2 = 25;
Astop1 = 65;
Apass  = 0.5;
Astop2 = 65;
Fs = 250;

d = designfilt('bandpassfir', ...
  'StopbandFrequency1',Fstop1,'PassbandFrequency1', Fpass1, ...
  'PassbandFrequency2',Fpass2,'StopbandFrequency2', Fstop2, ...
  'StopbandAttenuation1',Astop1,'PassbandRipple', Apass, ...
  'StopbandAttenuation2',Astop2, ...
  'DesignMethod','equiripple','SampleRate',Fs);

fvtool(d)

impz(d,50)
figure
y2 = filtfilt(d,ecg);  % zero-phase implementation - delay compensation
plot(t,ecg);
hold on
% plot(t,y1,'r','linewidth',1.5);
plot(t,y2,'r','linewidth',1.5);
[pxx,f] = pwelch(y2,250);
figure;plot(f,pxx)
hold on
[pxx2,f] = pwelch(ecg,250);
plot(f,pxx2)
%%
dbutter = designfilt('bandpassiir','FilterOrder',10, ...
    'HalfPowerFrequency1',2,'HalfPowerFrequency2',6, ...
    'SampleRate',250,'DesignMethod','butter');

%%
dbutter2 = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',2,'HalfPowerFrequency2',6, ...
    'SampleRate',250,'DesignMethod','butter');
dbutter3 = designfilt('bandpassiir','FilterOrder',6, ...
    'HalfPowerFrequency1',300,'HalfPowerFrequency2',5000, ...
    'SampleRate',24414.0625,'DesignMethod','butter');
y3 = filtfilt(dbutter3,double(results.waveform(1:200000))); 
fvtool(dbutter,dbutter2,dbutter3)