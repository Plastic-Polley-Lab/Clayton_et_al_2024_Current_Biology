[FileName,PathName] = uigetfile('F:\KeChen\RawData\Rach_Ephys\*.mat','Select the FRA Impale file')
path = [PathName,FileName];
cd(PathName)
chan = 1;
[fra, freqs, spls] = fra_impale(path);
figure(100)
for i = 1:16
    subplot(4,4, i)
    fra_plot(fra, freqs, spls, i)
end
set(gcf,'position',[100,100,1000,800])
set(gcf, 'Color', 'w')
export_fig('FRA_multiunits-1',  '-pdf')

figure
for i = 17:32
    subplot(4,4, i-16)
    fra_plot(fra, freqs, spls, i)
end
set(gcf,'position',[100,100,1000,800])
set(gcf, 'Color', 'w')
export_fig('FRA_multiunits-2',  '-pdf')
%%
figure
for i = 33:48
    subplot(4,4, i-32)
    fra_plot(fra, freqs, spls, i)
end
set(gcf,'position',[100,100,1000,800])
set(gcf, 'Color', 'w')
export_fig('FRA_multiunits-3',  '-pdf')

figure
for i = 49:64
    subplot(4,4, i-48)
    fra_plot(fra, freqs, spls, i)
end
set(gcf,'position',[100,100,1000,800])
set(gcf, 'Color', 'w')
export_fig('FRA_multiunits-4',  '-pdf')
%%