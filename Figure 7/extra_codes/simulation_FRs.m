 tau1 = 1/20;
 tau2 = 1/40;
 ti   = 1;
 tf   = 250;
 t        = linspace(ti,tf,250);
 yinitial = exp(-tau1*t);
 yfinal   = exp(-tau2*t);
 figure
 plot(t,yinitial,'r*') 
 hold on
 plot(t,yfinal,'c*')
 hold of

figure
k = 1
y = exp(k*(1:1:20))/max(exp(k*(1:1:20)));
plot(y)

data1 = [y,yinitial];
data2 = [y, yfinal];

figure
plot(data1,'-r*')
hold on
plot(data2,'-c*')
hold of
[normalizedACF1, lags] = autocorr(data1,249);
[normalizedACF2, lags] = autocorr(data2,249);
figure
subplot(2,1,2)
h1 = plot(lags, normalizedACF1,'-r*')
hold on
h2 = plot(lags,normalizedACF2,'-c*')
plot([0, 250], [0,0])

zerocross1 = min(find(normalizedACF1<0));
zerocross2 = min(find(normalizedACF2<0));
h3 = plot([zerocross1, zerocross1],[-0.1, 1], '-r')
h4 = plot([zerocross2, zerocross2],[-0.1, 1], '-c')
xlim([0,250])
ylabel('ACF')
xlabel('Time (ms)')
legend([h1, h2, h3, h4], {'\tau = 20 ms', '\tau = 40 ms', 'cross zeros after 61 ms', 'cross zeros after 78 ms'})

subplot(2,1,1)
h1 = plot(data1,'-r*')
hold on
h2 = plot(data2,'-c*')
h3 = plot([zerocross1, zerocross1]+21,[-0.1, 1], '-r')
h4 = plot([zerocross2, zerocross2]+21,[-0.1, 1], '-c')
h5 = plot([0, 250], [mean(data1),mean(data1)], '-m')
h6 = plot([0, 250], [mean(data2),mean(data2)], '-b')
xlim([0,250])
legend([h1, h2, h3, h4, h5, h6], {'simulation FR with \tau = 20 ms', 'simulation FR with \tau = 40 ms', 'cross zeros after 61 ms + 21 ms(peak firing rate)', 'cross zeros after 78 ms + 21 ms(peak firing rate)', ...
    'mean firing rate for \tau = 20 ms', 'mean firing rate for \tau = 40 ms' })
ylabel('Firing rate')
xlabel('Time (s)')

[~,ind]=max(diff(find(data1<mean(data1))));

find(data2<mean(data2))
