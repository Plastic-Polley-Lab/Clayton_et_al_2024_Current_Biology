i =1  % 1, 3, 7, 8, 10, 14
figure; plot(data(1).bin5.time, data(i).bin5.psth_avg)
hold on;
j = 1
onset  = summary.ko_psth_pop.exc_cp.firstOnset(j)*1000;
offset = summary.ko_psth_pop.exc_cp.firstOffset(j)*1000;
lastOffset = summary.ko_psth_pop.exc_cp.lastOffset(j)*1000;
y_lim = 10
plot([onset, onset], [0, y_lim], '--r' )
plot([offset, offset], [0, y_lim], '--r')
plot([lastOffset, lastOffset], [0, y_lim], '--b')
xlabel('Time (ms)')
ylabel('FR (Hz)')
%%
i = 9 % 2, 4, 5, 6, 9
figure; plot(data(1).bin5.time, data(i).bin5.psth_avg)
hold on;
j = 5
onset  = summary.ko_psth_pop.inh_cp.firstOnset(j)*1000;
offset = summary.ko_psth_pop.inh_cp.firstOffset(j)*1000;
lastOffset = summary.ko_psth_pop.inh_cp.lastOffset(j)*1000;
y_lim = 10
plot([onset, onset], [0, y_lim], '--r' )
plot([offset, offset], [0, y_lim], '--r')
plot([lastOffset, lastOffset], [0, y_lim], '--b')
xlabel('Time (ms)')
ylabel('FR (Hz)')