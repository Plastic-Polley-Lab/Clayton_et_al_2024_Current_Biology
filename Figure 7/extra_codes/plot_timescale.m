function plot_timescale(wt_psth_pop, ko_psth_pop, wt_exc, ko_exc, wt_inh, ko_inh)

%% plot tau
wt_fs_exc_tau = wt_psth_pop.exc_tau(wt_exc);
ko_fs_exc_tau = ko_psth_pop.exc_tau(ko_exc);

wt_fs_inh_tau = wt_psth_pop.inh_tau(wt_inh);
ko_fs_inh_tau = ko_psth_pop.inh_tau(ko_inh);

ecdf_bar_plot(wt_fs_exc_tau, ko_fs_exc_tau)

ecdf_bar_plot([wt_fs_exc_tau, wt_fs_inh_tau],...
    [ko_fs_exc_tau,ko_fs_inh_tau] )

ecdf_bar_plot(wt_fs_inh_tau,ko_fs_inh_tau)

figure(100)
xlabel('Tau (ms)')
xlim([0,100])
title('Putative FS Neurons')
% export_fig('Ecdf_Exc_FS_TimeScale', '-png')

figure(101)
ylabel('Tau (ms)')
ylim([0,50])
title('Putative FS Neurons')
% export_fig('Bar_Exc_TimeScale', '-png')

%%
wt_fs_exc_timescale = wt_psth_pop.timescale(wt_exc);
ko_fs_exc_timescale = ko_psth_pop.timescale(ko_exc);

wt_fs_inh_timescale = wt_psth_pop.inh_timescale(wt_inh);
ko_fs_inh_timescale = ko_psth_pop.inh_timescale(ko_inh);

ecdf_bar_plot(wt_fs_exc_timescale, ko_fs_exc_timescale)

% ecdf_bar_plot([wt_fs_exc_timescale, wt_fs_inh_timescale],...
%     [ko_fs_exc_timescale,ko_fs_inh_timescale] )

ecdf_bar_plot(wt_fs_inh_timescale,ko_fs_inh_timescale)

figure(100)
xlabel('Time scale (ms)')
xlim([0,60])
title('Putative FS Neurons')
% export_fig('Ecdf_Exc_FS_TimeScale', '-png')

figure(101)
ylabel('Time scale (ms)')
ylim([0,50])
title('Putative FS Neurons')
% export_fig('Bar_Exc_TimeScale', '-png')

%% plot the results Lags of FS neurons

wt_fs_exc_lag = wt_psth_pop.lags(wt_exc);
ko_fs_exc_lag = ko_psth_pop.lags(ko_exc);

wt_fs_inh_lag = wt_psth_pop.inh_lags(wt_inh);
ko_fs_inh_lag = ko_psth_pop.inh_lags(ko_inh);

ecdf_bar_plot(wt_fs_exc_lag, ko_fs_exc_lag)
ecdf_bar_plot([wt_fs_exc_lag, wt_fs_inh_lag],...
    [ko_fs_exc_lag,ko_fs_inh_lag] )
ecdf_bar_plot( wt_fs_inh_lag, ko_fs_inh_lag)

figure(100)
xlabel('Lags (ms)')
xlim([0,150])
title('Putative FS Neurons')
% export_fig('Ecdf_FS_lags', '-png')

figure(101)
ylabel('Lags (ms)')
ylim([0,100])
title('Putative FS Neurons')
% export_fig('Bar_Inh_lags', '-png')

%%
%% plot the results Onset of FS neurons

wt_fs_exc_Onset = wt_psth_pop.exc_cp.firstOnset(wt_exc);
ko_fs_exc_Onset = ko_psth_pop.exc_cp.firstOnset(ko_exc);

wt_fs_inh_Onset = wt_psth_pop.inh_cp.firstOnset(wt_inh);
ko_fs_inh_Onset = ko_psth_pop.inh_cp.firstOnset(ko_inh);

ecdf_bar_plot(wt_fs_exc_Onset, ko_fs_exc_Onset)
ecdf_bar_plot([wt_fs_exc_Onset, wt_fs_inh_Onset],...
    [ko_fs_exc_Onset,ko_fs_inh_Onset] )

ecdf_bar_plot(wt_fs_inh_Onset, ko_fs_inh_Onset)

figure(100)
xlabel('Response Onset (s)')
xlim([0,0.1])
title('Putative FS Neurons')
% export_fig('Ecdf_FS_Onset', '-png')

figure(101)
ylabel('Response Onset (s)')
ylim([0,0.1])
title('Putative FS Neurons')
% export_fig('Bar_Inh_lags', '-png')
%% plot the results Offset of FS neurons


wt_fs_exc_Offset = wt_psth_pop.exc_cp.firstOffset(wt_exc);
ko_fs_exc_Offset = ko_psth_pop.exc_cp.firstOffset(ko_exc);

% wt_fs_exc_duration = wt_psth_pop.exc_cp.duration(wt_exc);
% ko_fs_exc_duration = ko_psth_pop.exc_cp.duration(ko_exc);

wt_fs_exc_duration = wt_fs_exc_Offset-wt_fs_exc_Onset;
ko_fs_exc_duration = ko_fs_exc_Offset-ko_fs_exc_Onset;


wt_fs_inh_Offset = wt_psth_pop.inh_cp.firstOffset(wt_inh);
ko_fs_inh_Offset = ko_psth_pop.inh_cp.firstOffset(ko_inh);

% wt_fs_inh_duration = wt_psth_pop.inh_cp.duration(wt_inh);
% ko_fs_inh_duration = ko_psth_pop.inh_cp.duration(ko_inh);

wt_fs_inh_duration = wt_fs_inh_Offset-wt_fs_inh_Onset;
ko_fs_inh_duration = ko_fs_inh_Offset-ko_fs_inh_Onset;



ecdf_bar_plot(wt_fs_exc_duration, ko_fs_exc_duration)
ecdf_bar_plot([wt_fs_exc_duration, wt_fs_inh_duration],...
    [ko_fs_exc_duration,ko_fs_inh_duration] )
ecdf_bar_plot(wt_fs_inh_duration, ko_fs_inh_duration)

figure(100)
xlabel('Response Duration (s)')
xlim([0,0.3])
title('Putative FS Neurons')
% export_fig('Ecdf_excFS_duration', '-png')

figure(101)
ylabel('Response Duration (s)')
ylim([0,0.1])
title('Putative FS Neurons')
% export_fig('Bar_excFS_Duration', '-png')
