function video_analysis_align_roi(filename)
load(filename)
Mtrace_anterior = Tosca.orofacial.Mtrace_anterior;
Mtrace_posterior = Tosca.orofacial.Mtrace_posterior;



%% plotting the data
% idx = find(Tosca.Video.frame<0);
% event_time = Tosca.Video.time(idx) + 0.1 ;% there is 100 ms delay
% for i = 1:length(event_time)
%     indx_event(i)= min(find(Tosca.Video.time>event_time(i)));
%     
% end

figure;
for i = 1: length(Tosca.Video.eventFrame)
    eventIdx(i) = find(Tosca.Video.Frame_real == Tosca.Video.eventFrame(i));
end
plot(Mtrace_anterior)
hold on
plot(Mtrace_posterior)
hold on
for i = 1:length(eventIdx)
    plot([eventIdx(i),eventIdx(i)], [0, 20], 'r')
end
ylim([0,30])
%% trying to store the 

for i = 1:length(eventIdx)-1
    video_summary.Mtrace_anterior(i,:) = Mtrace_anterior(eventIdx(i)-10 : eventIdx(i)+16);
    video_summary.Mtrace_Posterior(i,:) = Mtrace_posterior(eventIdx(i)-10 : eventIdx(i)+16);
end
t = -0.033 * 10 :0.033:0;
t_post = 0.033:0.033:0.033*16;
t = [t, t_post];

% filename = Tosca.filename;
Tosca.video_summary = video_summary;

%%
baseline_bin = 9:13; % the impale event is between 10 and 11th frame, there was 100 ms delay
f = fieldnames(Tosca.video_summary);
for i = 1:length(f)
    temp = Tosca.video_summary.(f{i});
    for j = 1:size(temp, 1)
        baseline = temp(j,baseline_bin);
        Tosca.video_summary.([f{i}, '_norm'])(j,:) = (temp(j,:) - mean(baseline))/mean(baseline);
    end
    Tosca.video_summary.([f{i}, '_avg']) = mean(Tosca.video_summary.([f{i}, '_norm']),1);
end
%%
cd('E:\Ke_Chen\Processed Data\PTCHD1-Project\Video_analysis\ROIbased\Summary')
save(filename, 'Tosca')