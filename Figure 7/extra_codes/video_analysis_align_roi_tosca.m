function video_analysis_align_roi_tosca(filename)
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
h1 = plot(Mtrace_anterior);
hold on
h2 = plot(Mtrace_posterior);
hold on

for i = 1:length(eventIdx)
    h3 = plot([eventIdx(i),eventIdx(i)], [0, 20], 'k');
end
ylim([0,30])
legend([h1, h2, h3], {'Anterior', 'Posterior', 'Sound'})

% add the time information
n_timepoint = length(Mtrace_posterior);
t_interval = 0.033;  % 30 Hz
time = 0:t_interval:((n_timepoint-1) * t_interval);
event_time = time(eventIdx);

%% Fix the bug when oscassionally the trial number recorded from vidoe and tosca does not match
if length(eventIdx) < length(Tosca.StimulusData.Level)
    warning('Less markers in the videos relative to tosca')
    event_interval = diff(Tosca.Video.eventFrame)*0.033;
    outlier_frame = find(event_interval>22); % two events are 22s apart
    if length(outlier_frame) ~=1
        error('Something is wrong')
    else
        Tosca.StimulusData.Level(outlier_frame+1) =[];
    end
elseif length(eventIdx) > length(Tosca.StimulusData.Level)
    error('Something is wrong')
    
end
%% trying to store the 
level = unique(Tosca.StimulusData.Level);
for i = 1: length(Tosca.Video.eventFrame)
    eventIdx(i) = find(Tosca.Video.Frame_real == Tosca.Video.eventFrame(i));
end
for j = 1:length(level)
    indx = find(Tosca.StimulusData.Level ==level(j));
%     hold on
%     for i = 1:length(indx)
%         plot([eventIdx(indx(i)),eventIdx(indx(i))], [20, 40], 'r')
%     end
iter = 0;
    for i = 1:length(indx)
        if indx(i) ==1
            iter = 1; % skip the first sound event
        else

            video_summary.Mtrace_anterior.(['level', num2str(level(j))])(i-iter,:) = Mtrace_anterior(eventIdx(indx(i))-60 : eventIdx(indx(i))+300 );
            video_summary.Mtrace_posterior.(['level', num2str(level(j))])(i-iter,:) = Mtrace_posterior(eventIdx(indx(i))-60 : eventIdx(indx(i))+300 );
        
        end
    end
end

for i = 3:(length(eventIdx)-1)
    video_summary.sponta_anter(i-2,:) = Mtrace_anterior((eventIdx(i)-270) : (eventIdx(i)-1) );
    video_summary.sponta_post(i-2,:) = Mtrace_posterior((eventIdx(i)-270) : (eventIdx(i)-1) );
end

% filename = Tosca.filename;
Tosca.video_summary = video_summary;
% add time information and stimuli information
Tosca.time = time;
Tosca.event_time = event_time;
Tosca.level = Tosca.StimulusData.Level;

%% trying to store the

% filename = Tosca.filename;

%%
baseline_bin = 56:60; % the impale event is between 60 and 61th frame,
% f = fieldnames(Tosca.video_summary);
f = {'Mtrace_anterior', 'Mtrace_posterior'};
for i = 1:length(f)
    f_level  = fieldnames(Tosca.video_summary.(f{i}));
    for k = 1:length(f_level)
        temp = Tosca.video_summary.(f{i}).(f_level{k});
        for j = 1:size(temp, 1)
            baseline = temp(j,baseline_bin);
            Tosca.video_summary.(f{i}).([f_level{k}, '_norm'])(j,:) = (temp(j,:) - mean(baseline))/mean(baseline);
        end
        Tosca.video_summary.(f{i}).([f_level{k}, '_avg']) = mean(Tosca.video_summary.(f{i}).([f_level{k}, '_norm']),1);
    end
end
%%
cd('E:\Ke_Chen\Processed Data\PTCHD1-Project\Video_analysis\ROIbased_Tosca\Summary')
Tosca.filename = filename;
save(filename, 'Tosca')