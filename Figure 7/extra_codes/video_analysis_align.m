function video_analysis_align(filename)
load(filename)
pupil_diameter = Tosca.Pupil.pupil_size;
% pupil_diameter = Tosca.Pupil.pupil_diametere;

% m_dia = mean(pupil_diameter);
% std_dia = std(pupil_diameter);
% lowBound = m_dia - 5 * std_dia;
% highBound = m_dia + 5 * std_dia;
% figure; 
% histogram(pupil_diameter)
% hold on
% plot([lowBound, lowBound], [0, 100], '-r')
% plot([highBound, highBound], [0, 100], '-r')
% pupil_likelihood = Tosca.Pupil.likelihood;
% idx = find(pupil_likelihood < 0.7);
% pupil_diameter_new = pupil_diameter;
% pupil_diameter_new(idx) = NaN;
% pupil_diameter_new=inpaint_nans(pupil_diameter_new);
% a = diff(pupil_diameter_new);
% indx = find(a>5 | a<-5);
% indx_diff = diff(indx); % interpolate a bigger window
% step = find(indx_diff>10);
% indx_nan = [];
% for i = 1:length(step)
%     if i == 1;
%         indx_nan = [indx_nan, indx(1):indx(step(i))];
%     else
%         indx_nan = [indx_nan, indx(step(i-1)+1):indx(step(i))];
%     end
%     
% end
% pupil_diameter_new(indx_nan+1)=NaN;
% pupil_diameter_new=inpaint_nans(pupil_diameter_new);
indx = find(isnan(pupil_diameter));
seq = [];
for i = 1:length(indx)
    seq(i,:) = (indx(i)-10):(indx(i)+10); % here 10 frame around the NaN were set to NaN
    
end

seq = unique(seq(:));

seq(find(seq<=0)) =[];
seq(find(seq>length(pupil_diameter))) =[];
pupil_diameter_new = pupil_diameter;
pupil_diameter_new(seq) = NaN;
a = diff(pupil_diameter_new);
idx = find(a>120 | a<-120);
indx_diff = diff(idx); % interpolate a bigger window
step = find(indx_diff>5);
indx_nan = [];
for i = 1:length(step)
    if i == 1;
        indx_nan = [indx_nan, idx(1):idx(step(i))];
    else
        indx_nan = [indx_nan, idx(step(i-1)+1):idx(step(i))];
    end
end
    
pupil_diameter_new(indx_nan) = NaN;
pupil_diameter_new=inpaint_nans(pupil_diameter_new, 5);

%% plotting the data

figure;
for i = 1: length(Tosca.Video.eventFrame)
    eventIdx(i) = find(Tosca.Video.Frame_real == Tosca.Video.eventFrame(i));
end
plot(pupil_diameter)
hold on
plot(pupil_diameter_new)
hold on
for i = 1:length(eventIdx)
    plot([eventIdx(i),eventIdx(i)], [600, 1000], 'r')
end
% add the time information
n_timepoint = length(pupil_diameter_new);
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
        if eventIdx(indx(i))+300 > length(pupil_diameter_new)
        else
            video_summary.(['level', num2str(level(j))])(i-iter,:) = pupil_diameter_new(eventIdx(indx(i))-60 : eventIdx(indx(i))+300 );
        end
    end
end
end

for i = 3:(length(eventIdx)-1)
    video_summary.sponta(i-2,:) = pupil_diameter_new((eventIdx(i)-270) : (eventIdx(i)-1) );
end

filename = Tosca.filename;
Tosca.video_summary = video_summary;
Tosca.pupil_area_correct = pupil_diameter_new;
% add time information and stimuli information
Tosca.time = time;
Tosca.event_time = event_time;
Tosca.level = Tosca.StimulusData.Level;
%%
baseline_bin = 46:60; % the event is at 61th frame
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

cd('E:\Ke_Chen\Processed Data\PTCHD1-Project\Video_analysis\Pupil_analysis\Summary')
save(filename, 'Tosca')