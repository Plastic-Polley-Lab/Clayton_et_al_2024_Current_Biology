%% load the files
function DLC_Video_Align(folder)
addpath(genpath('E:\Ke_Chen\MATLAB\Tosca\MATLAB'));
% folder = pwd;
save_path = 'E:\Ke_Chen\Processed Data\PTCHD1-Project\Video_analysis\Pupil_analysis'; % where you save the data
%%
% Step 1: load the tosca files
cd([folder, '\', 'Tosca Files'])
% tosca_path = uigetdir('F:\KeChen\RawData\Tosca_Ephys\','Select the data tank folder');
% cd(tosca_path)
file = dir('*.trace.txt');
file_parts = split(file.name, '-');
file_lastparts = split(file_parts(end), '.');
file_run = strjoin({char(file_parts(1)), char(file_parts(2)), char(file_lastparts(1))}, '-');
[Data, Params] = tosca_read_run([file_run,'.txt']);

for iTrial=1:length(Data)
    Tosca.StimulusData.Level(iTrial) = Data{1, iTrial}.Sound.Freefield.Level.dB_SPL;
%     Tosca.trial(iTrial) =tosca_read_trial(Params, Data, iTrial);
%     idxStateChange=find((trial.State_Change(2:end)-trial.State_Change(1:end-1))==1)+1;
%     tOnset(iTrial)=trial.Time_s(idxStateChange(1))-trial.Time_s(1);
%     t(iTrial) = trial.Time_s(end)-trial.Time_s(1);
end
%%
% Step 2: load the trigger event in the videos files
% clear
cd([folder, '\', 'Video Events'])
file_videotext = dir('*.avi.txt');
time =[];
diam =[];
frame =[];
for i = 1:length(file_videotext)
    filename=file_videotext(i).name;
    table=table2array(readtable(filename));
    time=[time table(:,3)'];
    frame=[frame table(:,2)'];
end
Tosca.Video.time = time;
Tosca.Video.frame = frame;
Tosca.Video.event = find(Tosca.Video.frame <0) + 1;
Tosca.Video.eventFrame = frame(Tosca.Video.event);
frame(Tosca.Video.event - 1) =[];
Tosca.Video.Frame_real = frame;
find(diff(diff(Tosca.Video.Frame_real))>0) % sanity check
find(diff(diff(Tosca.Video.Frame_real))<0) % sanity check

%%
% Step 3: load the pupil size
cd(folder)
load('Summary_data.mat')
pupil_size =[];
pupil_diameter =[];
pupil_error =[];
pupil_raw =[];
snout = [];
for i = 1: length(data)
    pupil_size = [pupil_size, data(i).pupil_fit.pupil_area];
    pupil_diameter = [pupil_diameter; data(i).distance_ap'];
    pupil_raw = [pupil_raw, data.pupil];
    snout     = [snout, data(i).snout];
end
Tosca.Pupil.pupil_size     = pupil_size;
Tosca.Pupil.pupil_diameter = pupil_diameter;
Tosca.Pupil.pupil_raw      = pupil_raw;
Tosca.Snout                = snout;

%% Step 4: save the data
cd(save_path)
filename = [file_run, '.mat'];
% filename = 'KeC24-Session1.mat';
Tosca.filename = filename;
save(filename, 'Tosca')

%% Step 5:Let's remove artifacts and summarize the pupil size changes
video_analysis_align(filename)

