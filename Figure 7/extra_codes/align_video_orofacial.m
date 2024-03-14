
function align_video_orofacial(folder)
%% load the files
%% load the package from extracting Tosca events

addpath(genpath('E:\Ke_Chen\MATLAB\Tosca\MATLAB'));
% folder = pwd;
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

%% load the orofacial quantified by roi ////you need to check the ROI analysis before you run this session
cd(folder)
file = dir('*_Orofacial.mat');
load(file.name)
Tosca.orofacial = Orofacial;
file_parts = split(filename, '.');
file_tosave = char(file_parts(1));
% cd('E:\Ke_Chen\Processed Data\PTCHD1-Project\Video_analysis\ROIbased')
cd('E:\Ke_Chen\Processed Data\PTCHD1-Project\Video_analysis\ROIbased_Tosca')
save([file_tosave, '.mat'], 'Tosca')
% cd('Z:\KeChen\Analysis\Noiseburst')
filename = [file_tosave, '.mat'];
video_analysis_align_roi_tosca(filename)


