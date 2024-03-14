function varargout = spike_sorting(varargin)
% SPIKE_SORTING MATLAB code for spike_sorting.fig
%      SPIKE_SORTING, by itself, creates a new SPIKE_SORTING or raises the existing
%      singleton*.
%
%      H = SPIKE_SORTING returns the handle to a new SPIKE_SORTING or the handle to
%      the existing singleton*.
%
%      SPIKE_SORTING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPIKE_SORTING.M with the given input arguments.
%
%      SPIKE_SORTING('Property','Value',...) creates a new SPIKE_SORTING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before spike_sorting_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to spike_sorting_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help spike_sorting

% Last Modified by GUIDE v2.5 18-Oct-2020 22:14:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @spike_sorting_OpeningFcn, ...
                   'gui_OutputFcn',  @spike_sorting_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before spike_sorting is made visible.
function spike_sorting_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to spike_sorting (see VARARGIN)

% Choose default command line output for spike_sorting
handles.output = hObject;
chanMapPath = 'E:\Ke_Chen\MATLAB\chanMap32x1.mat';
load(chanMapPath)
set(handles.PenNum_tg, 'String', 2);
set(handles.block_num, 'String', '1:10');
set(handles.digitizer_1, 'Value', 1)
set(handles.digitizer_2, 'Value', 0)
set(handles.channel_map_edit, 'String', chanMapPath)
set(handles.chan_num, 'String', length(chanMap))
set(handles.region_tg, 'String', 'ACtx')

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes spike_sorting wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = spike_sorting_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function data_directory_tg_Callback(hObject, eventdata, handles)
% hObject    handle to data_directory_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of data_directory_tg as text
%        str2double(get(hObject,'String')) returns contents of data_directory_tg as a double


% --- Executes during object creation, after setting all properties.
function data_directory_tg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to data_directory_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in loading_tg.
function loading_tg_Callback(hObject, eventdata, handles)
% hObject    handle to loading_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path = uigetdir('Z:\KeChen\Rach Recordings\','Select the data tank folder');
handles.path = path;
PenNum = handles.PenNum_tg.String;
% set the parameter for extracting the waveforms
% eval(handles.block_num.String)
tankName = path;
PenNum = str2num(PenNum);
blockNums = eval(handles.block_num.String);
chans_Num = str2num(get(handles.chan_num, 'String'));
if handles.digitizer_1.Value
    whichdigitizer = 'RSn1'
    chans = 1:chans_Num;
else
    whichdigitizer = 'RSn2'
    chans = 1:chans_Num;
end
% extract waveform 
% set the raw waveform path
region = handles.region_tg.String;
raw_waveform_path = [path,'\', region, '-Raw-Waveform','-', num2str(PenNum)];
set(handles.data_directory_tg, 'HorizontalAlignment', 'left')
messages = {path,'Save raw waveform to...', raw_waveform_path};
set(handles.data_directory_tg, 'String', messages);
direc = raw_waveform_path;

extractWaveforms_Rach(tankName,PenNum,blockNums,chans,direc,whichdigitizer)

new_messages = {messages{:}, 'Finish extracting the waveforms', 'Start sorting the spikes'};
set(handles.data_directory_tg, 'String', new_messages);

% sorting with kilosort
manual = handles.manual_kilosort.Value;
rawDataPath = path;
rawDataPrefix = split(path, '\');
rawDataPrefix = [char(rawDataPrefix(end)), '-', num2str(PenNum)];
chanMapPath = handles.channel_map_edit.String;

NChan = chans_Num;
region = handles.region_tg.String;
process_pipeline_Rach(rawDataPath, rawDataPrefix, PenNum, chanMapPath, NChan, chans, region,blockNums, raw_waveform_path, manual)

% Update handles structure
guidata(hObject, handles);


function PenNum_tg_Callback(hObject, eventdata, handles)
% hObject    handle to PenNum_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PenNum_tg as text
%        str2double(get(hObject,'String')) returns contents of PenNum_tg as a double
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function PenNum_tg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PenNum_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in digitizer_1.
function digitizer_1_Callback(hObject, eventdata, handles)
% hObject    handle to digitizer_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of digitizer_1


% --- Executes on button press in digitizer_2.
function digitizer_2_Callback(hObject, eventdata, handles)
% hObject    handle to digitizer_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of digitizer_2



function block_num_Callback(hObject, eventdata, handles)
% hObject    handle to block_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of block_num as text
%        str2double(get(hObject,'String')) returns contents of block_num as a double
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function block_num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to block_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function chan_num_Callback(hObject, eventdata, handles)
% hObject    handle to chan_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of chan_num as text
%        str2double(get(hObject,'String')) returns contents of chan_num as a double
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function chan_num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chan_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function channel_map_edit_Callback(hObject, eventdata, handles)
% hObject    handle to channel_map_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of channel_map_edit as text
%        str2double(get(hObject,'String')) returns contents of channel_map_edit as a double


% --- Executes during object creation, after setting all properties.
function channel_map_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to channel_map_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chan_map.
function chan_map_Callback(hObject, eventdata, handles)
% hObject    handle to chan_map (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile('E:\Ke_Chen\MATLAB\*.mat','Select Channel Map File');
chanMapPath = [pathname, filename];
handles.channel_map_edit.String = chanMapPath;
set(handles.channel_map_edit, 'HorizontalAlignment', 'left')
set(handles.channel_map_edit, 'String', chanMapPath);
load(chanMapPath)
set(handles.chan_num, 'String', length(chanMap))
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in sorting_tg.
function sorting_tg_Callback(hObject, eventdata, handles)
% hObject    handle to sorting_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% offline sorting with Phy
%%% Phy2
%%%% Anaconda prompt
%%%% cd /d E:\Ke_Chen\RawData\Rach_Ephys
%%%% Activate the phy enviroment: conda activate phy2
%%%% conda activate phy2
%%% Phy template-gui params.py



function region_tg_Callback(hObject, eventdata, handles)
% hObject    handle to region_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of region_tg as text
%        str2double(get(hObject,'String')) returns contents of region_tg as a double


% --- Executes during object creation, after setting all properties.
function region_tg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to region_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in manual_kilosort.
function manual_kilosort_Callback(hObject, eventdata, handles)
% hObject    handle to manual_kilosort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of manual_kilosort


% --- Executes on button press in assemble_tg.
function assemble_tg_Callback(hObject, eventdata, handles)
% hObject    handle to assemble_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path = uigetdir('Z:\KeChen\','Select the folder with the sorted spikes');
analdir = path;
chanMapPath = handles.channel_map_edit.String;
kilosort_spiketimes_to_clusterData(analdir,chanMapPath)


% --- Executes on button press in Extract_waveforms_tosca.
function Extract_waveforms_tosca_Callback(hObject, eventdata, handles)
% hObject    handle to Extract_waveforms_tosca (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path = uigetdir('Z:\KeChen\Tosca Recording\','Select the data tank folder');
handles.path = path;
% set the parameter for extracting the waveforms
% eval(handles.block_num.String)
tankName = path;
blockNums = eval(handles.block_num.String);
chans_Num = str2num(get(handles.chan_num, 'String'));
if handles.digitizer_1.Value
    whichdigitizer = 'RSn1'
    chans = 1:chans_Num;
else
    whichdigitizer = 'RSn2'
    chans = 1:chans_Num;
end
% extract waveform 
% set the raw waveform path
region = handles.region_tg.String;
% path_seg = split(path, '\');
% path_tosave = char(strjoin(path_seg(1:end-1), '\'));
% raw_waveform_path = [path_tosave,'\', region, '-Raw-Waveform'];
raw_waveform_path = [path,'\', region, '-Raw-Waveform'];
set(handles.data_directory_tg, 'HorizontalAlignment', 'left')
messages = {path,'Save raw waveform to...', raw_waveform_path};
set(handles.data_directory_tg, 'String', messages);
direc = raw_waveform_path;

% extract the tosca data
addpath(genpath('E:\Ke_Chen\MATLAB\Tosca\MATLAB'))
% tosca_path = [char(strjoin(path_seg(1:end-2), '\')), '\Tosca Files'];
tosca_path = [path, '\Tosca Files'];

cd(tosca_path)
file = dir('*.trace.txt');
file_parts = split(file.name, '-');
file_lastparts = split(file_parts(end), '.');
file_run = strjoin({char(file_parts(1)), char(file_parts(2)), char(file_lastparts(1))}, '-');
[Data, Params] = tosca_read_run([file_run,'.txt']);

for iTrial=1:length(Data)
    Tosca.StimulusData.Level(iTrial) = Data{1, iTrial}.Sound.Freefield.Level.dB_SPL;
    Tosca.trial(iTrial) =tosca_read_trial(Params, Data, iTrial);
%     idxStateChange=find((trial.State_Change(2:end)-trial.State_Change(1:end-1))==1)+1;
%     tOnset(iTrial)=trial.Time_s(idxStateChange(1))-trial.Time_s(1);
%     t(iTrial) = trial.Time_s(end)-trial.Time_s(1);
end
% extract all data
data_path = [path, '\', file_run];
data = TDTbin2mat(data_path);
Tosca.event = data.epocs.ToTr;
waveforms = data.streams.(whichdigitizer).data * 1e6;
waveforms = int16(waveforms);
mkdir(raw_waveform_path)
save([raw_waveform_path, '\tosca_event.mat'], 'Tosca')
fid = fopen([raw_waveform_path, '\data.bin'], 'w');
fwrite(fid, waveforms, 'int16');
fclose(fid)
cd(raw_waveform_path)
applyCARtoDat('data.bin', chans_Num);
kilosort


% --- Executes on button press in tosca_tg.
function tosca_tg_Callback(hObject, eventdata, handles)
% hObject    handle to tosca_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% tosca_path = uigetdir('F:\KeChen\RawData\Tosca_Ephys\','Select the data tank folder');
% handles.tosca_path = tosca_path;
% addpath(genpath('E:\Ke_Chen\MATLAB\Tosca\MATLAB'));
% cd(tosca_path)
% file = dir('*.trace.txt');
% file_parts = split(file.name, '-');
% file_lastparts = split(file_parts(end), '.');
% file_run = strjoin({char(file_parts(1)), char(file_parts(2)), char(file_lastparts(1))}, '-');
% [Data, Params] = tosca_read_run([file_run,'.txt']);
% 
% for iTrial=1:length(Data)
%     Tosca.StimulusData.Level(i) = Data{1, iTrial}.Sound.Freefield.Level.dB_SPL;
%     Tosca.trial(iTrial) =tosca_read_trial(Params, Data, iTrial);
% %     idxStateChange=find((trial.State_Change(2:end)-trial.State_Change(1:end-1))==1)+1;
% %     tOnset(iTrial)=trial.Time_s(idxStateChange(1))-trial.Time_s(1);
% %     t(iTrial) = trial.Time_s(end)-trial.Time_s(1);
% end


% --- Executes on button press in assemble_tosca.
function assemble_tosca_Callback(hObject, eventdata, handles)
% hObject    handle to assemble_tosca (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path = uigetdir('Z:\KeChen\Tosca Recording\','Select the folder with the sorted spikes');
analdir = path;
chanMapPath = handles.channel_map_edit.String;
kilosort_spiketimes_Tosca_clusterData(analdir)
