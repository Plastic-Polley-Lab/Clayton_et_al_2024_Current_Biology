function varargout = analysis_toolbox(varargin)
% ANALYSIS_TOOLBOX MATLAB code for analysis_toolbox.fig
%      ANALYSIS_TOOLBOX, by itself, creates a new ANALYSIS_TOOLBOX or raises the existing
%      singleton*.
%
%      H = ANALYSIS_TOOLBOX returns the handle to a new ANALYSIS_TOOLBOX or the handle to
%      the existing singleton*.
%
%      ANALYSIS_TOOLBOX('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ANALYSIS_TOOLBOX.M with the given input arguments.
%
%      ANALYSIS_TOOLBOX('Property','Value',...) creates a new ANALYSIS_TOOLBOX or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before analysis_toolbox_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to analysis_toolbox_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help analysis_toolbox

% Last Modified by GUIDE v2.5 15-Apr-2021 11:14:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @analysis_toolbox_OpeningFcn, ...
                   'gui_OutputFcn',  @analysis_toolbox_OutputFcn, ...
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


% --- Executes just before analysis_toolbox is made visible.
function analysis_toolbox_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to analysis_toolbox (see VARARGIN)

% Choose default command line output for analysis_toolbox
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes analysis_toolbox wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = analysis_toolbox_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in drc_analysis_tg.
function drc_analysis_tg_Callback(hObject, eventdata, handles)
% hObject    handle to drc_analysis_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
edit DRC_analysis


% --- Executes on button press in fra_analysis_tg.
function fra_analysis_tg_Callback(hObject, eventdata, handles)
% hObject    handle to fra_analysis_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
edit fra_analysis_singleUnits


% --- Executes on button press in raster_pattern_tg.
function raster_pattern_tg_Callback(hObject, eventdata, handles)
% hObject    handle to raster_pattern_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path = uigetdir('E:\Ke_Chen\Processed Data\Rach Recording\','Select the folder that stores sorted spikes');
raster_data_pattern(path)


% --- Executes on button press in rate_level_tg.
function rate_level_tg_Callback(hObject, eventdata, handles)
% hObject    handle to rate_level_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
edit ratelevel_spike


% --- Executes on button press in pattern_analysis_tg.
function pattern_analysis_tg_Callback(hObject, eventdata, handles)
% hObject    handle to pattern_analysis_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
edit pattern_analysis


% --- Executes on button press in append_pdf_tg.
function append_pdf_tg_Callback(hObject, eventdata, handles)
% hObject    handle to append_pdf_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uigetfile('E:\Ke_Chen\Processed Data\Rach Recording\*.pdf','Select the folder that stores sorted spikes','Select the pdf file');
f_split = split(FileName, '_');
% neuron_id = split(f_split(end), '.');
% neuron_id = str2num(neuron_id{1});
filenames = strjoin(f_split(1:end-1), '_');
files = struct2cell(dir([PathName, '\', char(filenames), '*.pdf']))
file_list = files(1,:)';
addpath('E:\Ke_Chen\MATLAB\Export_fig')
cd(PathName)
append_pdfs([char(strjoin(f_split(1:end-1), '_')), '_summary.pdf'], file_list)


% --- Executes on button press in fra_impale_tg.
function fra_impale_tg_Callback(hObject, eventdata, handles)
% hObject    handle to fra_impale_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fra_impale_allChan


% --- Executes on button press in noiseburst_tg.
function noiseburst_tg_Callback(hObject, eventdata, handles)
% hObject    handle to noiseburst_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
edit noiseburst_analysis


% --- Executes on button press in paired_pulse_tg.
function paired_pulse_tg_Callback(hObject, eventdata, handles)
% hObject    handle to paired_pulse_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
edit paired_pulse_analysis


% --- Executes on button press in raster_tosca_tg.
function raster_tosca_tg_Callback(hObject, eventdata, handles)
% hObject    handle to raster_tosca_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Tosca_raster


% --- Executes on button press in rlf_tosca.
function rlf_tosca_Callback(hObject, eventdata, handles)
% hObject    handle to rlf_tosca (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
edit rlf_tosca


% --- Executes on button press in allen_cc_tg.
function allen_cc_tg_Callback(hObject, eventdata, handles)
% hObject    handle to allen_cc_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
addpath(genpath('E:\Ke_Chen\MATLAB\allenCCF'))
allen_ccf_npx


% --- Executes on button press in video_analysis_tg.
function video_analysis_tg_Callback(hObject, eventdata, handles)
% hObject    handle to video_analysis_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
edit video_analysis


% --- Executes on button press in preprocess_tg.
function preprocess_tg_Callback(hObject, eventdata, handles)
% hObject    handle to preprocess_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

path = uigetdir('Z:\KeChen\Analysis', 'Select the folder where the DLC analyzed results store')
tic
DLC_Video_PreProcess(path)
toc


% --- Executes on button press in abr_analysis_tg.
function abr_analysis_tg_Callback(hObject, eventdata, handles)
% hObject    handle to abr_analysis_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
edit ABR_analysis


% --- Executes on button press in get_roi_tg.
function get_roi_tg_Callback(hObject, eventdata, handles)
% hObject    handle to get_roi_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
edit get_orofacial_ROIs


% --- Executes on button press in create_video_tg.
function create_video_tg_Callback(hObject, eventdata, handles)
% hObject    handle to create_video_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
edit Create_Videos


% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
edit signal_noise_correlation


% --- Executes on button press in manually_tg.
function manually_tg_Callback(hObject, eventdata, handles)
% hObject    handle to manually_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uigetfile('*.mat','Select the rastered spike files')
path = [PathName,FileName];
manually_correct_spikes_tosca(path)


% --- Executes on button press in extract_orofacial_tg.
function extract_orofacial_tg_Callback(hObject, eventdata, handles)
% hObject    handle to extract_orofacial_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path = uigetdir('Z:\KeChen\Analysis', 'Select the folder where the video and orofacial ROIs are')
tic
fprintf('Processing Video # %d of %d\n', 1, 1)
video_analysis_basedROI(path)
toc



% --- Executes on button press in align_tg.
function align_tg_Callback(hObject, eventdata, handles)
% hObject    handle to align_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
folder = uigetdir('Z:\KeChen\Analysis', 'Select the folder where the analyzed orofacial are')

align_video_orofacial(folder)


% --- Executes on button press in pupil_align_tg.
function pupil_align_tg_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_align_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
folder = uigetdir('Z:\KeChen\Analysis', 'Select the folder where the analyzed DLC results are')
DLC_Video_Align(folder)
