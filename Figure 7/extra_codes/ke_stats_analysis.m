function varargout = ke_stats_analysis(varargin)
% KE_STATS_ANALYSIS MATLAB code for ke_stats_analysis.fig
%      KE_STATS_ANALYSIS, by itself, creates a new KE_STATS_ANALYSIS or raises the existing
%      singleton*.
%
%      H = KE_STATS_ANALYSIS returns the handle to a new KE_STATS_ANALYSIS or the handle to
%      the existing singleton*.
%
%      KE_STATS_ANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in KE_STATS_ANALYSIS.M with the given input arguments.
%
%      KE_STATS_ANALYSIS('Property','Value',...) creates a new KE_STATS_ANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ke_stats_analysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ke_stats_analysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ke_stats_analysis

% Last Modified by GUIDE v2.5 27-Mar-2021 10:50:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ke_stats_analysis_OpeningFcn, ...
                   'gui_OutputFcn',  @ke_stats_analysis_OutputFcn, ...
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


% --- Executes just before ke_stats_analysis is made visible.
function ke_stats_analysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ke_stats_analysis (see VARARGIN)

% Choose default command line output for ke_stats_analysis
handles.output = hObject;
set(handles.window_tg, 'String', '1:1:1000');
set(handles.n_repeats, 'String', 10);
set(handles.n_kfold, 'String', 10);
set(handles.set_num_tg, 'String', 1)
set(handles.cycle_length_tg, 'String', 4)
set(handles.output_directory_tg,'String','E:\Ke_Chen\Behavior\Patterns\test')
set(handles.jitter_tg, 'String', 5)
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ke_stats_analysis wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ke_stats_analysis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%% Statistical Test %%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in paired_t_test.
function paired_t_test_Callback(hObject, eventdata, handles)
% hObject    handle to paired_t_test (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uigetfile({'*.xlsx'},'Select the Excel file');
[num, ~, raw] = xlsread([PathName,FileName]);
Pairline_plot(num)
% tab = cell2table(raw(2:end,:),'VariableNames', raw(1,:))
% figure
% % re-organize the data
% for i = 1:size(num,2)
%     data{i} = num(~isnan(num(:,i)),1)
% end
% plotSpread(data,'distributionColors',{'r','b'})
xticks([1,2])
xticklabels(raw(1,:))
% ylabel(input('Label the Y axis\n'))


% --- Executes on button press in t_test.
function t_test_Callback(hObject, eventdata, handles)
% hObject    handle to t_test (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uigetfile({'*.xlsx'},'Select the Excel file');
[num, ~, raw] = xlsread([PathName,FileName]);
raw_table = readtable([PathName,FileName]);
% raw_stack = stack(raw,1:size(raw, 2), 'IndexVariableName','Condition', 'NewDataVariableName','Measures');
% measures = table2array(raw_stack(:, {'Measures'}));
% condition =  cellstr(string(table2cell(raw_stack(:, {'Condition'}))));
% figure;
% violinplot(raw_table)
bar_plot(num)
xticks([1,2])
xticklabels(raw(1,:))
handles.ttest_num = num;
handles.ttest_raw = raw;
handles.ttest_rawTable = raw_table;
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in anova_one.
function anova_one_Callback(hObject, eventdata, handles)
% hObject    handle to anova_one (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uigetfile({'*.xlsx'},'Select the Excel file');
raw_table = readtable([PathName,FileName]);
figure;
violin_plot(raw_table)
% print(raw_table)

% --- Executes on button press in anova_two.
function anova_two_Callback(hObject, eventdata, handles)
% hObject    handle to anova_two (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in chi_test.
function chi_test_Callback(hObject, eventdata, handles)
% hObject    handle to chi_test (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in plot_style.
function plot_style_Callback(hObject, eventdata, handles)
% hObject    handle to plot_style (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plot_style contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plot_style
contents = cellstr(get(hObject,'String'));
style = contents{get(hObject,'Value')};
handles.style = style;
if isfield(handles, 'ttest_num')
    num = handles.ttest_num;
    raw = handles.ttest_raw;
    raw_table = handles.ttest_rawTable;
end
switch style
    case 'Bar'
        bar_plot(num)
        xticks(1:length(raw(1,:)))
        xticklabels(raw(1,:))
    case 'Box'
    case 'Violin'
        violin_plot(raw_table)
end

% Update handles structure
% guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function plot_style_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_style (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%% Analysis for Impale %%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in fra_impale.
function fra_impale_Callback(hObject, eventdata, handles)
% hObject    handle to fra_impale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uigetfile('F:\KeChen\RawData\Rach_Ephys\*.mat','Select the FRA Impale file')
path = [PathName,FileName];
chan = 1;
[fra, freqs, spls] = fra_impale(path);
figure(100)
fra_plot(fra, freqs, spls, chan)
handles.path = path;
handles.fra = fra;
handles.freqs = freqs;
handles.spls  = spls;
% Update handles structure
guidata(hObject, handles);


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
contents = cellstr(get(hObject,'String'));
chan = str2num(contents{get(hObject,'Value')});
figure(100)
fra = handles.fra;
freqs = handles.freqs;
spls = handles.spls;
fra_plot(fra, freqs, spls, chan)
figure(handles.figure1)




% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes on button press in levelrate.
function levelrate_Callback(hObject, eventdata, handles)
% hObject    handle to levelrate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uigetfile('F:\KeChen\RawData\Rach_Ephys\*.mat','Select the RateLevel Impale file')
path = [PathName,FileName];
[rlf, spls] = rlf_impale(path);
handles.rlf = rlf;
handles.rlf_spls = spls;
figure(101)
chan = 1;
rlf_plot(rlf, spls, chan)
% Update handles structure
guidata(hObject, handles);


% --- Executes on selection change in levelRate_channel.
function levelRate_channel_Callback(hObject, eventdata, handles)
% hObject    handle to levelRate_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns levelRate_channel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from levelRate_channel
contents = cellstr(get(hObject,'String'));
chan = str2num(contents{get(hObject,'Value')});
figure(101)
rlf = handles.rlf;
spls = handles.rlf_spls;
rlf_plot(rlf, spls, chan)
figure(handles.figure1)


% --- Executes during object creation, after setting all properties.
function levelRate_channel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to levelRate_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on anova_two and none of its controls.
function anova_two_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to anova_two (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on levelRate_channel and none of its controls.
function levelRate_channel_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to levelRate_channel (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
% get(gcf,'currentkey')
 % determine the key that was pressed 

function popupmenu1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to anova_two (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)         
         


% --- Executes on selection change in chan_both.
function chan_both_Callback(hObject, eventdata, handles)
% hObject    handle to chan_both (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns chan_both contents as cell array
%        contents{get(hObject,'Value')} returns selected item from chan_both
contents = cellstr(get(hObject,'String'));
chan = str2num(contents{get(hObject,'Value')});
figure(100)
fra = handles.fra;
freqs = handles.freqs;
spls = handles.spls;
fra_plot(fra, freqs, spls, chan)
figure(101)
rlf = handles.rlf;
spls = handles.rlf_spls;
rlf_plot(rlf, spls, chan)
figure(handles.figure1)

% --- Executes during object creation, after setting all properties.
function chan_both_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chan_both (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in psth.
function psth_Callback(hObject, eventdata, handles)
% hObject    handle to psth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uigetfile('F:\KeChen\RawData\Rach_Ephys\*.mat','Select the RateLevel Impale file')
path = [PathName,FileName];
window = eval(handles.window_tg.String);
handles.window = window;
psth = psth_impale(path, [0,window]);
handles.psth = psth;
figure(102)
chan = 1;
psth_plot(psth, chan, window,1, chan)
% Update handles structure
guidata(hObject, handles);


% --- Executes on selection change in psth_channel.
function psth_channel_Callback(hObject, eventdata, handles)
% hObject    handle to psth_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns psth_channel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from psth_channel
contents = cellstr(get(hObject,'String'));
chan = str2num(contents{get(hObject,'Value')});
figure(102)
psth = handles.psth;
window = eval(handles.window_tg.String);
psth_plot(psth, chan, window, 0, chan)
figure(handles.figure1)

% --- Executes during object creation, after setting all properties.
function psth_channel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to psth_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%% Frequency Domain Analysis %%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in csd.
function csd_Callback(hObject, eventdata, handles)
% hObject    handle to csd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotCSD_rach_KeC32x2


% --- Executes on button press in spike_sorting.
function spike_sorting_Callback(hObject, eventdata, handles)
% hObject    handle to spike_sorting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
spike_sorting

% --- Executes on button press in raster_tg.
function raster_tg_Callback(hObject, eventdata, handles)
% hObject    handle to raster_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path = uigetdir('E:\Ke_Chen\Processed Data\Rach Recording\','Select the folder that stores sorted spikes');
raster_data(path)

% --- Executes on button press in svm_2md.
function svm_2md_Callback(hObject, eventdata, handles)
% hObject    handle to svm_2md (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uigetfile('F:\KeChen\RawData\Rach_Ephys\*.mat','Select the RateLevel Impale file')
path = [PathName,FileName];
load(path);  % data should be organized into a structure ready for decoding, and being name data
n_kfold = str2double(handles.n_kfold.String);
n_repeats = str2double(handles.n_repeats.String);
[accuracy_decoding,Confusion_all]= svm_decoding(data,n_kfold, n_repeats);
axis square
fprintf('average decoding accuracy is %d\n', mean(accuracy_decoding))
% Compute validation predictions and scores
% confusionchart(testdata.Response, predictions)



function n_kfold_Callback(hObject, eventdata, handles)
% hObject    handle to n_kfold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of n_kfold as text
%        str2double(get(hObject,'String')) returns contents of n_kfold as a double
% handles.n_kfold = get(hObject,'String');
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function n_kfold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n_kfold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function n_repeats_Callback(hObject, eventdata, handles)
% hObject    handle to n_repeats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of n_repeats as text
%        str2double(get(hObject,'String')) returns contents of n_repeats as a double
% handles.n_repeats = get(hObject,'String');
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function n_repeats_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n_repeats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fs = 1000;
t = 0:1/fs:1-1/fs;
x = 1.8*cos(2*pi*100*t);
[pxx,f] = periodogram(x,[],length(x),fs);
figure;
plot(f, pxx, 'LineWidth', 1)
box off
ylabel('Power Spectral Density')
xlabel('Frequency (Hz)')
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,400,300])
hold off


% --- Executes on button press in patterns_tg.
function patterns_tg_Callback(hObject, eventdata, handles)
% hObject    handle to patterns_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set_num = str2num(handles.set_num_tg.String);
CycleLength = str2num(handles.cycle_length_tg.String);
outputFolder = handles.output_directory_tg.String;
jitters = str2num(handles.jitter_tg.String);
file_sound = fullfile(outputFolder, sprintf('final_cyc%d_regrand_set%d_jitter%d.wav',CycleLength, set_num, jitters));
if handles.classic.Value
    [result, output,reg_interval,t] = NoisePattern_ke_gui(set_num, CycleLength, jitters, outputFolder);
    
else
    idx = find([handles.jitter_fixed.Value, handles.jitter2_rand.Value, handles.scale_tg.Value, handles.swap_tg.Value])
    switch idx
        case 1
            jitter_type = 'jitter-fixed';
        case 2
            jitter_type = 'jitter-rand'
        case 3
            jitter_type = 'jitter-scale'
        case 4
            jitter_type = 'jitter-swap'
    end
    if handles.complex_tg.Value
        fprintf('Generate jittered Interval with the same reg_interval')
        %%%%%%%%%%%%%%%% you modify the jitters and intervals here %%%%%
        jitters = [0.1, 0.15, 0.20, 0.25, 0.3, 0.35, 0.4];
        intervals = 250:10:450;        % increase the inter-click interval, so that the jittered version is still above 100 ms
        Variable.CycleLength = CycleLength;
        Variable.NumCycles   = 25;
        Variable.Seed  = -1;
        reg_interval = randomSeq_noRepeats(intervals, Variable);
%         [t_reg, t_rand, t_randJitter, t_scale] = complex_sound(intervals, reg_interval, jitters, set_num, CycleLength, outputFolder);
        [~, ~, t] = complex_sound(intervals, reg_interval, jitters, set_num, CycleLength, outputFolder);
        % let me add some initial rand sequence to account for the
        % adaptation
        t_pre_rand = [];
        t_reg = unique(t(end-1,:));
        for i = 1:5     % random sample 5 frames intervals;
            t_pre_rand = [t_pre_rand, t_reg(randperm(CycleLength))]; % random sample 5 frames intervals;
        end
        t_jitter_asecending  = [t(15,:); t(1:7,:); t(end,:)];
        t_jitter_desecending = flipud(t_jitter_asecending);
        t_scale_asecending   = [t(15,:); t(8:14,:); t(end,:)];
        t_scale_desecending   = flipud(t_scale_asecending);
        answer = questdlg('Would you like to save the sequence and audio file?');
        switch answer
            case 'Yes'
                file = sprintf('final_cyc%d_reg_rand_set%d_%s',CycleLength, set_num, 'jitter');
                outputfile = [outputFolder,'\', file];
%                 savefile_writeaudio(t_sequencies, outputfile) % save t_reg
%                 save([outputfile, '.mat'],'t_sequencies', 't_sequence', 't', 't_pre_rand');

%                 plot_patternSeq_gui(t_jitter_asecending, 4, [0,jitters,1], 'jitter')
                plot_patternSeq_gui(t_jitter_desecending, 4, fliplr([0,jitters,1]), 'jitter')
                savefig([file, '_desecending.fig'])
                plot_patternSeq_gui(t_scale_desecending, 4, fliplr([0,jitters,1]), 'Scale')
                savefig([file(1:end-6),'scale', '_desecending.fig'])
                plot_patternSeq_gui(t_jitter_asecending, 4, [0,jitters,1], 'jitter')
                savefig([file, '_asecending.fig'])
                plot_patternSeq_gui(t_scale_asecending, 4, [0,jitters,1], 'Scale')
                savefig([file(1:end-6),'scale', '_asecending.fig'])
                
                save([outputfile, '.mat'],'t_jitter_asecending', 't_jitter_desecending', 't_scale_asecending', ...
                    't_scale_desecending', 't', 't_pre_rand');

%                 savefile_writeaudio(t_reg, outputfile) % save t_reg
%                 
%                 file = sprintf('final_cyc%d_reg_set%d_%s%s',CycleLength, set_num, 'jitter','-rand-rand');
%                 outputfile = [outputFolder,'\', file];
%                 savefile_writeaudio(t_rand, outputfile) % save t_rand
%                 
%                 for i = 1: size(t_randJitter,1)
%                     file = sprintf('final_cyc%d_reg_set%d_%s%d',CycleLength, set_num, 'jitter-rand', jitters(i)*100);
%                     outputfile = [outputFolder,'\', file];
%                     savefile_writeaudio(t_randJitter(i,:), outputfile) % save t_rand
%                 end
%                 
%                 for i = 1: size(t_scale,1)
%                     file = sprintf('final_cyc%d_reg_set%d_%s%d',CycleLength, set_num, 'jitter-scale', jitters(i)*100);
%                     outputfile = [outputFolder,'\', file];
%                     savefile_writeaudio(t_scale(i,:), outputfile) % save t_scale
%                 end
%                 
            case 'No'
            case 'Cancel'
        end
    else
        [t_jitter, t_reg, t_rand, jitters,reg_interval] = NoisePattern_ke_jitters_gui(set_num, CycleLength, jitter_type,outputFolder);
        t_jitter_all = [t_reg; flipud(t_jitter); t_rand];
        plot_patternSeq_gui(t_jitter_all, CycleLength, jitters, jitter_type)
        answer = questdlg('Would you like to save the sequence and audio file?');
        switch answer
            case 'Yes'
                saveas(gcf, ['noisePattern_Cycle', num2str(CycleLength), '_set', num2str(set_num), '.png'])
                save(['noisePattern_Cycle', num2str(CycleLength), '_set', num2str(set_num), '.mat'], 't_jitter_all', 'reg_interval')
            case 'No'
            case 'Cancel'
        end
    end

end

if handles.play_audio_tg.Value % play the sound after 
    [y, Fs] = audioread(file_sound);
    soundsc(y,Fs)
end


function output_directory_tg_Callback(hObject, eventdata, handles)
% hObject    handle to output_directory_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of output_directory_tg as text
%        str2double(get(hObject,'String')) returns contents of output_directory_tg as a double


% --- Executes during object creation, after setting all properties.
function output_directory_tg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to output_directory_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function set_num_tg_Callback(hObject, eventdata, handles)
% hObject    handle to set_num_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of set_num_tg as text
%        str2double(get(hObject,'String')) returns contents of set_num_tg as a double
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function set_num_tg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_num_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cycle_length_tg_Callback(hObject, eventdata, handles)
% hObject    handle to cycle_length_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cycle_length_tg as text
%        str2double(get(hObject,'String')) returns contents of cycle_length_tg as a double
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function cycle_length_tg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cycle_length_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in play_audio_tg.
function play_audio_tg_Callback(hObject, eventdata, handles)
% hObject    handle to play_audio_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of play_audio_tg
% Update handles structure
guidata(hObject, handles);



function jitter_tg_Callback(hObject, eventdata, handles)
% hObject    handle to jitter_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of jitter_tg as text
%        str2double(get(hObject,'String')) returns contents of jitter_tg as a double
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function jitter_tg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to jitter_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in classic.
function classic_Callback(hObject, eventdata, handles)
% hObject    handle to classic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of classic
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in separated_tg.
function separated_tg_Callback(hObject, eventdata, handles)
% hObject    handle to separated_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of separated_tg
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in jitter_fixed.
function jitter_fixed_Callback(hObject, eventdata, handles)
% hObject    handle to jitter_fixed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of jitter_fixed


% --- Executes on button press in jitter2_rand.
function jitter2_rand_Callback(hObject, eventdata, handles)
% hObject    handle to jitter2_rand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of jitter2_rand


% --- Executes on button press in scale_tg.
function scale_tg_Callback(hObject, eventdata, handles)
% hObject    handle to scale_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of scale_tg


% --- Executes on button press in swap_tg.
function swap_tg_Callback(hObject, eventdata, handles)
% hObject    handle to swap_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of swap_tg


% --- Executes on button press in complex_tg.
function complex_tg_Callback(hObject, eventdata, handles)
% hObject    handle to complex_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in tone_pattern.
function tone_pattern_Callback(hObject, eventdata, handles)
% hObject    handle to tone_pattern (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
repeats_num = 10;
trials_num  = 100
set_num = str2num(handles.set_num_tg.String);
CycleLength = str2num(handles.cycle_length_tg.String);
outputFolder = handles.output_directory_tg.String;
cd([outputFolder,'\frequence_patterns'])
generate_frequency_pattern(set_num, CycleLength, repeats_num, trials_num);


% --- Executes on button press in write_audio.
function write_audio_Callback(hObject, eventdata, handles)
% hObject    handle to write_audio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
edit generate_audio_noise_pattern.m



function window_tg_Callback(hObject, eventdata, handles)
% hObject    handle to window_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of window_tg as text
%        str2double(get(hObject,'String')) returns contents of window_tg as a double


% --- Executes during object creation, after setting all properties.
function window_tg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to window_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in manually_tg.
function manually_tg_Callback(hObject, eventdata, handles)
% hObject    handle to manually_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% [FileName,PathName] = uigetfile('E:\Ke_Chen\Processed Data\Rach Recording\*.mat','Select the rastered spike files')
[FileName,PathName] = uigetfile('*.mat','Select the rastered spike files')

path = [PathName,FileName];
manually_correct_spikes(path)


% --- Executes on button press in example_tg.
function example_tg_Callback(hObject, eventdata, handles)
% hObject    handle to example_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uigetfile('E:\Ke_Chen\Processed Data\Rach Recording\*.mat','Select the rastered spike files')
path = [PathName,FileName];
load(path)
load('E:\Ke_Chen\MATLAB\Ephys\EphysAnalysis\BW.mat')
n_spike = length(spikedata);
for i = 1:n_spike
    if spikedata(i).keep ==0
        spikedata(i).example =0;
    else
        figure(100)
        psth1 = spikedata(i).clusterData.psth;
        psth_plot(psth1,1, 1:size(psth1.scmatrix, 2),0)
        set(gcf,'position',[100,200,600,800])
        fr_avg= (sum(psth1.scmatrix(:))/size(psth1.scmatrix, 1))/(size(psth1.scmatrix, 2)/1000);
        fprintf('Mean firing rate:%4.2f\n', fr_avg)
        x = input(['Good example or not', ' ', num2str(i), '/', num2str(n_spike), ': 1 or 0\n'])
        switch x
            case 1
                spikedata(i).example = 1;
            case 0
                spikedata(i).example = 0;
        end
    end
end
save(path, 'spikedata')


% --- Executes on button press in analysis_tool_tg.
function analysis_tool_tg_Callback(hObject, eventdata, handles)
% hObject    handle to analysis_tool_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
analysis_toolbox


% --- Executes on button press in rename_tg.
function rename_tg_Callback(hObject, eventdata, handles)
% hObject    handle to rename_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
edit rename_files


% --- Executes on button press in simulation_tg.
function simulation_tg_Callback(hObject, eventdata, handles)
% hObject    handle to simulation_tg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
edit simulation_FRs
