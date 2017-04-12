function varargout = BP_annotate_GUI(varargin)
% BP_ANNOTATE_GUI MATLAB code for BP_annotate_GUI.fig
%      BP_ANNOTATE_GUI, by itself, creates a new BP_ANNOTATE_GUI or raises the existing
%      singleton*.
%
%      H = BP_ANNOTATE_GUI returns the handle to a new BP_ANNOTATE_GUI or the handle to
%      the existing singleton*.
%
%      BP_ANNOTATE_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BP_ANNOTATE_GUI.M with the given input arguments.
%
%      BP_ANNOTATE_GUI('Property','Value',...) creates a new BP_ANNOTATE_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BP_annotate_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BP_annotate_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BP_annotate_GUI

% Last Modified by GUIDE v2.5 06-Apr-2017 16:42:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BP_annotate_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @BP_annotate_GUI_OutputFcn, ...
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


% --- Executes just before BP_annotate_GUI is made visible.
function BP_annotate_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
addpath('/Users/alexandrel/Codes/ClinicLab/BP_annotate');
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BP_annotate_GUI (see VARARGIN)

% Choose default command line output for BP_annotate_GUI
handles.output = hObject;
handles.Colors = [51, 19, 228;...
    246,108,254;...
    230,1,128;...
    106,3,239;...
    161,242,255]./255;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes BP_annotate_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = BP_annotate_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in loadbutton.
function loadbutton_Callback(hObject, eventdata, handles)
% hObject    handle to loadbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% actually load the data
handles.infocentertext.String = 'Loading...';
if isfield(handles,'data')
    handles = rmfield(handles,'data');
end
[handles.FileName,handles.PathName,~] = uigetfile({'*.txt; *.mat','Compatible files'},'Find the data file.');
FilePath = fullfile(handles.PathName, handles.FileName);
handles.data.waveform = load(FilePath);

% ask about data parameters, set them
prompt = {'Enter acquisition frequency (Hz):','Is the data clean? (Y or N):','Enter data units','Very high notch? (Y or N)'};
dlg_title = 'Data parameters';
num_lines = 1;
defaultans = {'1000','Y','mmHg','N'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
handles.data.Fs = str2num(answer{1});
handles.data.isClean = 0;
if regexp(answer{2},'y','ignorecase')
    handles.data.isClean = 1;
end
handles.data.units = 'other';
if regexp(answer{3},'mmhg','ignorecase')
    handles.data.units = 'mmHg';
else
    handles.infocentertext.String = fprintf('Data format %s is unsupported, contact administrator to add support.',answer{3});
end
handles.data.highNotch = 0;
if regexp(answer{4},'y','ignorecase')
    handles.data.highNotch = 1;
end
handles.data.time = (0:length(handles.data.waveform)-1)./handles.data.Fs;

handles = updateAxes(hObject, handles);

handles.infocentertext.String = 'Loading...done.';
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in annotatebutton.
function annotatebutton_Callback(hObject, eventdata, handles)

handles.infocentertext.String = 'Performing annotations...';
if ~handles.data.highNotch
    [ footIndex, systolicIndex, notchIndex, dicroticIndex, time, waveform ] =...
        BP_annotate( handles.data.waveform, handles.data.Fs, 0, handles.data.units, handles.data.isClean );
else
    [ footIndex, systolicIndex, notchIndex, dicroticIndex, time, waveform ] =...
        BP_annotate_high_notch( handles.data.waveform, handles.data.Fs, 0, handles.data.units, handles.data.isClean );
end
handles.data.bp_annotation.footIndex = footIndex;
handles.data.bp_annotation.systolicIndex = systolicIndex;
handles.data.bp_annotation.notchIndex = notchIndex;
handles.data.bp_annotation.dicroticIndex = dicroticIndex;
handles.data.bp_annotation.time = time;
handles.data.bp_annotation.waveform = waveform;
handles.infocentertext.String = 'Performing annotations...done.';

handles = updateAxes(hObject, handles);
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in exportbutton.
function exportbutton_Callback(hObject, eventdata, handles)
% hObject    handle to exportbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in clearbutton.
function clearbutton_Callback(hObject, eventdata, handles)
% hObject    handle to clearbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function handles = updateAxes(hObject, handles)
axes(handles.axes1)
plot(handles.data.time, handles.data.waveform)
xlim([handles.data.time(1) handles.data.time(end)])
box off
xlabel('time (s)')
ylabel(['pressure (',handles.data.units,')'])
legend('waveform')
if isfield(handles.data,'bp_annotation')
    cla; hold on;
    plot(handles.data.bp_annotation.time, handles.data.bp_annotation.waveform,'Color',handles.Colors(1,:))
    plot(handles.data.bp_annotation.time(handles.data.bp_annotation.footIndex),...
        handles.data.bp_annotation.waveform(handles.data.bp_annotation.footIndex),...
        'v','Color',handles.Colors(2,:),'MarkerFaceColor',handles.Colors(2,:));
    plot(handles.data.bp_annotation.time(handles.data.bp_annotation.systolicIndex),...
        handles.data.bp_annotation.waveform(handles.data.bp_annotation.systolicIndex),...
        '^','Color',handles.Colors(3,:),'MarkerFaceColor',handles.Colors(3,:));
    plot(handles.data.bp_annotation.time(handles.data.bp_annotation.notchIndex),...
        handles.data.bp_annotation.waveform(handles.data.bp_annotation.notchIndex),...
        '>','Color',handles.Colors(4,:),'MarkerFaceColor',handles.Colors(4,:));
    legend({'waveform','foot','systole','notch'})
end
    

linkaxes([handles.axes1,handles.axes5,handles.axes4]);
% Update handles structure
guidata(hObject, handles);