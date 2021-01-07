function varargout = winPACT_simulate(varargin)
% WINPACT_SIMULATE MATLAB code for winPACT_simulate.fig
%      WINPACT_SIMULATE, by itself, creates a new WINPACT_SIMULATE or raises the existing
%      singleton*.
%
%      H = WINPACT_SIMULATE returns the handle to a new WINPACT_SIMULATE or the handle to
%      the existing singleton*.
%
%      WINPACT_SIMULATE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WINPACT_SIMULATE.M with the given input arguments.
%
%      WINPACT_SIMULATE('Property','Value',...) creates a new WINPACT_SIMULATE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before winPACT_simulate_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to winPACT_simulate_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help winPACT_simulate

% Last Modified by GUIDE v2.5 16-Aug-2018 12:56:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @winPACT_simulate_OpeningFcn, ...
                   'gui_OutputFcn',  @winPACT_simulate_OutputFcn, ...
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


% --- Executes just before winPACT_simulate is made visible.
function winPACT_simulate_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to winPACT_simulate (see VARARGIN)

% Choose default command line output for winPACT_simulate
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes winPACT_simulate wait for user response (see UIRESUME)
% uiwait(handles.winPACT_simulateFigure);


% --- Outputs from this function are returned to the command line.
function varargout = winPACT_simulate_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function lfoHzEdit_Callback(hObject, eventdata, handles)
% hObject    handle to lfoHzEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lfoHzEdit as text
%        str2double(get(hObject,'String')) returns contents of lfoHzEdit as a double


% --- Executes during object creation, after setting all properties.
function lfoHzEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lfoHzEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function hfoHzEdit_Callback(hObject, eventdata, handles)
% hObject    handle to hfoHzEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hfoHzEdit as text
%        str2double(get(hObject,'String')) returns contents of hfoHzEdit as a double


% --- Executes during object creation, after setting all properties.
function hfoHzEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hfoHzEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function noiseLevelEdit_Callback(hObject, eventdata, handles)
% hObject    handle to noiseLevelEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of noiseLevelEdit as text
%        str2double(get(hObject,'String')) returns contents of noiseLevelEdit as a double


% --- Executes during object creation, after setting all properties.
function noiseLevelEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to noiseLevelEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dataLengthEdit_Callback(hObject, eventdata, handles)
% hObject    handle to dataLengthEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dataLengthEdit as text
%        str2double(get(hObject,'String')) returns contents of dataLengthEdit as a double


% --- Executes during object creation, after setting all properties.
function dataLengthEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dataLengthEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in coloredNoisePopupmenu.
function coloredNoisePopupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to coloredNoisePopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns coloredNoisePopupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from coloredNoisePopupmenu


% --- Executes during object creation, after setting all properties.
function coloredNoisePopupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to coloredNoisePopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in generateDataPushbutton.
function generateDataPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to generateDataPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Obtain user input.
phaseLowFreqHz  = str2num(get(handles.lfoHzEdit, 'String'));
ampHighFreqHz   = str2num(get(handles.hfoHzEdit, 'String'));
noiseLevel      = str2num(get(handles.noiseLevelEdit, 'String'));
dataLengthInSec = str2num(get(handles.dataLengthEdit, 'String'));
coloredNoise    = get(handles.coloredNoisePopupmenu, 'value');
pacPresentWindowEdges = str2num(get(handles.pacPresentEdit, 'String'));

% Generate data.
switch coloredNoise
    case 1
        coloredNoise = 1;
    case 2
        coloredNoise = 0;
end
[lfoHfoNoise, SNR] = Oezkurt2011_synthesize_pac_modified(phaseLowFreqHz, ampHighFreqHz, noiseLevel, dataLengthInSec, coloredNoise);
disp(sprintf('\n\n SNR %.fdB\n\n', SNR))

% Construct data.
onlyNoise        = lfoHfoNoise(3,:);
simulatedPacData = onlyNoise;
lfoPlusNoise     = onlyNoise;
hfoPlusNoise     = onlyNoise;
pacPresentWindowEdgesFrames = round(pacPresentWindowEdges*1000);
for numWindows = 1:size(pacPresentWindowEdgesFrames,1);
    currentWindow = pacPresentWindowEdgesFrames(numWindows,1):pacPresentWindowEdgesFrames(numWindows,2);
    simulatedPacData(currentWindow) = sum(lfoHfoNoise(:,    currentWindow));
    lfoPlusNoise(    currentWindow) = sum(lfoHfoNoise([1 3],currentWindow));
    hfoPlusNoise(    currentWindow) = sum(lfoHfoNoise([2 3],currentWindow));
end
outputData = cat(1, simulatedPacData, lfoPlusNoise, hfoPlusNoise, onlyNoise);

% Export it to EEGLAB.
EEG = eeg_emptyset;
EEG.setname = sprintf('LFO%.1fHz HFO%.0fHz SNR%.fdB PAC present in %s(s)', phaseLowFreqHz, ampHighFreqHz, SNR, get(handles.pacPresentEdit, 'String'));
EEG.subject = 'Simulated data';
EEG.nbchan  = 4;
EEG.trials  = 1;
EEG.pnts    = size(outputData,2);
EEG.srate   = 1000;
EEG.xmin    = 0.001;
EEG.xmax    = EEG.pnts/1000;
EEG.times   = 0.001:0.001:EEG.xmax;
EEG.data    = outputData;
assignin('base', 'EEG', EEG);

% Close this GUI
close(handles.winPACT_simulateFigure)
eeglab redraw



function pacPresentEdit_Callback(hObject, eventdata, handles)
% hObject    handle to pacPresentEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pacPresentEdit as text
%        str2double(get(hObject,'String')) returns contents of pacPresentEdit as a double


% --- Executes during object creation, after setting all properties.
function pacPresentEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pacPresentEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
