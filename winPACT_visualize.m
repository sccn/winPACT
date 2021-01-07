% winPACT_visualize() - This is a plugin to visualize precomputed cross-frequency
%                       phase-amplitude coupling (PAC). 

% History:
% 08/06/2018 Makoto. Updated to be used as a general EEGLAB plugin.
% 05/04/2018 Makoto. Methods by Oezkurt and Tort supported. Phase histogram supported.
% 03/14/2017 Makoto. Window mean amplitude output.
% 03/10/2017 Makoto. Created.

% Copyright (C) 2017, Makoto Miyakoshi. SCCN, INC, UCSD.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1.07  USA


function varargout = winPACT_visualize(varargin)
% WINPACT_VISUALIZE MATLAB code for winPACT_visualize.fig
%      WINPACT_VISUALIZE, by itself, creates a new WINPACT_VISUALIZE or raises the existing
%      singleton*.
%
%      H = WINPACT_VISUALIZE returns the handle to a new WINPACT_VISUALIZE or the handle to
%      the existing singleton*.
%
%      WINPACT_VISUALIZE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WINPACT_VISUALIZE.M with the given input arguments.
%
%      WINPACT_VISUALIZE('Property','Value',...) creates a new WINPACT_VISUALIZE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before winPACT_visualize_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to winPACT_visualize_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help winPACT_visualize

% Last Modified by GUIDE v2.5 16-Aug-2018 14:54:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @winPACT_visualize_OpeningFcn, ...
                   'gui_OutputFcn',  @winPACT_visualize_OutputFcn, ...
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


% --- Executes just before winPACT_visualize is made visible.
function winPACT_visualize_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to winPACT_visualize (see VARARGIN)

% Initial setup
set(gcf, 'Name', 'winPACT_visualize()')

% Choose default command line output for winPACT_visualize
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes winPACT_visualize wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = winPACT_visualize_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on slider movement.
function timeSlider_Callback(hObject, eventdata, handles)
% hObject    handle to timeSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% Obtain the current time in second.
TIME  = evalin('base', 'TIME');
PHASE = evalin('base', 'PHASE');
timeSeriesData = evalin('base', 'timeSeriesData');
pValue = evalin('base', 'pValue');
sliderPosition = get(hObject,'Value');
currentTimeIdx = floor(length(TIME)*sliderPosition)+1;
if sliderPosition == 1
    currentTimeIdx = currentTimeIdx-1;
end
currentTimeSec = TIME(currentTimeIdx);

% Display the current time.
set(handles.timeInSecText, 'String', sprintf('%.1f s',currentTimeSec))

% Display the phaser.
currentChIdx       = str2num(get(handles.enterChannelIdxEdit, 'String'));
currentPhase       = PHASE(1+36*(currentChIdx-1):36*(currentChIdx), currentTimeIdx);
currentPhase0to2pi = currentPhase([18:end 1:18])';
maxValue = max(PHASE(:));
customPolarPlot(handles.phaserAxes, 0:2*pi/36:2*pi, currentPhase0to2pi, maxValue)

% Display the bar graph.
barData = [currentPhase0to2pi currentPhase0to2pi];
bar(handles.barAxes, barData, 1, 'facecolor', [0.93 0.96 1], 'edgecolor', [0.66 0.76 1])
set(handles.barAxes, 'xlim', [0.5 74.5], 'ylim', [0 maxValue],...
    'xtick', [1:18:74], 'xticklabel', {'0' '180' '360' '540' '720'}, 'box', 'off')
set(get(handles.barAxes, 'xlabel'), 'string', 'Phase (deg.)')
set(get(handles.barAxes, 'ylabel'), 'string', 'Normalized Amplitude')

% Display the vertical bar.
try
    delete(handles.timeCursor)
end
handles.timeCursor = line([currentTimeSec currentTimeSec], ylim, 'LineWidth', 2, 'Color', [0.76 0.76 0.76], 'linestyle', '--');

% Update handles structure
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function timeSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to timeSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end




% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function enterChannelIdxEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to enterChannelIdxEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function enterChannelIdxEdit_Callback(hObject, eventdata, handles)
% hObject    handle to enterChannelIdxEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of enterChannelIdxEdit as text
%        str2double(get(hObject,'String')) returns contents of enterChannelIdxEdit as a double

% Clear the axes.
cla(handles.phaserAxes)
cla(handles.barAxes)
cla(handles.timeSeriesAxes)

% Import variables from base workspace.
TIME           = evalin('base', 'TIME');
timeSeriesData = evalin('base', 'timeSeriesData');
pValue         = evalin('base', 'pValue');
currentChIdx   = str2num(get(handles.enterChannelIdxEdit, 'String'));

% Display significant frames.
significanceMask = find(pValue(currentChIdx,:)<0);
axes(handles.timeSeriesAxes)
hold on
for significantFrames = 1:length(significanceMask)
    line([TIME(significanceMask(significantFrames)) TIME(significanceMask(significantFrames))], [0 max(timeSeriesData(:))], 'color', [1 0.76 0.66], 'linewidth', 5)
end

% Display PAC time series.
plot(handles.timeSeriesAxes, TIME, timeSeriesData(currentChIdx,:), 'linewidth', 3, 'color', [0.66 0.66 1]);
xlim(handles.timeSeriesAxes, [min(TIME) max(TIME)])
ylim(handles.timeSeriesAxes, [0 max(timeSeriesData(:))])
set(handles.timeSeriesAxes, 'box', 'off')
hold off

% Update the the slider.
timeSlider_Callback(handles.timeSlider, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function selectTimeSeriesPopupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selectTimeSeriesPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in selectTimeSeriesPopupmenu.
function selectTimeSeriesPopupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to selectTimeSeriesPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns selectTimeSeriesPopupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from selectTimeSeriesPopupmenu

% Import EEG.
EEG = evalin('base', 'EEG');

% Clear the axes.
cla(handles.phaserAxes)
cla(handles.barAxes)
cla(handles.timeSeriesAxes)



%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load phase data. %%%
%%%%%%%%%%%%%%%%%%%%%%%%
phaseData = EEG.etc.winPACT.ampDistribAllChan;
TIME  = phaseData(1,:);
PHASE = phaseData(2:end,:);
assignin('base', 'TIME',  TIME);
assignin('base', 'PHASE', PHASE);



%%%%%%%%%%%%%%%%%%%%%%
%%% Load p-values. %%%
%%%%%%%%%%%%%%%%%%%%%%
pValue = EEG.etc.winPACT.pValueAllChan(2:end,:);
assignin('base', 'pValue', pValue);



%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load PAC measure. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
switch get(handles.selectTimeSeriesPopupmenu, 'Value')
    case 1
        timeSeriesData = EEG.etc.winPACT.windowMeanAmpAllChan(2:end,:);
        pacTypeDisplayHandle = findall(gcf, 'tag', 'pacTypeText');
        set(pacTypeDisplayHandle, 'String', 'HFO Amplitude (Red, gFWER-corrected p < 0.05 in Özkurt''s test across all frames & channels)')
    case 2
        timeSeriesData = EEG.etc.winPACT.canoltysMIAllChan(2:end,:);
        pacTypeDisplayHandle = findall(gcf, 'tag', 'pacTypeText');
        set(pacTypeDisplayHandle, 'String', 'Canolty''s MI (Red, gFWER-corrected p < 0.05 in Özkurt''s test across all frames & channels)')
    case 3
        timeSeriesData = EEG.etc.winPACT.MInormAllChan(2:end,:);
        pacTypeDisplayHandle = findall(gcf, 'tag', 'pacTypeText');
        set(pacTypeDisplayHandle, 'String', 'Özkurt''s normalized MI (Red, gFWER-corrected p < 0.05 in Özkurt''s test across all frames & channels)')
    case 4
        timeSeriesData = EEG.etc.winPACT.klDistAllChan(2:end,:);
        pacTypeDisplayHandle = findall(gcf, 'tag', 'pacTypeText');
        set(pacTypeDisplayHandle, 'String', 'Tort''s KL Divergence (Red, gFWER-corrected p < 0.05 in Özkurt''s test across all frames & channels)')
end
assignin('base', 'timeSeriesData', timeSeriesData);

% Update the channel index (which also updates the slider)
enterChannelIdxEdit_Callback(handles.enterChannelIdxEdit, eventdata, handles)
