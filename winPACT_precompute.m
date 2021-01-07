% winPACT_precompute() - This is a plugin to precompute cross-frequency
%                        phase-amplitude coupling (PAC). It applies windows at
%                        specified time and compute PAC.
%
% 

% History:
% 07/31/2018 Makoto. Updated to support selection of EEGLAB events.
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



function varargout = winPACT_precompute(varargin)
% WINPACT_PRECOMPUTE MATLAB code for winPACT_precompute.fig
%      WINPACT_PRECOMPUTE, by itself, creates a new WINPACT_PRECOMPUTE or raises the existing
%      singleton*.
%
%      H = WINPACT_PRECOMPUTE returns the handle to a new WINPACT_PRECOMPUTE or the handle to
%      the existing singleton*.
%
%      WINPACT_PRECOMPUTE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WINPACT_PRECOMPUTE.M with the given input arguments.
%
%      WINPACT_PRECOMPUTE('Property','Value',...) creates a new WINPACT_PRECOMPUTE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before winPACT_precompute_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to winPACT_precompute_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help winPACT_precompute

% Last Modified by GUIDE v2.5 03-Aug-2018 15:55:43

% 05/08/2018 Makoto. Modified.
% 

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @winPACT_precompute_OpeningFcn, ...
                   'gui_OutputFcn',  @winPACT_precompute_OutputFcn, ...
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


% --- Executes just before winPACT_precompute is made visible.
function winPACT_precompute_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to winPACT_precompute (see VARARGIN)

% Choose default command line output for winPACT_precompute
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Initial setup
set(gcf, 'Name', 'winPACT_precompute()');
set(handles.savePathEdit, 'String', '', 'BackgroundColor', [1 1 1]);
set(handles.loadXlsxPathEdit, 'BackgroundColor', [0.9 0.9 0.9]);

% Obtain EEG.
EEG = evalin('base', 'EEG');

% Extract unique event types.
if isempty(EEG.event)
    uniqueEvents = '(none)';
else
    uniqueEvents = unique({EEG.event.type}');
end

% Plot the event names in the popup menu.
set(handles.eventTypeListbox, 'string', uniqueEvents);
set(handles.eventTypeListbox, 'max', length(uniqueEvents));


% UIWAIT makes winPACT_precompute wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = winPACT_precompute_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function phaseFreqEdit_Callback(hObject, eventdata, handles)
% hObject    handle to phaseFreqEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of phaseFreqEdit as text
%        str2double(get(hObject,'String')) returns contents of phaseFreqEdit as a double


% --- Executes during object creation, after setting all properties.
function phaseFreqEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phaseFreqEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ampFreqEdit_Callback(hObject, eventdata, handles)
% hObject    handle to ampFreqEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ampFreqEdit as text
%        str2double(get(hObject,'String')) returns contents of ampFreqEdit as a double


% --- Executes during object creation, after setting all properties.
function ampFreqEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ampFreqEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function windowSizeEdit_Callback(hObject, eventdata, handles)
% hObject    handle to windowSizeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of windowSizeEdit as text
%        str2double(get(hObject,'String')) returns contents of windowSizeEdit as a double


% --- Executes during object creation, after setting all properties.
function windowSizeEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to windowSizeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in loadXlsxPushbutton.
function loadXlsxPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to loadXlsxPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Clear and pseudo-disable 'eventTypeListbox'.
set(handles.loadXlsxPathEdit,      'BackgroundColor', [1 1 1]);
set(handles.eventTypeListbox, 'BackgroundColor', [0.9 0.9 0.9]);

% Open GUI to obtain the file path.
[fileName, pathName] = uigetfile({'*.xls'; '*.xlsx'}, 'Select .xls file to import event times');
if any(fileName)
    set(handles.loadXlsxPathEdit, 'String', [pathName fileName])
else
    set(handles.loadXlsxPathEdit, 'String', '');
    disp('Cancelled.')
end



function loadXlsxPathEdit_Callback(hObject, eventdata, handles)
% hObject    handle to loadXlsxPathEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of loadXlsxPathEdit as text
%        str2double(get(hObject,'String')) returns contents of loadXlsxPathEdit as a double

% Clear and pseudo-disable 'eventTypeListbox'.
set(handles.loadXlsxPathEdit,      'BackgroundColor', [1 1 1]);
set(handles.eventTypeListbox, 'BackgroundColor', [0.9 0.9 0.9]);




% --- Executes during object creation, after setting all properties.
function loadXlsxPathEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loadXlsxPathEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in savePathPushbutton.
function savePathPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to savePathPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[FILENAME, PATHNAME] = uiputfile('.xls', 'Select a folder to save the file');
set(handles.savePathEdit, 'String', [PATHNAME FILENAME]);



function savePathEdit_Callback(hObject, eventdata, handles)
% hObject    handle to savePathEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of savePathEdit as text
%        str2double(get(hObject,'String')) returns contents of savePathEdit as a double


% --- Executes during object creation, after setting all properties.
function savePathEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to savePathEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function saveFileNameEdit_Callback(hObject, eventdata, handles)
% hObject    handle to saveFileNameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of saveFileNameEdit as text
%        str2double(get(hObject,'String')) returns contents of saveFileNameEdit as a double


% --- Executes during object creation, after setting all properties.
function saveFileNameEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to saveFileNameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in startPushbutton.
function startPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to startPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Obtain data.
EEG = evalin('base', 'EEG');
if isempty(EEG.event)
    EEG.event = struct('type', 'boundary', 'latency', 0);
end

% Set frequency bandpass edges.
lowFreqForPhase = str2num(get(handles.phaseFreqEdit, 'String'));
highFreqForAmp  = str2num(get(handles.ampFreqEdit,   'String'));

% Preserve the result container.
windowLengthS = str2num(get(handles.windowSizeEdit, 'String'));

% Check is the window size is longer than one cycle of low frequency
if 1/lowFreqForPhase(1) >= windowLengthS
    error('Phase frequency is too low for the window size. Must be 1/phase_freq < window_length.')
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Obtain window center time %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlsxPath = get(handles.loadXlsxPathEdit, 'String');
if any(xlsxPath)
    
    % Load .xlsx file.
    if strfind(xlsxPath, '.xls')
        eventTimesSec = xlsread(xlsxPath);
                
    % Generate number series.
    else
        eventTimesSec = eval(xlsxPath);
    end
    
    % Select event markers.
else
    allEventTypes = {EEG.event.type}';
    selectedConditionIdx     = get(handles.eventTypeListbox, 'value');
    conditionLabels          = get(handles.eventTypeListbox, 'string');
    selectedConditionLabels  = conditionLabels(selectedConditionIdx);
    selectedEventIdx         = cellfun(@(x) strcmp(allEventTypes, x), selectedConditionLabels, 'uniformoutput', false);
    selectedEventIdxCombined = cell2mat(cellfun(@(x) find(x), selectedEventIdx, 'uniformoutput', false));
    selectedEventFrame = [EEG.event(selectedEventIdxCombined).latency];
    eventTimesSec = EEG.times(selectedEventFrame)/1000;
end

% Make sure eventTimeSec is a column vector.
if size(eventTimesSec,2)>size(eventTimesSec,1)
    eventTimesSec = eventTimesSec';
end

% Check the window coverage.
eventTimesSec = sort(eventTimesSec);
currentWindowFrameEdges      = [round(eventTimesSec-windowLengthS/2) round(eventTimesSec+windowLengthS/2)]*EEG.srate;
currentWindowFrameEdges(:,1) = currentWindowFrameEdges(:,1)+1;
if EEG.times(currentWindowFrameEdges(1,1)) < EEG.times(1)
    error('The latency of the first analysis window is outside the data time range.')
end
if currentWindowFrameEdges(end,2) > EEG.pnts
    error('The latency of the last analysis window is outside the data time range.')
end


    
%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate PAC. %%%
%%%%%%%%%%%%%%%%%%%%%%
numPhaseBinsPerPi     = 18; % 36 bins for 2pi. Tort et al. (2010) used 18 bins for 2pi.
numIterations         = str2num(get(handles.surroStatIterEdit, 'String'));
canoltysMIAllChan     = zeros(EEG.nbchan, length(eventTimesSec));
windowMeanAmpAllChan  = zeros(EEG.nbchan, length(eventTimesSec));
MInormAllChan         = zeros(EEG.nbchan, length(eventTimesSec));
pValueAllChan         = zeros(EEG.nbchan, length(eventTimesSec));
klDistAllChan         = zeros(EEG.nbchan, length(eventTimesSec));
ampDistribAllChan     = zeros(EEG.nbchan, length(eventTimesSec), numPhaseBinsPerPi*2); % Because 2*pi.
genMaxStatAllChan     = zeros(EEG.nbchan, length(eventTimesSec));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Report designs of the band-pass filters. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Report low-frequency phase filter.
if length(lowFreqForPhase) == 3
    phaseTbw = lowFreqForPhase(3);
else
    phaseTbw = lowFreqForPhase(1)*0.5;
end
phaseFilterOrder = pop_firwsord('hamming', EEG.srate, phaseTbw);

disp(sprintf('\n'))
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(sprintf('Bandpass filter parameters for low-frequency phase:'))
disp(sprintf('Filter type: Zero-phase FIR filter (Hamming window)'))
disp(sprintf('Filter order: %.0f (Sampling rate: %.0f [Hz])', phaseFilterOrder, EEG.srate))
disp(sprintf('Pass-band edges: %.2f to %.2f [Hz]', lowFreqForPhase(1), lowFreqForPhase(2)))
disp(sprintf('Cutoff (-6dB) frequencies: %.2f to %.2f [Hz]', lowFreqForPhase(1)-phaseTbw/2, lowFreqForPhase(2)+phaseTbw/2))
disp(sprintf('Transition band width: %.2f [Hz]', phaseTbw))
disp(sprintf('\n'))

% Report high-frequency amplitude filter.
if length(highFreqForAmp) == 3
    ampTbw = highFreqForAmp(3);
else
    ampTbw = highFreqForAmp(1)*0.2;
end
ampFilterOrder = pop_firwsord('hamming', EEG.srate, ampTbw);
disp(sprintf('Bandpass filter parameters for high-frequency amplitude:'))
disp(sprintf('Filter type: Zero-phase FIR filter (Hamming window)'))
disp(sprintf('Filter order: %.0f (Sampling rate: %.0f [Hz])', ampFilterOrder, EEG.srate))
disp(sprintf('Pass-band edges: %.2f to %.2f [Hz]', highFreqForAmp(1), highFreqForAmp(2)))
disp(sprintf('Cutoff (-6dB) frequencies: %.2f to %.2f [Hz]', highFreqForAmp(1)-ampTbw/2, highFreqForAmp(2)+ampTbw/2))
disp(sprintf('Transition band width: %.2f [Hz]', ampTbw))
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(sprintf('\n'))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute PAC for all the channels. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Entering main loop.')
for currentChannelIdx = 1:EEG.nbchan
    
    % Start the timer.
    startedTime = tic;
    
    % Select the current channel. This is to use pop_firws.
    singleChannelEEG = pop_select(EEG, 'channel', currentChannelIdx);
    
    % Extract low-frequency phase.
    phaseEEG      = pop_firws(singleChannelEEG, 'fcutoff',  [lowFreqForPhase(1) lowFreqForPhase(2)], 'ftype', 'bandpass', 'wtype', 'hamming', 'forder', phaseFilterOrder, 'minphase', 0);
    analyticPhase = angle(hilbert(phaseEEG.data'));
    clear phaseEEG
    
    % Extract high-frequency amplitude.
    ampEEG        = pop_firws(singleChannelEEG, 'fcutoff',  [highFreqForAmp(1) highFreqForAmp(2)], 'ftype', 'bandpass', 'wtype', 'hamming', 'forder', ampFilterOrder, 'minphase', 0);
    analyticAmp   = abs(hilbert(ampEEG.data'));
    clear ampEEG
    
    % Compute low-freq phase coupled with high-freq amp. 
    Z = analyticAmp.* exp(1i*analyticPhase);
   
    % Comnpute PAC.
    for eventIdx = 1:length(eventTimesSec)
        
        % Prepare the current window phase, amp, and Z.
        currentWindowPhase = analyticPhase(currentWindowFrameEdges(eventIdx,1)+1:currentWindowFrameEdges(eventIdx,2));
        currentWindowAmp   = analyticAmp(currentWindowFrameEdges(eventIdx,1)+1:currentWindowFrameEdges(eventIdx,2));
        currentWindowZ     = Z(currentWindowFrameEdges(eventIdx,1)+1:currentWindowFrameEdges(eventIdx,2));

        % Compute Canolty's Modulation Index (MI). Using mean instead of sum.
        canoltysMI = abs(mean(currentWindowZ));
        
        % Compute Oezkurt's Normalized Modulation Index (MInorm). Taken from Oezkurt's code.
        surroFrameIdx = zeros(length(currentWindowZ), numIterations);
        for surroIdx = 1:numIterations
            surroFrameIdx(:,surroIdx) = circshift([1:length(currentWindowZ)]', [1+floor(rand(1)*length(currentWindowZ)-1) 1]);
        end
        surroAmp   = repmat(currentWindowAmp,   [1 numIterations]);
        surroPhase = repmat(currentWindowPhase, [1 numIterations]);
        surroAmp   = surroAmp(surroFrameIdx);
        surroMi    = abs(mean(surroAmp.*exp(1i*surroPhase),1));
        [surroMiMean, surroMiStd] = normfit(surroMi);
        MInorm = (canoltysMI - surroMiMean)/surroMiStd; % Converting to a z-score.
        pValue = 1 - normcdf(canoltysMI, surroMiMean, surroMiStd);
        
        % Compute generalized maximum statistics for multiple comparison correction.
        surroZ         = (surroMi - surroMiMean)/surroMiStd;
        [~, maxIdx]    = max(surroZ);
        surroZ(maxIdx) = nan;
        genMaxStat     = max(surroZ);        
        
        % Compute Tort's Kullbak-Leibler divergence. Normalize currentBinAmp (Tort et al., 2011)
        freqBins = -pi:pi/numPhaseBinsPerPi:pi; 
        [~, freqBinIdx] = histc(currentWindowPhase, freqBins);
        amplitudeDistributionOverPhaseBins = zeros(1,length(freqBins)-1);
        for binIdxIdx = 1:length(freqBins)-1
            currentBinIdx = find(freqBinIdx == binIdxIdx);
            amplitudeDistributionOverPhaseBins(binIdxIdx) = mean(currentWindowAmp(currentBinIdx));
        end
        amplitudeDistributionOverPhaseBins = amplitudeDistributionOverPhaseBins/sum(amplitudeDistributionOverPhaseBins);
        klDistance = KLDiv(amplitudeDistributionOverPhaseBins, ones(1,length(amplitudeDistributionOverPhaseBins)));
    
        % Update the output structure
        windowMeanAmpAllChan(currentChannelIdx,eventIdx) = mean(currentWindowAmp);
        canoltysMIAllChan(currentChannelIdx,eventIdx)    = canoltysMI;
        MInormAllChan(currentChannelIdx,eventIdx)        = MInorm;
        pValueAllChan(currentChannelIdx,eventIdx)        = pValue;
        klDistAllChan(currentChannelIdx,eventIdx)        = klDistance;
        ampDistribAllChan(currentChannelIdx,eventIdx,:)  = amplitudeDistributionOverPhaseBins;
        genMaxStatAllChan(currentChannelIdx,eventIdx)    = genMaxStat;
    end
    
    % Show time lapse.
    timeElapsed = toc(startedTime);
    expectedHoursToFinish = timeElapsed*(EEG.nbchan-currentChannelIdx)/3600;
    disp(sprintf('\n\n%.0f/%.0fch done. %.2f hours left.\n\n', currentChannelIdx, EEG.nbchan, expectedHoursToFinish));
end

% Compute generalized FWER correction. p < 0.05 is highlighted by *-1.
criticalZscore = prctile(genMaxStatAllChan(:), 95);
gFWER_mask     = MInormAllChan > criticalZscore;
correctionMask = pValueAllChan.*(-1*gFWER_mask);
pValueAllChan(gFWER_mask) = correctionMask(gFWER_mask);

% Add latency on the top raw
windowMeanAmpAllChan = [eventTimesSec'; windowMeanAmpAllChan];
canoltysMIAllChan    = [eventTimesSec'; canoltysMIAllChan];
MInormAllChan        = [eventTimesSec'; MInormAllChan];
pValueAllChan        = [eventTimesSec'; pValueAllChan];
klDistAllChan        = [eventTimesSec'; klDistAllChan];
ampDistribAllChan    = permute(ampDistribAllChan, [3 1 2]);
ampDistribAllChan    = reshape(ampDistribAllChan, [size(ampDistribAllChan,1)*size(ampDistribAllChan,2) size(ampDistribAllChan,3)]);
ampDistribAllChan    = [eventTimesSec'; ampDistribAllChan];

% Store the results under EEG.etc.winPACT
EEG.etc.winPACT.windowMeanAmpAllChan = windowMeanAmpAllChan;
EEG.etc.winPACT.canoltysMIAllChan    = canoltysMIAllChan;
EEG.etc.winPACT.MInormAllChan        = MInormAllChan;
EEG.etc.winPACT.pValueAllChan        = pValueAllChan;
EEG.etc.winPACT.klDistAllChan        = klDistAllChan;
EEG.etc.winPACT.ampDistribAllChan    = ampDistribAllChan;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Update EEG in the base workspace. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
assignin('base', 'EEG', EEG);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Export the results as multiple .xls files if requested. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(get(handles.savePathEdit, 'String'))
    saveFullpath = get(handles.savePathEdit, 'String');
    saveFullpath = saveFullpath(1:end-4);
    savePath = [saveFullpath '_AMP'];
    xlswrite(savePath, windowMeanAmpAllChan);
    savePath = [saveFullpath '_MI'];
    xlswrite(savePath, canoltysMIAllChan);
    savePath = [saveFullpath '_MInorm'];
    xlswrite(savePath, MInormAllChan);
    savePath = [saveFullpath '_pValue'];
    xlswrite(savePath, pValueAllChan);
    savePath = [saveFullpath '_KLD'];
    xlswrite(savePath, klDistAllChan);
    savePath = [saveFullpath '_PHASE'];
    xlswrite(savePath, ampDistribAllChan);
end



% Done.
disp('Done.')



function surroStatIterEdit_Callback(hObject, eventdata, handles)
% hObject    handle to surroStatIterEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of surroStatIterEdit as text
%        str2double(get(hObject,'String')) returns contents of surroStatIterEdit as a double



% --- Executes during object creation, after setting all properties.
function surroStatIterEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to surroStatIterEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dist = KLDiv(P,Q)
%  dist = KLDiv(P,Q) Kullback-Leibler divergence of two discrete probability
%  distributions
%  P and Q  are automatically normalised to have the sum of one on rows
% have the length of one at each 
% P =  n x nbins
% Q =  1 x nbins or n x nbins(one to one)
% dist = n x 1
% Taken from: https://www.mathworks.com/matlabcentral/fileexchange/20688-kullback-leibler-divergence

if size(P,2)~=size(Q,2)
    error('the number of columns in P and Q should be the same');
end

if sum(~isfinite(P(:))) + sum(~isfinite(Q(:)))
   error('the inputs contain non-finite values!') 
end

% normalizing the P and Q
if size(Q,1)==1
    Q = Q ./sum(Q);
    P = P ./repmat(sum(P,2),[1 size(P,2)]);
    temp =  P.*log(P./repmat(Q,[size(P,1) 1]));
    temp(isnan(temp))=0;% resolving the case when P(i)==0
    dist = sum(temp,2);
        
elseif size(Q,1)==size(P,1)
    
    Q = Q ./repmat(sum(Q,2),[1 size(Q,2)]);
    P = P ./repmat(sum(P,2),[1 size(P,2)]);
    temp =  P.*log(P./Q);
    temp(isnan(temp))=0; % resolving the case when P(i)==0
    dist = sum(temp,2);
end


                % % % % --- Executes on button press in selectEventPushbutton.
                % % % function selectEventPushbutton_Callback(hObject, eventdata, handles)
                % % % % hObject    handle to selectEventPushbutton (see GCBO)
                % % % % eventdata  reserved - to be defined in a future version of MATLAB
                % % % % handles    structure with handles and user data (see GUIDATA)
                % % % 
                % % % % Clear and disable 'load xls' edit.
                % % % set(handles.selectEventEdit, 'Enable', 'on');
                % % % set(handles.loadXlsxPathEdit,     'String', '');
                % % % set(handles.loadXlsxPathEdit,     'Enable', 'off');
                % % % 
                % % % selectedEvents = winPACT_eventListPopupmenu;
                % % % disp('test')
                % % % 
                % % % function selectEventEdit_Callback(hObject, eventdata, handles)
                % % % % hObject    handle to selectEventEdit (see GCBO)
                % % % % eventdata  reserved - to be defined in a future version of MATLAB
                % % % % handles    structure with handles and user data (see GUIDATA)
                % % % 
                % % % % Hints: get(hObject,'String') returns contents of selectEventEdit as text
                % % % %        str2double(get(hObject,'String')) returns contents of selectEventEdit as a double
                % % % 
                % % % 
                % % % 
                % % % % --- Executes during object creation, after setting all properties.
                % % % function selectEventEdit_CreateFcn(hObject, eventdata, handles)
                % % % % hObject    handle to selectEventEdit (see GCBO)
                % % % % eventdata  reserved - to be defined in a future version of MATLAB
                % % % % handles    empty - handles not created until after all CreateFcns called
                % % % 
                % % % % Hint: edit controls usually have a white background on Windows.
                % % % %       See ISPC and COMPUTER.
                % % % if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
                % % %     set(hObject,'BackgroundColor','white');
                % % % end


                
% --- Executes on selection change in eventTypeListbox.
function eventTypeListbox_Callback(hObject, eventdata, handles)
% hObject    handle to eventTypeListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns eventTypeListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from eventTypeListbox

set(handles.eventTypeListbox, 'BackgroundColor', [1 1 1]);
set(handles.loadXlsxPathEdit, 'string', '');
set(handles.loadXlsxPathEdit, 'BackgroundColor', [0.9 0.9 0.9]);



% --- Executes during object creation, after setting all properties.
function eventTypeListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eventTypeListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
