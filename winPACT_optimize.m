% winPACT_optimize() - This GUI environment allows users to optimize and
%                      validate the choice of window length.  

% History:
% 08/17/2018 Makoto. Modified for optimizing and validating the selection of the window length.
% 07/31/2018 Makoto. Updated to support selection of EEGLAB events.
% 05/04/2018 Makoto. Methods by Oezkurt and Tort supported. Phase histogram supported.
% 03/14/2017 Makoto. Window mean amplitude output.
% 03/10/2017 Makoto. Created.

% Copyright (C) 2018, Makoto Miyakoshi. SCCN, INC, UCSD.
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



function varargout = winPACT_optimize(varargin)
% WINPACT_OPTIMIZE MATLAB code for winPACT_optimize.fig
%      WINPACT_OPTIMIZE, by itself, creates a new WINPACT_OPTIMIZE or raises the existing
%      singleton*.
%
%      H = WINPACT_OPTIMIZE returns the handle to a new WINPACT_OPTIMIZE or the handle to
%      the existing singleton*.
%
%      WINPACT_OPTIMIZE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WINPACT_OPTIMIZE.M with the given input arguments.
%
%      WINPACT_OPTIMIZE('Property','Value',...) creates a new WINPACT_OPTIMIZE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before winPACT_optimize_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to winPACT_optimize_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help winPACT_optimize

% Last Modified by GUIDE v2.5 17-Aug-2018 15:31:45

% 05/08/2018 Makoto. Modified.
% 

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @winPACT_optimize_OpeningFcn, ...
                   'gui_OutputFcn',  @winPACT_optimize_OutputFcn, ...
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


% --- Executes just before winPACT_optimize is made visible.
function winPACT_optimize_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to winPACT_optimize (see VARARGIN)

% Choose default command line output for winPACT_optimize
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes winPACT_optimize wait for user response (see UIRESUME)
% uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.
function varargout = winPACT_optimize_OutputFcn(hObject, eventdata, handles) 
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


% --- Executes on button press in computeSnrLimitPushbutton.
function computeSnrLimitPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to computeSnrLimitPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Number of windows to generate. Half with PAC, half with noise.
numWindows = 40; 

% Generate simulated PAC data.
phaseLowFreqHz  = str2num(get(handles.phaseFreqEdit, 'string'));
ampHighFreqHz   = str2num(get(handles.ampFreqEdit, 'string'));
dataLengthInSec = str2num(get(handles.windowSizeEdit, 'string'))*numWindows;
coloredNoise    = get(handles.useColoredNoisePopupmenu, 'value');
switch coloredNoise
    case 1
        coloredNoise = 1;
    case 2
        coloredNoise = 0;
end

% Set 100 noise levels from -10 to 10 dB SNR (Power ratio: 3dB = x2, 6dB = x4, 10dB = x10).
noiseLevelList = logspace(log10(1/sqrt(20)), log10(sqrt(5)), 100);

% Generate one simulated data.
[lfoHfoNoise, SNR] = Oezkurt2011_synthesize_pac_modified(phaseLowFreqHz, ampHighFreqHz, 1, dataLengthInSec, coloredNoise);

% Generate data with 50 different SNRs
pacData            = sum(lfoHfoNoise([1 2],:), 1);
pacDataHalfPresent = pacData;
pacDataHalfPresent(1:ceil(length(pacData)/2)) = 0;
snrList          = zeros(length(noiseLevelList), 1);
snrLimitTestData = zeros(length(noiseLevelList), size(lfoHfoNoise,2));
for snrIdx = 1:length(noiseLevelList)
    currentNoiseNoise          = lfoHfoNoise(3,:)*noiseLevelList(snrIdx);
    snrList(snrIdx)            = 10*log10((pacData*pacData') / (currentNoiseNoise*currentNoiseNoise'));
    snrLimitTestData(snrIdx,:) = pacDataHalfPresent + currentNoiseNoise;
end

% Export it to EEGLAB.
EEG = eeg_emptyset;
EEG.setname = 'SNR Limit test';
EEG.subject = 'SNR Limit test';
EEG.nbchan  = size(snrLimitTestData,1);
EEG.trials  = 1;
EEG.pnts    = size(snrLimitTestData,2);
EEG.srate   = 1000;
EEG.xmin    = 0.001;
EEG.xmax    = EEG.pnts/1000;
EEG.times   = 0.001:0.001:EEG.xmax;
EEG.data    = snrLimitTestData;

% Set frequency bandpass edges.
lowFreqForPhase = [1 phaseLowFreqHz*2];
highFreqForAmp  = [ampHighFreqHz-phaseLowFreqHz ampHighFreqHz+phaseLowFreqHz];

% Preserve the result container.
windowLengthS = str2num(get(handles.windowSizeEdit, 'String'));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Obtain window center time %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eventTimesSec = windowLengthS/2:windowLengthS:dataLengthInSec; % No overlap.

% Make sure eventTimeSec is a column vector.
if size(eventTimesSec,2)>size(eventTimesSec,1)
    eventTimesSec = eventTimesSec';
end

% Check the window coverage.
eventTimesSec = sort(eventTimesSec);
currentWindowFrameEdges      = [(eventTimesSec-windowLengthS/2) (eventTimesSec+windowLengthS/2)]*EEG.srate;
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
numIterations         = 200;
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



%%%%%%%%%%%%%%%%%%%%%%%
%%% Perform t-test. %%%
%%%%%%%%%%%%%%%%%%%%%%%
measureList = zeros(length(noiseLevelList),3);
pvalList    = zeros(length(noiseLevelList),3);
for snrIdx = 1:length(noiseLevelList)
    
    [~,P1,~,STATS1] = ttest2(canoltysMIAllChan(snrIdx,20:40), canoltysMIAllChan(snrIdx,1:20));
    [~,P2,~,STATS2] = ttest2(MInormAllChan(snrIdx,20:40),     MInormAllChan(snrIdx,1:20));
    [~,P3,~,STATS3] = ttest2(klDistAllChan(snrIdx,20:40),     klDistAllChan(snrIdx,1:20));
    
    measureList(snrIdx,:) = [STATS1.tstat STATS2.tstat STATS3.tstat];
    pvalList(   snrIdx,:) = [P1 P2 P3];
end

% Apply FDR correction.
fdrPvalList(:,1) = fdr(pvalList(:,1));
fdrPvalList(:,2) = fdr(pvalList(:,2));
fdrPvalList(:,3) = fdr(pvalList(:,3));

% Store the results under EEG.etc.winPACT
EEG.etc.winPACT.optimizeMiNormmiKld = measureList;
EEG.etc.winPACT.optimizeFdrPvalList = fdrPvalList;
EEG.etc.winPACT.optimizeSnrList     = snrList;
assignin('base', 'EEG', EEG);

% Plot the initial result.
plotMeasurePopupmenu_Callback(hObject, eventdata, handles)

% Done.
disp('Done. Try other measures to plot.')



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


% --- Executes on selection change in useColoredNoisePopupmenu.
function useColoredNoisePopupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to useColoredNoisePopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns useColoredNoisePopupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from useColoredNoisePopupmenu


% --- Executes during object creation, after setting all properties.
function useColoredNoisePopupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to useColoredNoisePopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function plotMeasurePopupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plotMeasurePopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in plotMeasurePopupmenu.
function plotMeasurePopupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to plotMeasurePopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plotMeasurePopupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plotMeasurePopupmenu

EEG = evalin('base', 'EEG');

if ~isfield(EEG.etc.winPACT, 'optimizeMiNormmiKld')
    error('Compute SNR Limit first.')
end

% Obtain the current measure.
switch get(handles.plotMeasurePopupmenu, 'value')
    case 1
        currentTscore = EEG.etc.winPACT.optimizeMiNormmiKld(:,1);
        currentPval   = EEG.etc.winPACT.optimizeFdrPvalList(:,1);
        
    case 2
        currentTscore = EEG.etc.winPACT.optimizeMiNormmiKld(:,2);
        currentPval   = EEG.etc.winPACT.optimizeFdrPvalList(:,2);
        
    case 3
        currentTscore = EEG.etc.winPACT.optimizeMiNormmiKld(:,3);
        currentPval   = EEG.etc.winPACT.optimizeFdrPvalList(:,3);
end

% Reverse the order.
currentTscore = currentTscore(end:-1:1);
currentPval    = currentPval(end:-1:1);
currentSnr     = EEG.etc.winPACT.optimizeSnrList(end:-1:1);

% Detect the p = 0.01 line.
significanceIdx = find(currentPval < 0.01, 1);

% Plot results.
cla(handles.resultAxes);
plot(handles.resultAxes, currentSnr, currentTscore, '.', 'markersize', 14, 'color', [0 0 1]);
xlim([-10 10])
ylim([min(currentTscore) max(currentTscore)])
xlabel('SNR [dB]')
ylabel('t-score [from two-sample t-test of PAC present vs. absent]')

hold on
textMargin = (max(currentTscore)-min(currentTscore))*0.05;
if any(significanceIdx)
    line([currentSnr(significanceIdx) currentSnr(significanceIdx)], [min(currentTscore) max(currentTscore)], 'linewidth', 2, 'color', [1 0 0]);
    textHandle = text(currentSnr(significanceIdx)+(20*0.02), max(currentTscore)-textMargin, sprintf('p = 0.01 (FDR)\nSNR = %.1f',currentSnr(significanceIdx)), 'color', [1 0 0]);
else
    textHandle = text(currentSnr(1)+(20*0.02), max(currentTscore)-textMargin, sprintf('n.s.'), 'color', [1 0 0]);
end