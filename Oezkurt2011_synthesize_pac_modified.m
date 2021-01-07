% Oezkurt2011_synthesize_pac_modified()--A function to generate EEGLAB data
%                                        structure to generate phase-amplitude
%                                        coupled simulation data.
%                                        Data sampling rate is fixed to be
%                                        1000 Hz.
%
% Output:
%         lfoHfoNoise   -- [3 x dataLengthInSec*1000] The first row is LFO, 
%                          the second row HFO, and the third row noise.
%                          Take sum(lfoHfoNoise) to obtain PAC data.
%         SNR           -- Signal-to-noise ratio. 10*log10(x)
%
%%%%%%%%%%%%%%%%%%%%%%%%%% From original function. 
%
% adapted from Kramer et al. (2008), Jrn. Nrsc. Methds.  
%
% this functions simulates a signal containing phase-amplitude cooupling (PAC)
% between the phase of a low frequency (15-16 Hz) and amplitude of a
% high-freqency band (60-90 Hz) for a chosen noise level. Hamming tapered
% high frequency signal is added at each cycle of low frequency component.
% sampling frequency is chosen as 1000 Hz.
%
% additionally 5 separate sinusoids (having NO coupling) to evaluate PAC
% estimation methods' discrimination capability
%
% noise lev : parameter describing the noise power
% snr : signal-to-noise ratio
% s_final: the synthesized signal
%
% For details of the simulation, please see (Ozkurt and Schnitzler, 2011)

% History
% 08/15/2018 Makoto. Copied and modified.

function [lfoHfoNoise, SNR] = Oezkurt2011_synthesize_pac_modified(phaseLowFreqHz, ampHighFreqHz, noiseLevel, dataLengthInSec, coloredNoise)

% Add path to the colored noise generator by Hristo Zhivomirov.
thisMfilePath = which('Oezkurt2011_synthesize_pac_modified');
filesepIdx = strfind(thisMfilePath, filesep);
addpath([thisMfilePath(1:filesepIdx(end)) 'external' filesep 'coloredNoiseByHristoZhivomirov'])

% Set the amount of fluctuation in LFO.
lfoFluctuation = 0.1; % Unit: Octave.

% Generate HFO.
continuousHfo = randn(1,dataLengthInSec*1000);                                                           % Generate white noise..
continuousHfo = eegfiltnew(continuousHfo, 1000, ampHighFreqHz-phaseLowFreqHz*2, ampHighFreqHz+phaseLowFreqHz*2); % Band-pass fitler the white noise.

% Compute the number of cycles needed.
deltaTime = 0.001; % Delta = 1 ms.
totalDataLengthInMs = dataLengthInSec*1000;
oneCycleLengthInMs  = floor(1000/phaseLowFreqHz);
numCyclesNeeded     = ceil(totalDataLengthInMs/oneCycleLengthInMs);
numCyclesNeeded     = ceil(numCyclesNeeded*1.1);

% Generate simulated PAC cycle by cycle.
simulatedPac = cell(numCyclesNeeded,1);
simulatedPac = cell(numCyclesNeeded,1);
for cycleIdx = 1:numCyclesNeeded
    
    % Generate sinunoidal LFO with center freq +/- lfoFluctuation [Hz].
    currentFreq = phaseLowFreqHz + rand()*lfoFluctuation*phaseLowFreqHz; % Center freq + fluctuation.
    currentLfo  = cos(2.0*pi*(0:deltaTime*currentFreq:1-deltaTime*currentFreq));
    
    % Identify the trough index with trough phase +/- 1/2 pi fluctuation.
    unwrappedPhase      = unwrap(angle(hilbert(currentLfo)));
    wholeTroughIdx      = find(((1-1/2)*pi < unwrappedPhase) & (unwrappedPhase < (1+1/2)*pi));
    quarterTroughLength = floor(length(wholeTroughIdx)/4);
    threeQuarterTroughLength = length(wholeTroughIdx)-quarterTroughLength;
    beginningIdx     = randi(quarterTroughLength);
    currentTroughIdx = wholeTroughIdx(beginningIdx:beginningIdx+threeQuarterTroughLength-1);
    
    % Obtain the HFO for the current cycle.
    currentHfoOnset   = randi(length(continuousHfo)-length(currentTroughIdx));
    currentHfo        = continuousHfo(currentHfoOnset:currentHfoOnset+length(currentTroughIdx)-1);
    currentHfoHanning = hanning(length(currentHfo))'.*currentHfo; % Hanning tapered.
    
    % Combine LFO and HFO.
    singlePacCycle = currentLfo;
    singlePacCycle(currentTroughIdx) = singlePacCycle(currentTroughIdx) + currentHfoHanning;
    
    % Store LFO and Pac.
    lfo{cycleIdx}          = currentLfo;
    simulatedPac{cycleIdx} = singlePacCycle;
end

% Concatenate the cycle-by-cycle PAC singal.
lfo          = cat(2, lfo{:});
simulatedPac = cat(2, simulatedPac{:});
hfo          = simulatedPac-lfo;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set the noise level. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if coloredNoise == 0
    generatedNoise = noiseLevel*randn(1,length(simulatedPac)); % White noise.
else
    generatedNoise = noiseLevel*pinknoise(length(simulatedPac)); % Pink noise.
end
SNR = 10*log10((simulatedPac*simulatedPac') / (generatedNoise*generatedNoise'));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Concatenate LFO, HFO, and noise. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lfoHfoNoise = cat(1, lfo, hfo, generatedNoise);
lfoHfoNoise = lfoHfoNoise(:,1:dataLengthInSec*1000);



% Support Andreas's filter function.
function filteredData = eegfiltnew(data, srate, locutoff, hicutoff)
tmpFiltData.data   = data;
tmpFiltData.srate  = srate;
tmpFiltData.trials = 1;
tmpFiltData.event  = [];
tmpFiltData.pnts   = length(data);
tmpFiltData_done   = pop_eegfiltnew(tmpFiltData, locutoff, hicutoff);
filteredData       = tmpFiltData_done.data;