% 08/06/2018 Makoto. Created. Prayer for Hiroshima.

thisMfilePath = which('pacSimulation');
filesepIdx = strfind(thisMfilePath, '/');
addpath([thisMfilePath(1:filesepIdx(end)) 'external/coloredNoiseByHristoZhivomirov'])

% 180 s data, assuming sampling rate == 1000Hz.
samplingRate = 1000;
channelData = pinknoise(180*samplingRate);
channelData = repmat(channelData, [3 1]);

% Generate 3Hz LFO sine wave.
lfoFreq = 3;
lfoOneCycle = sin(0:pi*2/round(samplingRate/lfoFreq):2*pi);
lfoOneCycle = lfoOneCycle(1:end-1);
lfoTenSec   = repmat(lfoOneCycle, [1 round(10*samplingRate/length(lfoOneCycle))]);

% Generate 100 Hz HFO sine wave.
hfoFreq = 100;
hfoOneCycle = sin(0:pi*2/round(samplingRate/hfoFreq):2*pi);
hfoOneCycle = hfoOneCycle(1:end-1);
hfoTenSec   = repmat(hfoOneCycle, [1 round(10*samplingRate/length(hfoOneCycle))]);

    % % Normalize using Pink noise curve i.e., half the power per octave.
    % % Kramer (2008) and Ozkurt and Schnitzler (2011) did not do this.
    % octabeDifference   = log(hfoFreq)/log(lfoFreq);
    % normalizationCoeff = 1/2^octabeDifference;
    % hfoTenSec          = hfoTenSec*normalizationCoeff;

% Amplitude-modulate the HFO so that PAC happens centered at -0.5pi.
centerAngleRad = -1;
phaseMargin    = 2*pi/8; % +/-45 degrees.
instPhaseLFO   = angle(hilbert(lfoTenSec));
targetPhaseIdx = (instPhaseLFO > wrapTo2Pi(centerAngleRad*pi-phaseMargin)) | (instPhaseLFO < centerAngleRad*pi+phaseMargin);
hfoMasked      = hfoTenSec(1:length(targetPhaseIdx)).*targetPhaseIdx;
pacSignal      = lfoTenSec + hfoMasked;
pacSignalNorm  = pacSignal/std(pacSignal);
lfoTenSecNorm  = lfoTenSec/std(lfoTenSec);

%{
figure
plot(pacSignal)
hold on
plot(instPhaseLFO, 'r')
%}
    
channelDataNorm = bsxfun(@rdivide, channelData, std(channelData,0,2));

% Set the SNR. Use 10*log10(x)
%snFactor = 5; % 10*log10(5) = SNR 7dB.
snFactor = 0.1; % 10*log10(10) = SNR -10dB.

    
% Implement the signal.
repetitionFactor = 1;

signalTime = 100*samplingRate+1:100*samplingRate+length(pacSignal)*repetitionFactor;

channelDataNorm(1,signalTime) = channelDataNorm(1,signalTime) + repmat(pacSignal,     [1 repetitionFactor])*snFactor; % Signal + noise
channelDataNorm(2,signalTime) = channelDataNorm(2,signalTime) + repmat(lfoTenSecNorm, [1 repetitionFactor])*snFactor; % LFO + noise
channelDataNorm(3,signalTime) = channelDataNorm(3,signalTime);                                                        % Noise

% Import it to EEGLAB.
EEG = pop_importdata('dataformat', 'array', 'nbchan', 0, 'data','channelDataNorm', 'srate', 1000,'pnts', 0, 'xmin', 0);
EEG = pop_eegfiltnew(EEG, [],1,3300,1,[],0);

[spectra,freqs,speccomp,contrib,specstd] = ...
                     spectopo(EEG.data, EEG.pnts, EEG.srate, 'freqrange', [2 500], 'freqfac', 8);

eeglab redraw