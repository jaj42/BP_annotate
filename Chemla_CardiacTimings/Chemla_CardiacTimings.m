function [ heartPeriod, LVET, diastolicTime, footIndex, systolicIndex, notchIndex, dicroticIndex, time, bpwaveform, missingIndices ] = Chemla_CardiacTimings( inWaveform, inFs, varargin )
%%[ heartPeriod, LVET, diastolicTime, time, bpwaveform ] = Chemla_CardiacTimings( inWaveform, inFs, verbose )
% An implementation of heart period, left ventricular ejection time (LVET),
% and diastolic time as described by Chemla et al. in:
%
% Chemla, Denis, et al. "Short-term variability of pulse pressure and 
% systolic and diastolic time in heart transplant recipients." American 
% Journal of Physiology-Heart and Circulatory Physiology 279.1 (2000): 
% H122-H129. 
%
% Depends on the BP_annotate software from http://fr.mathworks.com/matlabcentral/fileexchange/60172-bp-annotate
%% Inputs
% inWaveform : countinuous arterial blood pressure time-series
% inFs : sampling frequency (Hz) of the time-series
% verbose : boolean, should be true if figures are wanted
%
%% Outputs
% heartPeriod : the heart period (ms)
% LVET : left ventricular ejection time (ms)
% diastolicTime : diastolic time (ms)
% time : time vector (s) of the resampled 200 Hz time-series
% bpwaveform : resampled filtered 200 Hz time-series
%
%% Authors
% Alexandre Laurin, PhD, ?cole Polytechnique, France, alexandre.laurin@inria.fr
% Jona Joachim, MD, H?pital Lariboisi?re, France, jona.joachim@inria.fr
%%

verbose = 1;
Units = 'unkown';
isClean = 0;
if nargin >= 3
    verbose = varargin{1};
end
if nargin >=4
    Units = varargin{2};
end
if nargin >=5
    isClean = varargin{3};
end
% All blood pressure waveforms are resampled at 200Hz. Annotation indices
% reflect this choice.
[ footIndex, systolicIndex, notchIndex, dicroticIndex, time, bpwaveform ] = BP_annotate( inWaveform, inFs, verbose, Units, isClean );
heartPeriod = 1000 .* (time(footIndex(2:end)) - time(footIndex(1:end - 1)));
LVET = 1000 .* (time(notchIndex(1:end-1)) - time(footIndex(1:end-1)));
diastolicTime = heartPeriod - LVET;

% figure(69); clf; plot(LVET,'+-'); hold on; plot(diastolicTime,'x-')
if strcmp(Units, 'mmHg')
    % clean up the indices a little when aberrations occur
    
    medianFoot = median(bpwaveform(footIndex));
    medianNotch = median(bpwaveform(notchIndex));
    medianSystolic = median(bpwaveform(systolicIndex));
    footDist = abs(bpwaveform(footIndex) - medianFoot);
    notchDist = abs(bpwaveform(notchIndex) - medianNotch);
    maxDist = max([footDist; notchDist]);
    missingIndices = find(maxDist > (medianSystolic - medianFoot)/4);
    missingIndices(missingIndices == length(footIndex)) = [];
    
%     missingIndices = find(LVET > diastolicTime);
%     if (max(diastolicTime - LVET) > 0)
        heartPeriod(missingIndices ) = [];
        LVET( missingIndices ) = [];
        diastolicTime( missingIndices ) = [];
        footIndex( missingIndices ) = [];
        systolicIndex( missingIndices ) = [];
        notchIndex( missingIndices ) = [];
        dicroticIndex( missingIndices ) = [];
%     end
else
    missingIndices = [];
end

if verbose
    figure;
    hold on;
    plot(time(footIndex(1:end - 1)), heartPeriod,'*-');
    plot(time(footIndex(1:end - 1)), LVET,'+-');
    plot(time(footIndex(1:end - 1)), diastolicTime,'x-');
    xlabel('time (s)')
    ylabel('cardiac timing (ms)')
    legend({'heart period','LVET','diastolic time'},'box','off')
end
