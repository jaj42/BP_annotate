function [ newWaveform, newx, oldx ] = BP_resample( waveform, fs )
%BP_RESAMPLE Resample to 200 Hz
newfs = 200;

duration = length(waveform) / fs;

oldx = linspace(0, duration, length(waveform));
newx = linspace(0, duration, newfs * duration);

newWaveform = interp1(oldx, waveform, newx, 'pchip');