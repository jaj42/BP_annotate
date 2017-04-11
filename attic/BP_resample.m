function [ newWaveform, newx, oldx ] = BP_resample( waveform, origFs )
    %BP_RESAMPLE Resample to 200 Hz
    Fs = 200;
    duration = length(waveform) / origFs;
    oldx = linspace(0, duration, length(waveform));
    newx = linspace(0, duration, Fs * duration);

    newWaveform = interp1(oldx, waveform, newx, 'pchip');