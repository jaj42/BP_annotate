function [ waveformDDPlus, waveformDD ] = doubleDerive( waveform )
%Perform the first derivative
waveformD = diff(waveform);
waveformD = [waveformD NaN];

%Perform the second derivative
waveformDD = diff(waveformD);
waveformDD = [waveformDD NaN];

%Low-pass filter to remove noise (assumes fs = 200)
%waveformD  = BP_Lowpass(waveformD);
%waveformDD = BP_Lowpass(waveformDD);

%Perform the switch for positive and negative first-derivative values
waveformDDPlus = waveformDD .* (waveformD > 0 & waveformDD > 0);
waveformDDPlus = waveformDDPlus .^ 2;
