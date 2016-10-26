function [ waveformDDPlus, waveformDD ] = doubleDerive( waveform )
%Returns the second derivative of the waveform if the first derivative is
%positive; otherwise returns

%Perform the first derivative
waveformD = diff(waveform);
waveformD = [waveformD(1), waveformD];

%Perform the second derivative
waveformDD = diff(waveformD);
waveformDD = [0, waveformDD];

%Low-pass filter to remove noise (assumes fs = 200)
waveformD  = BP_Lowpass(waveformD);
waveformDD = BP_Lowpass(waveformDD);

%Perform the switch for positive and negative first-derivative values
% XXX
waveformDDPlus = zeros( size(waveformDD) );
waveformDDPlus(waveformD > 0 & waveformDD > 0) = waveformDD(waveformD > 0 & waveformDD > 0);
waveformDDPlus = waveformDDPlus .^ 2;
waveformDDPlus(1: 2) = 0;