function [ waveformDDPlus ] = doubleDerive( waveform, fs )
%[ waveformDD ] = doubleDerive( waveform, fs )
%
%Returns the second derivative of the waveform if the first derivative is
%positive; otherwise returns
%%
    %Perform the first derivative
    waveformD = diff(waveform);
    waveformD = waveformD .* fs;
    waveformD = [waveformD(1); waveformD];
    
    %Perform the second derivative
    waveformDD = diff(waveformD);
    waveformDD = waveformDD.*fs;
    waveformDD = [0; waveformDD];
    
    %Perform the switch for positive and negative first-derivative values
    waveformDDPlus = zeros( size(waveformDD) );
    waveformDDPlus(waveformD > 0) = waveformDD(waveformD > 0);
    waveformDDPlus = waveformDDPlus .^ 2;
    waveformDDPlus(1: 2) = 0;
end

