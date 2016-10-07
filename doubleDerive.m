function [ waveformDDPlus, waveformDD ] = doubleDerive( waveform, fs )
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
     %Low-pass filter to remove noise (assumes fs = 200)
    b = [1 0 0 0 0 0 -2 0 0 0 0 0 1];
    a = [1 -2 1];
    h_l = filter(b,a,[1 zeros(1,12)]); 
    waveformDD = conv (waveformDD ,h_l);
    waveformDD = [waveformDD(13).*ones(7,1); waveformDD(13:end-13); waveformDD(end-13).*ones(20, 1)];
    
    %Perform the switch for positive and negative first-derivative values
    waveformDDPlus = zeros( size(waveformDD) );
    waveformDDPlus(waveformD > 0) = waveformDD(waveformD > 0);
    waveformDDPlus = waveformDDPlus .^ 2;
    waveformDDPlus(1: 2) = 0;
end
