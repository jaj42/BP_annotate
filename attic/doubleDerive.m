function [ waveformDDPlus, waveformDD, waveformD ] = doubleDerive( waveform )
    waveformD = diff(waveform);
    waveformD = [waveformD NaN];

    waveformDD = diff(waveformD);
    waveformDD = [waveformDD NaN];
    waveformDD = BP_lowpass(waveformDD);

    %Perform the switch for positive and negative first-derivative values
    waveformDDPlus = waveformDD .* (waveformD > 0 & waveformDD > 0);
    waveformDDPlus = waveformDDPlus .^ 2;
end