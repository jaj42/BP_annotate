function [ waveformDD ] = doubleDerive( waveform, fs )
%[ waveformDD ] = doubleDerive( waveform, fs )

waveformDD = diff(diff(waveform));
waveformDD = waveformDD.*fs;

end

