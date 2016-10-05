function [ newWaveform, newFs ] = BP_downsample( waveform, fs )
%[ waveform, fs ] = BP_downsample( waveform, fs )
%
%Downsamples to 200Hz if fs is higher
%%
if fs > 200
    oldTime = (0: length(waveform) -1) ./ fs;
    newTime = (oldTime(1): 1/200: oldTime(end));
    newWaveform = interp1(oldTime, waveform, newTime, 'pchip');
    newWaveform = newWaveform(:);%Column
    
    b = [1 0 0 0 0 0 -2 0 0 0 0 0 1];
    a = [1 -2 1];
    h_l = filter(b,a,[1 zeros(1,12)]); 
    newWaveform = conv (newWaveform ,h_l);

    newFs = 200;
else
    newWaveform = waveform;
    newFs = fs;
end

end

