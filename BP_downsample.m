function [ newWaveform, newFs ] = BP_downsample( waveform, fs )
% Downsamples to 200Hz if fs is higher

if fs > 200
    oldTime = (0: length(waveform) -1) ./ fs;
    newTime = (oldTime(1): 1/200: oldTime(end));
    newWaveform = interp1(oldTime, waveform, newTime, 'pchip');%Interpolation for downsampling
    newWaveform = newWaveform(:);%Column
    
    %Low-pass filter to remove noise (assumes fs = 200)
    b = [1 0 0 0 0 0 -2 0 0 0 0 0 1];
    a = [1 -2 1];
    h_l = filter(b,a,[1 zeros(1,12)]); 
    newWaveform = conv (newWaveform ,h_l);
    newWaveform = [newWaveform(13).*ones(7,1); newWaveform(13:end-13); newWaveform(end-13).*ones(20, 1)];
    newFs = 200;
else
    newWaveform = waveform;
    newFs = fs;
end

end
