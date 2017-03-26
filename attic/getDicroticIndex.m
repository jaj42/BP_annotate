function [ dicroticIndex, notchIndex ] = getDicroticIndex( waveformDD, waveformD, bpwaveform, footIndex, systolicIndex )
    Fs = 200;
    RR = median(footIndex(2:end) - footIndex(1:end-1)) ./ Fs;% This assumes steady heartrate
    Down = 1; Up = ~Down;
    minWavelength = round(RR/5 .* Fs);
    notchIndex = FixIndex(systolicIndex + minWavelength, waveformDD, Up, minWavelength);
    dicroticIndex = FixIndex(notchIndex + round(0.5*minWavelength), waveformDD, Down, minWavelength);
    systolicIndex = systolicIndex(1:length(dicroticIndex));
    
    % if a local minimum and maximum exist, move the dicrotic indices to
    % these
    for i = 1 : length(systolicIndex)
        Start = systolicIndex(i) + round(minWavelength/4);
        End = min([dicroticIndex(i) + round(minWavelength*5/3), length(waveformD)]);
        ZOI = waveformD(Start : End);
        ZOI = ZOI(2:end).*ZOI(1:end-1);
        extrema = find(ZOI<0);
        if length(extrema) >=2
            notchIndex(i) = FixIndex(notchIndex(i), bpwaveform, Down, 4);
            dicroticIndex(i) = FixIndex(min(notchIndex(i) + round(0.5*minWavelength), length(bpwaveform)), bpwaveform, Up, 4);
        end
    end
end

