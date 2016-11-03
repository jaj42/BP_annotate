function [ dicroticIndex, notchIndex ] = getDicroticIndex( waveformDD, fs, footIndex, systolicIndex )
%[ dicroticIndex ] = getDicroticIndex( waveform, fs, footIndex, systolicIndex )
%
%
%%
    RR = median(footIndex(2:end) - footIndex(1:end-1)) ./ fs;% This assumes steady heartrate

    Down = 1;
    minWavelength = round(RR/4 .* fs);
    notchIndex = FixIndex(systolicIndex + minWavelength, waveformDD, ~Down, minWavelength);
    dicroticIndex = FixIndex(notchIndex + round(0.5*minWavelength), waveformDD, Down, minWavelength);

end

