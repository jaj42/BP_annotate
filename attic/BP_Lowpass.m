function [ filtwaveform ] = BP_Lowpass( waveform )
%LOWPASS Filter input signal (200 Hz assumed)
    b = [1 0 0 0 0 0 -2 0 0 0 0 0 1];
    a = [1 -2 1] * 36; % Multiply by 36 to remove DC gain
    unit = [1 zeros(1, 12)];

    h_l = filter(b, a, unit);

    filtwaveform = conv(waveform, h_l);
    filtwaveform = circshift(filtwaveform, -5, 2);
    filtwaveform(1:5) = NaN;
    filtwaveform = filtwaveform(1 : length(waveform));