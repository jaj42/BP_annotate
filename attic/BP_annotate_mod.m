function [ footIndex, systolicIndex, notchIndex ] = BP_annotate( waveform, fs, verbose )
newFs = 200;
origwaveform = waveform;

[bpwaveform, time, origtime] = BP_resample(waveform, fs);
bpwaveform = BP_Lowpass(bpwaveform);

[ waveformDDPlus, waveformDD ] = doubleDerive(bpwaveform );

integwinsize = floor(newFs / 4);
threswinsize = floor(newFs * 3);

% Moving sum to increase SNR
integralWindow = rollingWindow(waveformDDPlus, integwinsize);
BP_integral = winsum(integralWindow);
% Center the integral
BP_integral = circshift(BP_integral, -floor(integwinsize / 2), 2);

thresholdWindow = rollingWindow(BP_integral, threswinsize);
threshold = winquant(thresholdWindow, .7);

[ zoneOfInterest ] = getZoneOfInterest( BP_integral, threshold );
footIndex = getFootIndex( waveformDDPlus, zoneOfInterest );

Down = 0;
systolicIndex = FixIndex(footIndex + integwinsize, bpwaveform, Down, integwinsize);
[ dicroticIndex, notchIndex ] = getDicroticIndex( waveformDD, fs, footIndex, systolicIndex );  

if verbose
    Colors = get(gca, 'ColorOrder');
    axs(1) = subplot(2, 1, 1);
    hold on;
    plot(time, bpwaveform);
    plot(origtime, origwaveform);
    plot(time(footIndex), bpwaveform(footIndex),         '<', 'color', Colors(4,:), 'markerfacecolor', Colors(4,:))
    plot(time(systolicIndex), bpwaveform(systolicIndex), '^', 'color', Colors(5,:), 'markerfacecolor', Colors(5,:))
    plot(time(notchIndex), bpwaveform(notchIndex),       '^', 'color', Colors(6,:), 'markerfacecolor', Colors(6,:))
    plot(time(dicroticIndex), bpwaveform(dicroticIndex), '^', 'color', Colors(7,:), 'markerfacecolor', Colors(7,:))
    legend({'Filtered','Waveform','Foot','Systole', 'Notch', 'Dicrotic Peak'},'box','off')

    axs(2) = subplot(2, 1, 2);
    hold on;
    plot(time, waveformDDPlus);
    plot(time, BP_integral);
    plot(time, threshold);
    plot(time, zoneOfInterest .* .1);
    legend({'2nd Derivative','Integral', 'Threshold', 'ZOI'}, 'box','off');

    linkaxes(axs, 'x')
end