function [ footIndex, systolicIndex, notchIndex ] = BP_annotate( waveform, fs, verbose )

newFs = 200;
waveform = BP_resample(waveform, fs);
waveform = BP_Lowpass(waveform);

[ waveformDDPlus, waveformDD ] = doubleDerive(waveform );

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
systolicIndex = FixIndex(footIndex + integwinsize, waveform, Down, integwinsize);
[ dicroticIndex, notchIndex ] = getDicroticIndex( waveformDD, fs, footIndex, systolicIndex );  

if verbose
    Colors = get(gca, 'ColorOrder');
    time = (0: length(waveform) - 1) ./ fs;
    newTime = (0: length(waveformDDPlus) - 1) ./ newFs;
    figure
    axs(1) = subplot(4, 1, 1); hold on;
    plot(time, waveform);
    plot(time(footIndex), waveform(footIndex),'<','color', Colors(4,:), 'markerfacecolor',Colors(4,:))
    plot(newTime(systolicIndex), waveform(systolicIndex),'^','color', Colors(5,:), 'markerfacecolor',Colors(5,:))
    plot(newTime(notchIndex), waveform(notchIndex),'^','color', Colors(6,:), 'markerfacecolor',Colors(6,:))
    plot(newTime(dicroticIndex), waveform(dicroticIndex),'^','color', Colors(7,:), 'markerfacecolor',Colors(7,:))
    legend({'Waveform','Foot','Systole', 'Notch', 'Dicrotic Peak'},'box','off')


    axs(2) = subplot(4, 1, 2);
    hold on;
    plot(newTime, waveformDDPlus./max(waveformDDPlus));
    plot(newTime, BP_integral./max(BP_integral));
    plot(newTime, threshold./max(threshold));
    plot(newTime, zoneOfInterest);
    legend({'2nd Derivative','Integral', 'Threshold', 'ZOI'}, 'box','off');

    axs(3) = subplot(4, 1, 3);

    axs(4) = subplot(4, 1, 4);

    linkaxes(axs, 'x')
end