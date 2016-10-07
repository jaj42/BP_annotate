function [ footIndex, systolicIndex, notchIndex ] = BP_annotate( waveform, fs, verbose )
%[ output_args ] = BP_annotate( waveform, fs )
%
%
%%
    footIndex = 0;
    systolicIndex = 0;
    notchIndex = 0;
    
    [ waveformDDPlus, newFs ] = doubleDerive( waveform, fs );
    
    integwinsize = floor(newFs / 4);
    threswinsize = floor(newFs * 3);

    integwindow = window(waveformDDPlus, integwinsize);
    integral = winsum(integwindow);
    % Center the integral
    integral = circshift(integral, -floor(integwinsize / 2), 2);

    threswindow = window(integral, threswinsize);
    threshold = winquant(threswindow, .7);
    
    zois = integral > threshold;
    
    if verbose
        time = (0: length(waveform) - 1) ./ fs;
        newTime = (0: length(waveformDDPlus) - 1) ./ newFs;
        figure
        axs(1) = subplot(4, 1, 1);
        plot(time, waveform);
        
        axs(2) = subplot(4, 1, 2);
        plot(newTime, waveformDDPlus);
        
        axs(3) = subplot(4, 1, 3);
        plot(newTime, integral);
        hold on;
        plot(newTime, threshold);
        legend('Integral', 'Threshold');

        axs(4) = subplot(4, 1, 4);
        plot(newTime, zois);
        legend('ZOIs');
        
        linkaxes(axs, 'x')
    end
end