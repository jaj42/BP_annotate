function [ footIndex, systolicIndex, notchIndex ] = BP_annotate( waveform, fs, verbose )
%[ output_args ] = BP_annotate( waveform, fs )
    %footIndex = 0;
    systolicIndex = 0;
    notchIndex = 0;

    waveform = waveform(:); %Column

    % Downsample the input signal to filter noise.
    [waveformDS, newFs] = BP_downsample(waveform, fs);

    waveformDDPlus = doubleDerive( waveformDS, fs );
    
    integrWinSize = floor(newFs / 4);
    threshWinSize = floor(newFs * 3);

    % Calculate a rolling sum of the 2nd derivative
    integwindow = rollingWindow(waveformDDPlus, integrWinSize);
    integral = winsum(integwindow);
    % Center the integral
    integral = circshift(integral, -floor(integrWinSize / 2), 2);

    threswindow = rollingWindow(integral, threshWinSize);
    threshold = winquant(threswindow, .7);
    
    zois = integral > threshold;
    
    footIndex = getFootIndex( waveformDDPlus, zois );

    if verbose
        time    = (0: length(waveformDS) - 1) ./ newFs;
        figure
        axs(1) = subplot(4, 1, 1);
        plot(time, waveformDS);
        
        axs(2) = subplot(4, 1, 2);
        plot(time, waveformDDPlus);
        
        axs(3) = subplot(4, 1, 3);
        plot(time, integral);
        hold on;
        plot(time, threshold);
        legend('Integral', 'Threshold');

        axs(4) = subplot(4, 1, 4);
        plot(time, zois);
        legend('ZOIs');
        
        linkaxes(axs, 'x')
    end
end