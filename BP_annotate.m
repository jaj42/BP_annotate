function [ footIndex, systolicIndex, notchIndex ] = BP_annotate( waveform, fs, verbose )
%[ output_args ] = BP_annotate( waveform, fs )
%
%
%%
    footIndex = 0;
    systolicIndex = 0;
    notchIndex = 0;
    
    [ waveformDDPlus, newFs ] = doubleDerive( waveform, fs );
    


    if verbose
        time = (0: length(waveform) - 1) ./ fs;
        newTime = (0: length(waveformDDPlus) - 1) ./ newFs;
        figure
        axs(1) = subplot(3, 1, 1);
        plot(time, waveform);
        
        axs(2) = subplot(3, 1, 2);
        plot(newTime, waveformDDPlus)
        
        linkaxes(axs, 'x')
    end
end

