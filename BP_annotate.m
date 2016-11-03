function [ footIndex, systolicIndex, notchIndex, dicroticIndex ] = BP_annotate( inWaveform, inFs, verbose )
    
    global Fs;
    Fs = 200;
    integwinsize = floor(Fs / 4);
    threswinsize = floor(Fs * 3);
    
    disp('Annotating blood pressure...');
    inWaveform = inWaveform(:);
    [bpwaveform, time, origtime] = BP_resample(inWaveform, inFs);
    bpwaveform = BP_lowpass(bpwaveform);
    [ waveformDDPlus, waveformDD ] = doubleDerive(bpwaveform );
      
    %Deal with very large data sets
    sizeLimit = 3*10^5;
    if length(bpwaveform) > sizeLimit
        %Deal with very large data sets
        fprintf('     Data set exceeds %d size limit, performing sub-windowing.',sizeLimit)
        BP_integral = zeros(size(bpwaveform));
        threshold = zeros(size(bpwaveform));
        
        numSubParts = ceil(length(bpwaveform)/sizeLimit);
        overlap = round( ( numSubParts * sizeLimit - length(bpwaveform) ) / (numSubParts - 1) );
%         subParts = struct;
%         subVerbose = 0;
%         footIndex = [];
%         systolicIndex = [];
%         notchIndex = [];
%         dicroticIndex = [];
        for i = 1 : numSubParts
            fprintf('.')
            Start = ( (i-1) * sizeLimit ) - (i-1) * overlap + 1;
            End = min( Start + sizeLimit - 1, length(bpwaveform));
            subIntegralWindow = rollingWindow(waveformDDPlus(Start : End), integwinsize);
            subIntegral = winsum(subIntegralWindow);
            subIntegral = circshift(subIntegral, -floor(integwinsize / 2), 2);
            
            subThresholdWindow = rollingWindow(subIntegral, threswinsize);
            subThreshold = winmean(subThresholdWindow, 1.5);
        
        
            BP_integral(Start + overlap : End) = subIntegral(1+overlap : end);
            threshold(Start + overlap : End) = subThreshold(1+overlap : end);
            if i > 1
                BP_integral(Start : Start + overlap) = mean([BP_integral(Start : Start + overlap); subIntegral(1 : 1+ overlap)]);
                threshold(Start : Start + overlap) = mean([threshold(Start : Start + overlap); subThreshold(1 : 1+ overlap)]);
            end
%             [ subFootIndex, subSystolicIndex, subNotchIndex, subDicroticIndex ] = BP_annotate( bpwaveform(Start : End), Fs, subVerbose );
%             footIndex = unique([footIndex, subFootIndex + Start - 1]);
%             systolicIndex = unique([systolicIndex, subSystolicIndex + Start - 1]);
%             notchIndex = unique([notchIndex, subNotchIndex + Start - 1]);
%             dicroticIndex = unique([dicroticIndex, subDicroticIndex + Start - 1]);
        end
        fprintf('\n')
    else

        % Moving sum to increase SNR
        integralWindow = rollingWindow(waveformDDPlus, integwinsize);
        BP_integral = winsum(integralWindow);
        % Center the integral
        BP_integral = circshift(BP_integral, -floor(integwinsize / 2), 2);

        thresholdWindow = rollingWindow(BP_integral, threswinsize);
        threshold = winmean(thresholdWindow, 1.5);
    end
    
    zoneOfInterest =  BP_integral > threshold ;
    footIndex = getFootIndex( waveformDDPlus, zoneOfInterest );

    Up = 1;
    Down = ~Up;
    systolicIndex = FixIndex(footIndex + floor(integwinsize/2), bpwaveform, Down, floor(integwinsize/2));
    [ dicroticIndex, notchIndex ] = getDicroticIndex( waveformDD, footIndex, systolicIndex );  
        

    if verbose
        figure;
        Colors = get(gca, 'ColorOrder');
        axs(1) = subplot(2, 1, 1);
        hold on;
        plot(time, bpwaveform);
        plot(origtime, inWaveform);
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
    disp('Done.');
end

function [ filtwaveform ] = BP_lowpass( waveform )
    %LOWPASS Filter input signal (200 Hz assumed)
    b = [1 0 0 0 0 0 -2 0 0 0 0 0 1];
    a = [1 -2 1] * 36; % Multiply by 36 to remove DC gain
    unit = [1 zeros(1, 12)];

    h_l = filter(b, a, unit);

    filtwaveform = conv(waveform, h_l);
    filtwaveform = circshift(filtwaveform, -5, 2);
    filtwaveform(1:5) = NaN;
    filtwaveform = filtwaveform(1 : length(waveform));
end

function [ newWaveform, newx, oldx ] = BP_resample( waveform, inFs )
    %BP_RESAMPLE Resample to 200 Hz
    global Fs;

    duration = length(waveform) ./ inFs;

    oldx = linspace(0, duration, length(waveform));
    newx = linspace(0, duration, Fs * duration);

    newWaveform = interp1(oldx, waveform, newx, 'pchip');
end

function [ waveformDDPlus, waveformDD ] = doubleDerive( waveform )
    waveformD = diff(waveform);
    waveformD = [waveformD NaN];

    waveformDD = diff(waveformD);
    waveformDD = [waveformDD NaN];

    %Low-pass filter to remove noise (assumes fs = 200)
    %waveformD  = BP_Lowpass(waveformD);
    %waveformDD = BP_Lowpass(waveformDD);

    %Perform the switch for positive and negative first-derivative values
    waveformDDPlus = waveformDD .* (waveformD > 0 & waveformDD > 0);
    waveformDDPlus = waveformDDPlus .^ 2;
end

function [fixedIndex] = FixIndex(BrokeIndex, BrokeSCG, Down, minWavelength)
    fixedIndex = BrokeIndex;
    for N = 1:length(BrokeIndex)
       if BrokeIndex+2 < length(BrokeSCG)
           OldIndex = BrokeIndex(N);
           if Down
               TempSCG = BrokeSCG(max(OldIndex-round(minWavelength/4),1):min(OldIndex+round(minWavelength/4),end));
               [~, NewIndex] = min(TempSCG);
               NewIndex = NewIndex + OldIndex - round(minWavelength/4) -1;
               while ~(OldIndex==NewIndex)
                   OldIndex = NewIndex;
                   TempSCG = BrokeSCG(max(OldIndex-10,1):min(OldIndex+10,end));
                   [~, NewIndex] = min(TempSCG);
                   NewIndex = NewIndex + max(OldIndex-10,1) -1;
               end
           else
               [~, NewIndex] = max(BrokeSCG(max(OldIndex-round(minWavelength/4),1):min(OldIndex+round(minWavelength/4),length(BrokeSCG))));
               NewIndex = NewIndex + OldIndex - round(minWavelength/4) -1;
               while ~(OldIndex==NewIndex)
                   OldIndex = NewIndex;
                   TempSCG = BrokeSCG(max(OldIndex-10,1):min(OldIndex+10,end));
                   [~, NewIndex] = max(TempSCG);
                   NewIndex = NewIndex + max(OldIndex-10,1) -1;
               end
           end
           Index = NewIndex;
           fixedIndex(N) = Index;
       end
    end
    fixedIndex( fixedIndex > length(BrokeSCG) ) = [];
    fixedIndex( fixedIndex < 1 ) = [];
end

function [ footIndex ] = getFootIndex( waveformDDPlus, zoneOfInterest )
    zoneWall = diff( zoneOfInterest );
    BP_start = find(zoneWall == 1);
    BP_stop = find(zoneWall == -1);

    % Remove leading falling edges
    while BP_stop(1) < BP_start(1)
        BP_stop = BP_stop(2: end);
    end

    nfeet = min(numel(BP_start), numel(BP_stop));
    footIndex = zeros(1, nfeet);
    for i = 1 : nfeet
        [~, footIndex(i)] = max(waveformDDPlus(BP_start(i) : BP_stop(i)));
        footIndex(i) = footIndex(i) + BP_start(i) - 1;
    end
end

function [ dicroticIndex, notchIndex ] = getDicroticIndex( waveformDD, footIndex, systolicIndex )
    global Fs;
    RR = median(footIndex(2:end) - footIndex(1:end-1)) ./ Fs;% This assumes steady heartrate
    Down = 1;
    minWavelength = round(RR/4 .* Fs);
    notchIndex = FixIndex(systolicIndex + minWavelength, waveformDD, ~Down, minWavelength);
    dicroticIndex = FixIndex(notchIndex + round(0.5*minWavelength), waveformDD, Down, minWavelength);
end


% Rolling window functions
function [ rwin ] = rollingWindow( vector, winsize )
    vector = vector(:);
    vecsize = length(vector);
    rwin = NaN(winsize, vecsize);
    for i = 1 : winsize
        tmp = vector(1 : end - i + 1);
        rwin(i, i : end) = tmp;
    end
end

function [ res ] = winmean( rwin, scale )
    shape = size(rwin);
    vecsize = shape(2);
    res = zeros(1, vecsize);
    for i = 1:vecsize
        %res(i) = quantile(rwin(:,i), quant);
        res(i) = scale * mean(rwin(:,i));
    end
end

function [ res ] = winsum( rwin )
    shape = size(rwin);
    vecsize = shape(2);
    res = zeros(1, vecsize);
    for i = 1:vecsize
        res(i) = sum(rwin(:,i));
    end
end