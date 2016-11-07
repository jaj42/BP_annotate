function [ footIndex, systolicIndex, notchIndex, dicroticIndex, time, bpwaveform ] = BP_annotate( inWaveform, inFs, verbose )
%% function [ footIndex, systolicIndex, notchIndex, dicroticIndex, time, bpwaveform ] = BP_annotate( inWaveform, inFs, verbose )
% Implementation of a feature detection algorithm for arterial blood
% pressure in humans. The foot of the wave, systolic peak, dicrotic notch,
% and dicrotic peaks are identified. The blood pressure time series is
% always resampled at 200 Hz to allow standardisation.
%
% The technique was largely inpired by the derivatives and thresholds 
% described in Pan-Tompkins:
% Pan, Jiapu, and Willis J. Tompkins. "A real-time QRS detection algorithm." 
% IEEE transactions on biomedical engineering 3 (1985): 230-236.

%% Inputs
% inWaveform : countinuous arterial blood pressure time-series
% inFs : sampling frequency of the time-series
% verbose : boolean, should be true if figures are wanted
%% Outputs
% footIndex : index of the foot of each systolic wave
% systolicIndex : index of each sytolic peak
% notchIndex : index of each dicrotic notch
% dicroticIndex : index of each dicrotic peak
% time : time vector of the resampled 200 Hz time-series
% bpwaveform : resampled filtered 200 Hz time-series
%
%% Methods

%% Foot index
% The foot index is defined as the point where the second derivative of 
% the time-series is the highest in each interval where a moving average
% of the second derivative was bigger than a adaptative threshold. This 
% criterion was prefered over others because of its rubustness.
%% Systolic peak index
% The systolic peak is defined as the maximum of the waveform following 
% the foot index, relative to a window of radius 1/8 s around itself. 
%% Dicrotic notch and peak indices
% The dicrotic notch and peak indices were defined as the successive 
% maximum and minimum of the second derivative of the time-series 
% following the systolic peak, relative to a window of radius RR/5 s 
% around itself. (RR is the median heartbeat interval computed form the 
% foot indices). These indices are moved to waveform minima and maxima if
% these existed.
%
%% Authors
% Alexandre Laurin, PhD, ?cole Polytechnique, France, alexandre.laurin@inria.fr
% Jona Joachim, MD, H?pital Lariboisi?re, France, jona.joachim@inria.fr
%%
    global Fs;
    Fs = 200;
    integwinsize = floor(Fs / 4);
    threswinsize = floor(Fs * 3);
    
    disp('Annotating blood pressure...')
    inWaveform = inWaveform(:);
    % resample the time-series to allow standardisation
    [bpwaveform, time, origtime] = BP_resample(inWaveform, inFs);
    % filter
    bpwaveform = BP_lowpass(bpwaveform);
    % derivatives
    [ waveformDDPlus, waveformDD, waveformD ] = doubleDerive(bpwaveform );
    %Deal with very large data sets
    sizeLimit = 3*10^5;
    if length(bpwaveform) > sizeLimit
        
        fprintf('     Data set exceeds %d size limit, performing sub-windowing.',sizeLimit)
        BP_integral = zeros(size(bpwaveform));
        threshold = zeros(size(bpwaveform));
        
        numSubParts = ceil(length(bpwaveform)/sizeLimit);
        overlap = round( ( numSubParts * sizeLimit - length(bpwaveform) ) / (numSubParts - 1) );
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
            else
                BP_integral(Start : Start + overlap) = subIntegral(1 : 1+ overlap);
                threshold(Start : Start + overlap) = subThreshold(1 : 1+ overlap);
            end
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
    
    % each zone of interest corresponds to a heart beat
    zoneOfInterest =  BP_integral > threshold ;
    footIndex = getFootIndex( waveformDDPlus, zoneOfInterest );

    Down = 1; Up = ~Down;
    systolicIndex = FixIndex(footIndex + floor(integwinsize/2), bpwaveform, Up, floor(integwinsize/2));
    [ dicroticIndex, notchIndex ] = getDicroticIndex( waveformDD, waveformD, bpwaveform, footIndex, systolicIndex );  


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
        plot(time, waveformDD);
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

function [ newWaveform, newx, oldx ] = BP_resample( waveform, origFs )
    %BP_RESAMPLE Resample to 200 Hz
    global Fs;
    duration = length(waveform) / origFs;
    oldx = linspace(0, duration, length(waveform));
    newx = linspace(0, duration, Fs * duration);

    newWaveform = interp1(oldx, waveform, newx, 'pchip');
end

function [ waveformDDPlus, waveformDD, waveformD ] = doubleDerive( waveform )
    waveformD = diff(waveform);
    waveformD = [waveformD NaN];

    waveformDD = diff(waveformD);
    waveformDD = [waveformDD NaN];
    waveformDD = BP_lowpass(waveformDD);

    %Perform the switch for positive and negative first-derivative values
    waveformDDPlus = waveformDD .* (waveformD > 0 & waveformDD > 0);
    waveformDDPlus = waveformDDPlus .^ 2;
end


function [fixedIndex] = FixIndex(BrokeIndex, BrokeSCG, Down, minWavelength)
    % follows a slope either up or down in a given window until an
    % extremum is attained
    fixedIndex = BrokeIndex;
    Radius = round(minWavelength/4);
    for N = 1:length(BrokeIndex)
       if BrokeIndex+round(minWavelength) < length(BrokeSCG)
           OldIndex = BrokeIndex(N);
           if Down
               TempSCG = BrokeSCG(max(OldIndex - Radius, 1):min(OldIndex + Radius, end));
               [~, NewIndex] = min(TempSCG);
               NewIndex = NewIndex + OldIndex - Radius - 1;
               while ~(OldIndex==NewIndex)
                   OldIndex = NewIndex;
                   TempSCG = BrokeSCG(max(OldIndex - Radius, 1):min(OldIndex + Radius,end));
                   [~, NewIndex] = min(TempSCG);
                   NewIndex = NewIndex + max(OldIndex - Radius, 1) -1;
               end
           else
               TempSCG = BrokeSCG(max(OldIndex - Radius, 1):min(OldIndex + Radius, end));
               [~, NewIndex] = max(TempSCG);
               NewIndex = NewIndex + OldIndex - Radius - 1;
               while ~(OldIndex==NewIndex)
                   OldIndex = NewIndex;
                   TempSCG = BrokeSCG(max(OldIndex - Radius, 1):min(OldIndex + Radius, end));
                   [~, NewIndex] = max(TempSCG);
                   NewIndex = NewIndex + max(OldIndex - Radius, 1) -1;
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

function [ dicroticIndex, notchIndex ] = getDicroticIndex( waveformDD, waveformD, bpwaveform, footIndex, systolicIndex )
    global Fs;
    RR = median(footIndex(2:end) - footIndex(1:end-1)) ./ Fs;% This assumes steady heartrate
    Down = 1; Up = ~Down;
    minWavelength = round(RR/5 .* Fs);
    notchIndex = FixIndex(systolicIndex + minWavelength, waveformDD, Up, minWavelength);
    dicroticIndex = FixIndex(notchIndex + round(0.5*minWavelength), waveformDD, Down, minWavelength);
    
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
            dicroticIndex(i) = FixIndex(notchIndex(i) + round(0.5*minWavelength), bpwaveform, Up, 4);
        end
    end
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

function [ res ] = winmean( rwin, quant )
    shape = size(rwin);
    vecsize = shape(2);
    res = zeros(1, vecsize);
    parfor i = 1:vecsize
        res(i) = quant * mean(rwin(:,i));
    end
end

function [ res ] = winsum( rwin )
    shape = size(rwin);
    vecsize = shape(2);
    res = zeros(1, vecsize);
    parfor i = 1:vecsize
        res(i) = sum(rwin(:,i));
    end
end