function [ footIndex, systolicIndex, notchIndex, dicroticIndex, time, bpwaveform ] = BP_annotate( inWaveform, inFs, verbose, varargin )
%% [ footIndex, systolicIndex, notchIndex, dicroticIndex, time, bpwaveform ] = BP_annotate( inWaveform, inFs, verbose, Units, isClean )
% Implementation of a feature detection algorithm for arterial blood
% pressure in humans. The foot of the wave, systolic peak, dicrotic notch,
% and dicrotic peaks are identified. The blood pressure time series is
% always resampled at 200 Hz to allow standardisation.
%
% The technique was largely inspired by the derivatives and thresholds 
% described in Pan-Tompkins:
% Pan, Jiapu, and Willis J. Tompkins. "A real-time QRS detection algorithm." 
% IEEE transactions on biomedical engineering 3 (1985): 230-236.
%
% We also use criteria described by Sun, J. X., A. T. Reisner, and R. G. Mark. in
% "A signal abnormality index for arterial blood pressure waveforms."
% Computers in Cardiology, 2006. IEEE, 2006.
% These criteria have not been extensively tested by us.
%% Inputs
% inWaveform : countinuous arterial blood pressure time-series
% inFs : sampling frequency (Hz) of the time-series
% verbose : boolean, should be true if figures are wanted
% 
%% Otional Inputs
% Units : string, 'mmHg' or other, only mmHg are dealt with so far
% isClean : boolean, should be true if confidence in the entire signal is
% high. It allows annotation before the threshold window has initiated.
%
%% Outputs
% footIndex : index of the foot of each systolic wave
% systolicIndex : index of each sytolic peak
% notchIndex : index of each dicrotic notch
% dicroticIndex : index of each dicrotic peak
% time : time vector (s) of the resampled 200 Hz time-series
% bpwaveform : resampled filtered 200 Hz time-series
%
%% Methods
%% Foot index
% The foot index is defined as the point where the second derivative of 
% the time-series is the highest in each interval where a moving average
% of the second derivative was bigger than a adaptative threshold. This 
% criterion was preferred over others because of its robustness. It is
% possible to move this index to the minimum of the signal in this region
% by uncommenting line 302 in the source code.
%% Systolic peak index
% The systolic peak is defined as the maximum of the waveform following 
% the foot index, relative to a window of radius 1/8 s around itself. 
%% Dicrotic notch and peak indices
% The dicrotic notch is defined as the minimum of the subtraction of the signal and the staight
% line going from systole to diastole. Dicrotic peak indices were defined as the minimum of the
% second derivative of the time-series following the dicrotic notch, relative to a window of 
% radius RR/5 s around itself. (RR is the median heartbeat interval computed form the 
% foot indices). These indices are moved to waveform minima and maxima if
% these exist.
%
%% Authors
% Alexandre Laurin, PhD, ?cole Polytechnique, France, alexandre.laurin@inria.fr
% Jona Joachim, MD, H?pital Lariboisi?re, France, jona.joachim@inria.fr
%%
% varargin handling
    Units = 'unkown';
    isClean = 0;
    if nargin >= 4
        Units = varargin{1};
    end
    if nargin >=5
        isClean = varargin{2};
    end
    
% setting global parameters    
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
    
    % isClean forces the identification of indices before the window has
    % had time to initiate
    if isClean
        firstNotNan = find(~isnan(threshold));
        firstNotNan = mean(threshold(firstNotNan(1) : firstNotNan(1) + integwinsize ));
        threshold(isnan(threshold)) = firstNotNan.*ones(size(threshold(isnan(threshold))));
    end
    
    % each zone of interest corresponds to a heart beat
    zoneOfInterest =  BP_integral > threshold ;
    
    % implementation of Sun et al. paper A Signal Abnormality Index for
    % Arterial Blood Pressure Waveform
    if strcmp(Units, 'mmHg')
        zoneOfInterest(bpwaveform > 300) = 0;
        zoneOfInterest(bpwaveform < 20) = 0;
        zoneOfInterest(noiseLevel < -40) = 0;
    end
    
    footIndex = getFootIndex( waveformDDPlus, zoneOfInterest );

    Down = 1; Up = ~Down;
    systolicIndex = FixIndex(footIndex + floor(integwinsize/2), bpwaveform, Up, floor(integwinsize/2));

    % implementation of Sun et al. paper A Signal Abnormality Index for
    % Arterial Blood Pressure Waveform
    if strcmp(Units, 'mmHg')
        meanPressure = (bpwaveform(systolicIndex) + bpwaveform(footIndex))./2;
        pulsePressure = (bpwaveform(systolicIndex) - bpwaveform(footIndex));
        footIndex( (meanPressure < 30) | (meanPressure > 200) | (pulsePressure < 20) ) = [];
        systolicIndex( (meanPressure < 30) | (meanPressure > 200) | (pulsePressure < 20) ) = [];
    end
    
    [ dicroticIndex, notchIndex ] = getDicroticIndex( waveformDD, waveformD, bpwaveform, footIndex, systolicIndex );  
    footIndex = footIndex(1:length(notchIndex));
    systolicIndex = systolicIndex(1:length(notchIndex));
    
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
        ylabel('arterial pressure')

        axs(2) = subplot(2, 1, 2);
        hold on;
        plot(time, waveformDD);
        plot(time, BP_integral);
        plot(time, threshold);
        plot(time, zoneOfInterest .* .1);
        legend({'2nd Derivative','Integral', 'Threshold', 'ZOI'}, 'box','off');
        xlabel('time (s)')
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


function [fixedIndex] = FixIndex(BrokeIndex, Signal, Down, minWavelength)
    % follows a slope either up or down in a given window until an
    % extremum is attained
    fixedIndex = BrokeIndex;
    Radius = round(minWavelength/4);
    for N = 1:length(BrokeIndex)
       if BrokeIndex(N)+round(minWavelength) < length(Signal)
           OldIndex = BrokeIndex(N);
           if Down
               TempSignal = Signal(max(OldIndex - Radius, 1):min(OldIndex + Radius, end));
               [~, NewIndex] = min(TempSignal);
               NewIndex = NewIndex + OldIndex - Radius - 1;
               while ~(OldIndex==NewIndex)
                   OldIndex = NewIndex;
                   TempSignal = Signal(max(OldIndex - Radius, 1):min(OldIndex + Radius,end));
                   [~, NewIndex] = min(TempSignal);
                   NewIndex = NewIndex + max(OldIndex - Radius, 1) -1;
               end
           else
               TempSignal = Signal(max(OldIndex - Radius, 1):min(OldIndex + Radius, end));
               [~, NewIndex] = max(TempSignal);
               NewIndex = NewIndex + OldIndex - Radius - 1;
               while ~(OldIndex==NewIndex)
                   OldIndex = NewIndex;
                   TempSignal = Signal(max(OldIndex - Radius, 1):min(OldIndex + Radius, end));
                   [~, NewIndex] = max(TempSignal);
                   NewIndex = NewIndex + max(OldIndex - Radius, 1) -1;
               end
           end
           Index = NewIndex;
           fixedIndex(N) = Index;
       end
    end
    fixedIndex( fixedIndex > length(Signal) ) = [];
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
    
    Down = 1; Up = ~Down;
    % Depending on the morphology of the signal, the foot index may be
    % identified as the minimum of the waveform itself
    %{ uncomment for local signal minimum identification
    % footIndex = FixIndex(footIndex, bpwaveform, Down, 10);
    %}
end

function [ dicroticIndex, notchIndex ] = getDicroticIndex( waveformDD, waveformD, bpwaveform, footIndex, systolicIndex )
    global Fs;
    RR = median(footIndex(2:end) - footIndex(1:end-1)) ./ Fs;% This assumes steady heartrate
    Down = 1; Up = ~Down;
    minWavelength = round(RR/5 .* Fs);
    
    straightLines = zeros(size( bpwaveform) );
    notQuiteSystolicIndex = FixIndex(systolicIndex, bpwaveform, Up, minWavelength);
    
    for i = 1 : length(footIndex) - 1
        % compute the parameters of the straight line going from systole to
        % diastole
        slope = ( bpwaveform(footIndex(i+1)) - bpwaveform(notQuiteSystolicIndex(i)) ) / ( footIndex(i+1) - notQuiteSystolicIndex(i) );
        intercept = bpwaveform(footIndex(i+1)) - slope * footIndex(i+1);
        straightLines( footIndex(i) : notQuiteSystolicIndex(i) ) = bpwaveform( footIndex(i) : notQuiteSystolicIndex(i) );
        straightLines( notQuiteSystolicIndex(i) : footIndex(i+1) ) = slope * ( notQuiteSystolicIndex(i) : footIndex(i+1) ) + intercept;
    end
    eyeBallSignal = bpwaveform - straightLines;
    eyeBallSignal(footIndex(end) : end) = waveformDD(footIndex(end) : end);
    
    notchIndex = FixIndex(systolicIndex + round(minWavelength), eyeBallSignal, Down, minWavelength/4);
    dicroticIndex = FixIndex(notchIndex + round(0.25*minWavelength), waveformDD, Down, round(0.25*minWavelength));
    systolicIndex = systolicIndex(1:length(dicroticIndex));
    
    % if a local minimum and maximum exist, move the dicrotic indices to
    % these
    for i = 1 : length(systolicIndex)
        Start = systolicIndex(i) + round(minWavelength/2);
        End = min([dicroticIndex(i) + round(minWavelength/4), length(waveformD)]);
        ZOI = waveformD(Start : End);
        ZOI = ZOI(2:end).*ZOI(1:end-1);
        extrema = find(ZOI<0);
        if length(extrema) >=2
            notchIndex(i) = FixIndex(notchIndex(i), bpwaveform, Down, 4);
            dicroticIndex(i) = FixIndex(min(notchIndex(i) + round(0.25*minWavelength), length(bpwaveform)), bpwaveform, Up, 4);
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
    for i = 1:vecsize
        res(i) = quant * mean(rwin(:,i));
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