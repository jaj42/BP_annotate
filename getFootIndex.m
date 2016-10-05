function [ footIndex ] = getFootIndex( waveformDDPlus, zoneOfInterest )
%[ footIndex ] = getFootIndex( waveformDDPlus, zoneOfInterest )
%
%%
    zoneWall = diff( zoneOfInterest );
    BP_start = find(zoneWall == 1);
    BP_stop = find(zoneWall == -1);

    if BP_stop(1) < BP_start(1)
        BP_stop = BP_stop(2: end);
    end

    if BP_start(end) > BP_stop(end)
        BP_start = BP_start(1: end - 1);
    end

    footIndex = zeros( size(BP_start) );
    for i = 1 : length(BP_start)
        [~, footIndex(i)] = max(waveformDDPlus(BP_start(i) : BP_stop(i)));
        footIndex(i) = footIndex(i) + BP_start - 1;
    end


end

