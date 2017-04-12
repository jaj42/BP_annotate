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

