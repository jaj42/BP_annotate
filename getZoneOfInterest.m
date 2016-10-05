function [ zoneOfInterest ] = getZoneOfInterest( integralWindow, thresholdWindow )
%[ squareWave ] = getSquareWave( integralWindow, quantileWindow )
%
%Returns 1 if the windowed integral is bigger than the windowed threshold;
%returns 0 otherwise. This square wave describes the zones of interest
%%
zoneOfInterest = zeros( size(thresholdWindow) );
zoneOfInterest( integralWindow > thresholdWindow ) = 1;

end

