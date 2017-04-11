function [fixedIndex] = FixIndex(BrokeIndex, BrokeSCG, Down, minWavelength)
% function [fixedIndex] = FixIndex(BrokeIndex, BrokeSCG, Down, minWavelength)
    % follows a slope either up or down in a given window until an
    % extremum is attained
    fixedIndex = BrokeIndex;
    Radius = round(minWavelength/4);
    for N = 1:length(BrokeIndex)
       if BrokeIndex(N)+round(minWavelength) < length(BrokeSCG)
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