function [FixedIndex] = FixIndex(BrokeIndex, BrokeSCG, Down, minWavelength)
%function [FixedIndex] = FixIndex(BrokeIndex, BrokeSCG, Down, minWavelength)
%%
    FixedIndex = BrokeIndex;
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
           FixedIndex(N) = Index;
       end
    end
    FixedIndex( FixedIndex > length(BrokeSCG) ) = [];
    FixedIndex( FixedIndex < 1 ) = [];
end