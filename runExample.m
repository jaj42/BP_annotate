%Runs an example of BP_annotate
testPressures.P1.data = load('exampleWaveform_1kHz');
testPressures.P1.fs = 1000;


[ footIndex, systolicIndex, notchIndex ] = BP_annotate( testPressures.P1.data.waveform, testPressures.P1.fs, true );

%testPressures.P2.data = load('exampleWaveform_200Hz');
%testPressures.P2.data.waveform = testPressures.P2.data.waveform(1:1e5-1);
%testPressures.P2.fs = 200;
%
%verbose = 1;
%%%
%for i =2%: length(fields(testPressures))
%    [ footIndex, systolicIndex, notchIndex ] = BP_annotate( testPressures.(['P',num2str(i)]).data.waveform, testPressures.(['P',num2str(i)]).fs, verbose );
%end
