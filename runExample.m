%Runs an example of BP_annotate
load('exampleWaveform_1kHz');
fs = 1000;
verbose = 1;
%%
[ footIndex, systolicIndex, notchIndex ] = BP_annotate( waveform, fs, verbose );