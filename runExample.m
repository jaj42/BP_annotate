%Runs 2 examples of BP_annotate
signal1 = load('exampleWaveform_1kHz.mat');
[ footIndex1, systolicIndex1, notchIndex1, dicroticIndex1 ] = ...
    BP_annotate( signal1.waveform, 1000, 1, 'mmHg', 1);

signal2 = load('exampleWaveform_200Hz.mat');
[ footIndex2, systolicIndex2, notchIndex2, dicroticIndex2 ] = ...
    BP_annotate( signal2.waveform, 200, 1, 'other', 1);