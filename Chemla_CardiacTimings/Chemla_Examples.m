close all

Dir = '/Volumes/Data/alexandrel/BP_annotate/';
File1 = 'BI_4.txt';
Data1 = dlmread([Dir,File1],'\t',1,0);
Data1 = Data1(:,1);

[ heartPeriod, LVET, diastolicTime, time, bpwaveform ] = Chemla_CardiacTimings( Data1, 1000,1);
File2 = 'DC_3.txt';
Data2 = dlmread([Dir,File2],'\t',1,0);
Data2 = Data2(:,1);
[ heartPeriod, LVET, diastolicTime, time, bpwaveform ] = Chemla_CardiacTimings( Data2, 1000,1);

File3 = 'SM_1.txt';
Data3 = dlmread([Dir,File3],'\t',1,0);
Data3 = Data3(:,1);
[ heartPeriod, LVET, diastolicTime, time, bpwaveform ] = Chemla_CardiacTimings( Data3, 1000,1);