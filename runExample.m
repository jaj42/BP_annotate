%Runs an example of BP_annotate
testPressures = struct;
Dir = {'/Volumes/Data/alexandrel/BP_annotate/'};
Files = dir(fullfile(Dir{1}, '*.mat'));
Files = {Files.name};
disp('Loading files...')
for i = 1 : length(Files)
    testPressures.(['P',num2str(i)]).waveform = load([Dir{1},Files{i}]);
    Field = fields(testPressures.(['P',num2str(i)]).waveform);
    testPressures.(['P',num2str(i)]).waveform = testPressures.(['P',num2str(i)]).waveform.(Field{1});
    testPressures.(['P',num2str(i)]).fs = str2double(Files{i}(3:regexp(Files{i},'Hz')-1));
end
disp('Done.')
verbose = 1;
%%
Timings = zeros(length(fields(testPressures)),2);
for i = 1: length(fields(testPressures))
    i
    Timings(i,1) = length(testPressures.(['P',num2str(i)]).waveform)/testPressures.(['P',num2str(i)]).fs/60;
    tic
    [ footIndex, systolicIndex, notchIndex, dicroticIndex ] = BP_annotate( testPressures.(['P',num2str(i)]).waveform, testPressures.(['P',num2str(i)]).fs, verbose );
    Timings(i,2) = toc;
%     pause
    close gcf
end
%%
figure(1); clf; plot(Timings(:,1),Timings(:,2),'o')
xlabel('length of record (min)')
ylabel('computing time (s)')
