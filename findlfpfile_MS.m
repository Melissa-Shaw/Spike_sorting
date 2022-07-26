function [LFP] = findlfpfile_MS(lfpdir)
%% function to look inside the LFP directory and 

lfpdirectory = dir([lfpdir, '\*.dat']);

if size(lfpdirectory,1) == 1 % separate LFP directory
% find LFP files and median trace files
    lfpdata = dir([lfpdir, '\*.dat']);
    lfpmedian = dir([lfpdir, '\*medianTrace.mat']);
elseif size(lfpdirectory,1) > 1 % combined LFP and spikes directory (2 files)
    lfpdata = dir([lfpdir, '\*LFP*.dat']);
    lfpmedian = dir([lfpdir, '\*LFP_medianTrace*.mat']);
else
    disp(['No LFP files found in ' lfpdir])
end

LFP = [lfpdata.folder '\' lfpdata.name];
end