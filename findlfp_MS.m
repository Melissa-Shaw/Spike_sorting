function [lfpNOCAR1kHz, lfpNOCAR1kHzsaline, LFP] = findlfp_MS(exp, lfpdir, chans2load, totalchans, SR, need_stitching, NP)
% function to look inside the LFP directory and add back the median trace to the specified channels if it is present
% otherwise just outputs LFP trace of channels of interest

% Inputs:
% exp - experiment in db
% lfpdir - directory of lfp file
% chans to load  - channels to load
% totalchans - total channels in lfp file including

% Outputs:
% lfpNOCAR1kHz - cell array of LFP resampled t 1kHz (channels in cells)
% lfpNOCAR1kHzsaline - LFP trace outside brain (estimeted)
% LFP - cell array of file info structs with LFP file in position 1 and median trace in position 2 (if exists)


if exp > 15 & exp < 30 % concatenated dual probe experiments on martynas' animals with channels that needed adding back
    lfpdirectory = dir([lfpdir, '\*swappedCAR.dat']);
    lfpdata = dir([lfpdir, '\*swappedCAR.dat']);
    lfpmedian = dir([lfpdir, '\*medianTrace.mat']);
else
    lfpdirectory = dir([lfpdir, '\*.dat']);
end

if size(lfpdirectory,1) == 1 % separate LFP directory
    % find LFP files and median trace files
    lfpdata = dir([lfpdir, '\*.dat']);
    lfpmedian = dir([lfpdir, '\*medianTrace.mat']);
elseif size(lfpdirectory,1) > 1 % combined LFP and spikes directory (more than one .dat file present)
    lfpdata = dir([lfpdir, '\*LFP*.dat']);
    lfpmedian = dir([lfpdir, '\*LFP_medianTrace*.mat']);
end

if need_stitching == false % anaesthetised neuropixels use raw .bin files as they dont need stiching
    lfpdirectory = dir([lfpdir, '\*.bin']);
    lfpdata = dir([lfpdir, '\*.bin']);
    lfpmedian = dir([lfpdir, '\*medianTrace.mat']);
end

if strcmp(NP,'NP24') % NP24 recordings have one file for both spikes and lfp
    lfpdirectory = dir([lfpdir, '\*.ap.bin']);
    lfpdata = dir([lfpdir, '\*.ap.bin']);
    lfpmedian = dir([lfpdir, '\*medianTrace.mat']);
end

if size(lfpdirectory,1) == 0
disp(['No LFP files found in ' lfpdirectory.folder])
end

LFP = {[lfpdata.folder '\' lfpdata.name], [lfpmedian.folder '\' lfpmedian.name]};

% is median trace present in LFP folder?
if isfile(LFP{2})
    load(LFP{2})  % median trace
    count = 0;
    for lfpchan = chans2load
        count = count + 1;
        lfpCAR{count} = getContinuousDataFromDatfile(LFP{1}, totalchans, 0, +inf, lfpchan, SR); % load median subtracted LFP
        lfpNOCAR{count} = lfpCAR{count} + medianTrace; % Add back subtracted median
        disp(['Adding median trace for file ' LFP{1} ' chan: ' num2str(lfpchan)])
        lfpNOCAR1kHz{count}  = resample(lfpNOCAR{count}(1:end-1), 1000, SR);% from 2500Hz to 1KHz resolution (end-1 to correspond to raster)
    end
    lfpNOCAR1kHzsaline = [];
    
else % if no median trace file found
    count = 0;
    for lfpchan = chans2load % loop on lfp channels to analyse
        count = count + 1;
        lfpNOCAR{count}  = getContinuousDataFromDatfile(LFP{1}, totalchans, 0, +inf, lfpchan, SR); % load raw LFP
        lfpNOCAR1kHz{count}  = resample(lfpNOCAR{count}(1:end-1), 1000, SR);% from 2500Hz to 1KHz resolution
        disp(['Loading LFP trace for file ' LFP{1} ' chan: ' num2str(lfpchan)])
    end
    lfpNOCARsaline  = getContinuousDataFromDatfile(LFP{1}, totalchans, 0, +inf, 330, SR); % load raw LFP
    lfpNOCAR1kHzsaline  = resample(lfpNOCARsaline(1:end-1), 1000, SR);% from 2500Hz to 1KHz resolution
    disp(['Adding lfp in saline trace for file: ' LFP{1}])
end
end