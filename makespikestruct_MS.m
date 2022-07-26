function [spikestruct] = makespikestruct_MS(shared_drive,data_drive,db, exp, sp, method, need_stitching)
% saves spikestruct variable for quick use with subsequent analyse
% plots and saves specified LFP traces and visualises saturations into experiment directory

% Inputs:
% db - recording database (constructed by running makedb_TCB2.m)
% exp - experiment(s) to use
% sp - kilosort output from loadKSdir.m


% Output fields from makespikestruct.m function:
% spikestruct.date - date of recording
% spikestruct.meta = sub structure containing animal information such as: animal ID, age, weight, sex
% spikestruct.clusteridx - numbers of good units from kilosort output
% spikestruct.clusteridx1 - numbers of good units corresponding to phy
% spikestruct.raster - spiketrains of all units (1ms resolution, Neurons (rows) x Time (columns))
% spikestruct.condspikevector - spiketrains of all units split into conditions (cells)
% spikestruct.spikerate - spike rate (sp/s) split into conditions (columns) and neurons (rows)
% spikestruct.populationrate - sum of all raster rows
% spikestruct.populationrateall - population rate of all good units and MUA
% spikestruct.condpopspikevector - population rate split in to conditions (columns)
% spikestruct.LFP - LFP sampled (1Khz resolution) from channels specified in db(exp).lfp (rows)
% spikestruct.channel - the channel number on which the neuron was detected
% spikestruct.depth - the depth in uM of the channel on which the neuron was detected
% spikestruct.unitquality - the isolation distance (column 1) and contamination of refractory period (column 2)
% spikestruct.design - a character array specifying the nature of each condition
% spikestruct.timepointsms - time in milliseconds of the start of each condition and the finish time
% spikestruct.violationms - times in ms of LFP to exclude (channels in cells)
% spikestruct.waveforms - the maximum waveform of each unit (corresponding to clusteridx)
% spikestruct.spikeparams = column 1 - half amplitude duration, column 2 - trough to peak time, column 3 - trough vs peak height ratio
% spikestruct.StPR - StPR conditions in cells, units in rows (100ms lag)
% spikestruct.syringe_label = blind label on syringe
% spikestruct.syringe_contents = contents of syringe
% spikestruct.ref = referencing configuration
% spikestruct.gain = LFP gain setting used
% spikestruct.HTR = Head Twitch Resonse seen in 3x3 minute observations after recording
% spikestruct.dose - dose of TCB-2 given in mg/kg
% spikestruct.controldose - amount of control solution given in ml
% db.syringe_label = blind label on syringe
% db.syringe_contents = contents of syringe
% db.metadata = sub structure containing animal information such as: ID, age, weight, sex
% db.ref = referencing configuration
% db.LFPgain = LFP gain setting used
% db.HTR = Head Twitch Resonse seen in 3x3 minute observations after recording

% See R:\Neuropix\bd126\Code\spikestruct.txt for more info on fields
%% load experiment directory
%topDir = 'R:\Neuropix\bd126\Analysis\';
%expdir = [topDir 'Exp ' num2str(exp) ' ' db(exp).animal ' ' db(exp).injection{1:end-1}];
datafolder = [data_drive '\Data\' db(exp).animal '\' db(exp).date '\'];
topDir = [shared_drive '\cortical_dynamics\User\ms1121\Analysis Testing\'];
expdir = [topDir 'Exp_' num2str(exp) '_' db(exp).animal '_' db(exp).date];

%% set parameters
vStimSR = 9057.971014; % sampling rate for nidq
NIDQchan = 9;

%% Recording info from db
spikestruct.meta.date = db(exp).date;
spikestruct.meta.animal = db(exp).animal;
spikestruct.meta.age = db(exp).age;
spikestruct.meta.weight = db(exp).weight;
spikestruct.meta.sex = db(exp).sex;
spikestruct.dose = db(exp).dose;
spikestruct.controldose = db(exp).controldose;
spikestruct.syringe_label = db(exp).syringe_label;
spikestruct.syringe_contents = db(exp).syringe_contents;
spikestruct.ref = db(exp).ref;
spikestruct.LFPgain = db(exp).LFPgain;
spikestruct.HTR = db(exp).HTR;
spikestruct.design = db(exp).injection; % condition labels
spikestruct.timepoints = db(exp).timepointsms; % start time of each condition and end time

%% load LFP
if ~isfile([expdir '\LFP1kHz exp ' num2str(exp) ' Chan ' num2str(db(exp).lfp) '.mat']) % only load LFP if file doesn't exist already
    lfpchan = db(exp).lfp; % channels to save in variable
    SR = 2500;
    NP = 'NP1';
    if isfield(db(exp),'probe')
        if strcmp(db(exp).probe,'NP24')
            SR = 30000;
            NP = 'NP24';
        end
    end
    disp(['LFP SR = ' num2str(SR)]);
    lfpdir = db(exp).LFPdir;

    % find LFP file and add back median if subtracted
    disp('Finding lfp channels...');
    [lfpNOCAR1kHz, lfpNOCAR1kHzsaline, ~] = findlfp_MS(exp, lfpdir, lfpchan, db(exp).nChans{end}+db(exp).nChans{end-1}, SR, need_stitching, NP);

    disp('Saving lfp channels...');
    save([expdir '\LFP1kHz exp ' num2str(exp) ' Chan ' num2str(db(exp).lfp)], 'lfpNOCAR1kHz', 'lfpNOCAR1kHzsaline', '-v7.3');
else
    load([expdir '\LFP1kHz exp ' num2str(exp) ' Chan ' num2str(db(exp).lfp)])
end

spikestruct.LFP = lfpNOCAR1kHz;

if exist('lfpNOCAR1kHzsaline', 'var') & ~isempty(lfpNOCAR1kHzsaline)
    spikestruct.LFPsaline =  lfpNOCAR1kHzsaline;
end

if isfield(db(exp),'channelorder')
    spikestruct.channelorder = db(exp).channelorder;
end

%% find LFP saturations
for chan = 1:numel(spikestruct.LFP)
    disp(['Exp: ' num2str(exp) ', Channel: ' num2str(db(exp).lfp(chan))])
    
    % LFP threshold violations
    lfpthreshold{chan} = [];
    lfpthreshold{chan} = 10*std(spikestruct.LFP{chan});
    thresholdviolations{chan} = find(spikestruct.LFP{chan} > mean(spikestruct.LFP{chan})+lfpthreshold{chan}...
        | spikestruct.LFP{chan} < mean(spikestruct.LFP{chan}) - lfpthreshold{chan}); % times exceeding threshold (ms)
    
    % find LFP saturations using detection function
    LFPsaturations{chan} = [];
    [LFPsaturationslogical{chan}, time, ~, ~, ~, ~] =...
        detectLFPsaturationsMD_BD(spikestruct.LFP{chan}, 0.001, 'combined', 1, [0.2 0.25], ['Exp ' num2str(exp)], num2str(db(exp).lfp(chan)));
    LFPsaturations{chan} = single(time(logical(LFPsaturationslogical{chan}))*1000); % saturation times in ms
    close all
    
    % make saturation detection parameters less conservative for recordings with lower gain or internal referencing (look at total SD)
    
    % Add in padding for saturation points...
    
    % Threshold and saturation violations
    spikestruct.violationms{chan} = [];
    spikestruct.violationms{chan} = sort(unique([LFPsaturations{chan} thresholdviolations{chan}]));

end

%% get frame times for stimulus
if ~isempty(db(exp).NIDQ)
  disp('Finding frame times...')
  [frameTimes, ~] = detectViStim(db(exp).NIDQ{1}, NIDQchan, vStimSR); % gives frame times in seconds
  frameTimes = round(frameTimes*1000); % change to ms for 1kHz
  [frameTimes, frames_removed] = exclude_outlier_frames(frameTimes,6000); % 6000ms used as buffer for acceptable gap between frames
  if ~isempty(frames_removed)
    disp([num2str(numel(frames_removed)) ' frames needed to be removed. Frames removed were: ' num2str(frames_removed)])
  else
    disp('No frames needed to be removed.')
  end
  for cond = 1:numel(db(exp).injection)-1 % loop on conditions
    spikestruct.frameTimes{cond} = frameTimes(frameTimes > spikestruct.timepoints(cond) & frameTimes < spikestruct.timepoints(cond+1));
  end
end

%% get visual stim IDs and details
for c = 1:numel(db(exp).injection)
  if strcmp(db(exp).injection{c},'ViStimpre')
    load([datafolder db(exp).animal '_A.mat']); 
    spikestruct.ViStimID_natural = imagerandorder;
  end
  if strcmp(db(exp).injection{c},'ViStimpre_grating')
    stim_details = load([datafolder 'stim_details_' db(exp).animal '_' db(exp).date '_A.csv']);
    trial_IDs = load([datafolder 'trial_IDs_' db(exp).animal '_' db(exp).date '_A.csv']);
    spikestruct.ViStimID_grating{1} = trial_IDs;
    spikestruct.ViStimDetails_grating = stim_details;
  end
  if strcmp(db(exp).injection{c},'ViStimpost_grating')
    trial_IDs = load([datafolder 'trial_IDs_' db(exp).animal '_' db(exp).date '_B.csv']);
    spikestruct.ViStimID_grating{2} = trial_IDs;
  end
end

%% get rec map stimIDs and details
if ~isempty(db(exp).rec_map_files)
    [frameTimes, ~] = detectViStim(db(exp).rec_map_files{2}, NIDQchan, vStimSR); % gives frame times in seconds
    frameTimes = round(frameTimes*1000); % change to ms for 1kHz
    spikestruct.RecMap_frameTimes = frameTimes;

    stim_details = load([datafolder 'rec_map_stim_details_' db(exp).animal '_' db(exp).date '.csv']);
    trial_IDs = load([datafolder 'rec_map_stimID_' db(exp).animal '_' db(exp).date '.csv']);
    spikestruct.RecMapStimDetails = stim_details;
    spikestruct.RecMapStimID = trial_IDs;
    spikestruct.RecField = db(exp).rec_field;
end
  
%% find cluster IDs and make raster and population variables

if strcmp(method, 'auto') % in KS2.5 & KS3 sp.cids contains some units with 0 spikes so need to use clu
    clusteridx = sp.cids(sp.cgs == 2); % only good units
    clusteridxall = unique(sp.clu)'; % include MUA
    spikestruct.goodunits = ismember(clusteridxall, clusteridx)';
    % elseif strcmp(method, 'KS3')
    %       clusteridxall = sp.cids(sp.cgs < 3);
    %       spikestruct.goodunits = ismember(clusteridxall, clusteridx)';
elseif strcmp(method,'unsorted')
    disp('Reading unsorted KS unit classifications.');
    ksclusterfile = [db(exp).dir '\cluster_KSLabel.tsv'];
    [cids, cgs] = readClusterGroupsCSV(ksclusterfile);
    sp.cgs = cgs;
    clusteridxall = unique(sp.clu)';
    clusteridx = cids(sp.cgs == 2);
    spikestruct.goodunits = ismember(clusteridxall,clusteridx);
else
    clusteridx = sp.cids(sp.cgs == 2); % only good units
    clusteridxall = sp.cids(sp.cgs < 3);
    spikestruct.goodunits = (sp.cgs == 2)';
end

if isempty(clusteridx)
    spikestruct.clusteridx = []; % recording has no good units
    disp('Recording has no good units or is unsorted in KS 1 or 2')
    %     return
end

spikestruct.clusteridx = double(clusteridx)';
spikestruct.clusteridxall = double(clusteridxall)';

% raster for good units only
disp([num2str(numel(spikestruct.clusteridx)) ' good clusters in exp: ' num2str(exp)])
%  spikestruct.raster = zeros(numel(spikestruct.clusteridx), db(exp).timepointsms(end)+1); % Preset uint8

for clu = 1:numel(spikestruct.clusteridx) % loop on clusters
    spiketimes{clu} = sp.st(sp.clu == spikestruct.clusteridx(clu)); % time in s of every spike (units in rows)
    [spikestruct.raster(clu,:), edges] = histcounts(spiketimes{clu}*1000, 0:db(exp).timepointsms(end)+1); % number of spikes in each ms
    spikestruct.raster = uint8(spikestruct.raster);
    % disp(['Good cluster ' num2str(clu) '/' num2str(numel(spikestruct.clusteridx))])
end


comparefiles = 1:numel(db(exp).injection)-1; % number of conditions
disp([num2str(numel(spikestruct.clusteridxall)) ' total clusters in exp: ' num2str(exp)])
if ~isempty(clusteridx)
    spikestruct.raster(spikestruct.raster > 1) = 1; % Change all spikes to 1 if mislabelled as greater than one
    
    for  cond = comparefiles % loop on conditions
        spikestruct.condspikevector{cond} = spikestruct.raster(:,floor(db(exp).timepointsms(comparefiles(cond)))+1:...
            floor(db(exp).timepointsms(comparefiles(cond)+1))); % Split spike vectors into conditions (rows)
        for clu = 1:numel(spikestruct.clusteridx) % loop on clusters
            spikestruct.spikerate(clu,cond) = sum(spikestruct.condspikevector{cond}(clu,:))/db(exp).DurationS(comparefiles(cond));  % number of spikes per second
        end
    end
    
    % Population spike vector (number of spikes per ms)
    spikestruct.populationrate = uint8(sum(spikestruct.raster,1)); % Population rate made up of good units
end


% raster for all units
for clu = 1:numel(spikestruct.clusteridxall) % loop on clusters
    spikestruct.spiketimesall{clu,1} = sp.st(sp.clu == spikestruct.clusteridxall(clu)); % time in s of every spike (units in rows)
    [spikestruct.rasterall(clu,:), edges] = histcounts(spikestruct.spiketimesall{clu}*1000, 0:db(exp).timepointsms(end)+1); % number of spikes in each ms
    spikestruct.rasterall = uint8(spikestruct.rasterall);
    % disp(['MUA and good clusters: ' num2str(clu) '/' num2str(numel(spikestruct.clusteridxall))])
end
spikestruct.rasterall(spikestruct.rasterall > 1) = 1; % Change all spikes to 1 if mislabelled as greater than one

    for  cond = comparefiles % loop on conditions
        spikestruct.condspikevectorall{cond} = spikestruct.rasterall(:,floor(db(exp).timepointsms(comparefiles(cond)))+1:...
            floor(db(exp).timepointsms(comparefiles(cond)+1))); % Split spike vectors into conditions (rows)
        for clu = 1:numel(spikestruct.clusteridxall) % loop on clusters MUA
            spikestruct.spikerateall(clu,cond) = sum(spikestruct.condspikevectorall{cond}(clu,:))/db(exp).DurationS(comparefiles(cond));  % number of spikes per second
        end
    end
    
PopSpiketimems = sp.st*1000; % Population rate made up of all spikes (good units and MUA)
spikestruct.populationrateall = histcounts(PopSpiketimems, 0:db(exp).timepointsms(end)+1);

if isempty(spikestruct.clusteridx)    
    spikestruct.populationrate = histcounts(PopSpiketimems, 0:db(exp).timepointsms(end)+1);
    
end

for clu = 1:numel(comparefiles)
    spikestruct.condpopspikevector{1,clu} = uint8(spikestruct.populationrate(floor(db(exp).timepointsms(comparefiles(clu))+1):...
        floor(db(exp).timepointsms(comparefiles(clu)+1))));
end

%% Spike triggered population rate
xcorrms = 100;
tic
for cond = comparefiles % loop on conditions
    for clu = 1:numel(spikestruct.clusteridx)
        spikestruct.StPR{cond}(clu,:) = ... % StPR condition in cells, units in rows
            xcorr(spikestruct.condpopspikevector{cond}-spikestruct.condspikevector{cond}(clu,:),...
            spikestruct.condspikevector{cond}(clu,:), xcorrms)/sum(spikestruct.condspikevector{cond}(clu,:));
    end
    disp(['StPR condition: ' num2str(cond) '/' num2str(numel(comparefiles)) ' - ' db(exp).injection{cond}])
end
toc
%% Save struct without waveform and quality data
% Unit waveforms, channels and depth

% waveforms BD code needs testing (fast version to be used when post processing has not been run, channels may be out by a few)
[~, max_site] = max(max(abs(sp.temps),[],2),[],3); % the maximal site for each template (channels of largest magnitude point in waveform)
spike_ycoord = sp.ycoords(max_site(sp.spikeTemplates+1)); % find ycoordinates for each spike
channelind = max_site(sp.spikeTemplates+1); % find channel index  for each spike
[c, ia, ic] = unique(sp.clu); % find indices of unique cluster numbers
spike_ycoord = spike_ycoord(ia); % find ycoords of clusters
channelonprobe = spike_ycoord/10;
spiketemplates = sp.spikeTemplates+1;
spiketemplates = spiketemplates(ia); % find template of clusters
channelind = channelind(ia); % find channels of clusters

% (unfinished, sometimes needs sp.ycoords to add in cnannel depths that were removed in kilosort2 as bad channels)
if ~isempty(strfind(db(exp).dir, 'Neuropix')) && exp >= 62 % if neuropixels and kilosort 2
    probelayout = 20:10:3840;
    probelayout(2:2:end+1) = 20:20:3840;
    sp.ycoords = probelayout;
end

spikestruct.depthall = spike_ycoord;
spikestruct.channelall = channelonprobe;

if strcmp(method, 'manual')
    spikestruct.depth = spike_ycoord(sp.cgs==2); % good unit ycoords
    spikestruct.channel = channelonprobe(sp.cgs==2); % good unit channel
    qualityfilelist = dir([db(exp).dir '\*.qua.1.mat']);
    if isfile([db(exp).dir '\' qualityfilelist.name]) & isfile([db(exp).dir '\waveforms.mat'])
        load([db(exp).dir '\waveforms.mat'], 'maxWaveforms', 'chanMap', 'cluIDs');
        Channeldepth(:,1) = chanMap(:,2)-1; % channel on probe
        Channeldepth(:,2) = double(sp.ycoords(chanMap(:,2))); % depth
        noiseMUA = sum(spikestruct.clusteridx < 2); % number of units assigned indices 0 and 1 (these are used for noise and MUA in post processing)
        if noiseMUA > 0
            spikestruct.channel = circshift(Channeldepth(:,1), noiseMUA);
            spikestruct.depth = circshift(Channeldepth(:,2), noiseMUA);
            spikestruct.waveforms = circshift(maxWaveforms, noiseMUA, 1);
        else
            spikestruct.channel = Channeldepth(:,1);
            spikestruct.depth = Channeldepth(:,2);
            spikestruct.waveforms = maxWaveforms;
        end
        
        %     % waveform parameters
        %     for i = 1:size(maxWaveforms,1)
        %         spikestruct.spikeparams(i,:) = SpikeParams(maxWaveforms(i,:));
        %     end
        
        % Unit quality and contamination
        if isfile([db(exp).dir '\' qualityfilelist.name])  % load unit quality and contamination
            load([db(exp).dir '\' qualityfilelist.name], 'unitQ');
            stableunitquality = unitQ(ismember(unitQ(:,1), cluIDs),[1 2 6]); % Isolation distance (column 2) and quality (column 3) of units with waveforms
            spikestruct.unitquality = stableunitquality(:,2:3); % Isolation distance (column 1) and quality (column 2)
        end
    end
    
elseif strcmp(method, 'auto') || strcmp(method, 'unsorted')
    spikestruct.depth = spike_ycoord(spikestruct.goodunits); % good unit ycoords
    spikestruct.channel = channelonprobe(spikestruct.goodunits); % good unit channel
    %     spikestruct.depth = spike_ycoord(sp.cgs==2); % good unit ycoords
    %     spikestruct.channel = channelonprobe(sp.cgs==2); % good unit channel
end

% load waveforms based on the template and channel for each unit
% needs improving for good unit = 0
unitcount = 0;
for i = 1:numel(spikestruct.clusteridxall)
    unitcount = unitcount + 1;
    spikestruct.waveformsall(unitcount,:) = sp.temps(spiketemplates(i), :, channelind(i));
end

if isfield(spikestruct, 'waveforms')
    for i = 1:size(spikestruct.waveforms,1)
        spikestruct.spikeparams(i,:) = SpikeParamsBD(spikestruct.waveforms(i,:));
    end
end

%% sanity check waveforms.mat corresponds to sp
if isfield(spikestruct, 'depth') & ~(numel(spikestruct.clusteridxall) == numel(spikestruct.depthall))
    disp(['Error in exp ' num2str(exp) ' - parameters from sp and waveforms.mat do not match'])
end
end
%% test code
% makedb_TCB2
% sp = loadKSdir(db(48).dir);
% [spikestruct] = makespikestruct(db, 48, sp);

%% spare code
% if clusteridx(1) == 0; % Corrects for good unit labelled as 0
%     clusteridx(1) = max(clusteridx) + 1;
%     [clusteridx, idx]  = sort(clusteridx);
%     goodclusters = goodclusters(idx);
%     disp('Cluster number 0 in phy classified as good unit')
% end

% if clusteridxall(1) == 0;
%     clusteridxall(1) = max(clusteridxall) + 1;
%     [clusteridxall, idxall]  = sort(clusteridxall);
%     allclusters = allclusters(idxall);
% end

% find LFP saturations using BD detection function (for CAR data)
%[LFPsaturations{chan}, saturationTimes{chan}, nSaturations, fSaturations, meanSatDuration] =...
%detectLFPsaturationsBD(lfpCAR1kHz{chan}, 0.001, 'diff', ['Exp ' num2str(exp)], num2str(db(exp).lfp(chan)))
