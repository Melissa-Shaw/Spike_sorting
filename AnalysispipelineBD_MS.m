function  AnalysispipelineBD_MS(db, experiments, analysis, calc_psd)
% AnalysispipelineBD
% To be used after spikesorting (can be used to add LFP to spikestruct before spikesorting by running first loop)
% Makes spikestruct.mat experiment datafile
% Calculates PSD for all channels
% Runs post processing pipeline and adds output to spikestruct
% Finds LFP saturations in all channels and adds to spikestruct
% runs analysis scripts (optional)

addpath(genpath('S:\cortical_dynamics\Shared\Code\sortingPipeline'));
addpath(genpath('S:\cortical_dynamics\Shared\Code\github_cortex-lab_spikes'));
addpath('S:\cortical_dynamics\Shared\Code\github_kwikteam_npy-matlab');
addpath('S:\cortical_dynamics\Shared\Code\matlib\IO');

%% Make spikestruct
method = 'manual';
region  = 'PFC';

for exp = experiments
    % topDir = 'R:\Neuropix\bd126\Analysis\';
    topDir = 'E:\ms1121\Analysis Testing\';
    %expdir = [topDir 'Exp ' num2str(exp) ' ' db(exp).animal ' ' db(exp).injection{1:end-1}];
    expdir = [topDir 'Exp_' num2str(exp) '_' db(exp).animal '_' db(exp).date];
    mkdir(expdir);
    addpath(genpath(expdir))
    
    % load kilosort output struct SP
    switch method
        case 'manual'
            KSdir = db(exp).dir;
        case 'auto'
            if exp > 47  & exp < 62
                KSdir = [db(exp).dir(1:end-8) 'PFC\KS25']; % Automated sorting KS25 dir Batch1 PFC
            elseif exp < 100
                KSdir = [db(exp).dir(1:end-8) 'KS25']; % Automated sorting KS25 dir Batch2 PFC
            end
            
            % Adaptation
            if exp > 99 & exp < 110
                KSdir = [db(exp).dir(1:end-32) 'PFC\KS25']; % Adaptation1 KS25 sorting dir
            elseif exp > 109 & exp < 130
                KSdir = [db(exp).dir(1:end-8) 'KS25']; % Adaptation2/3 KS25 sorting dir
            end
            
            % Slice experiments or anaesthetised
            if exp > 129
                KSdir = db(exp).dir;
            end
                
        otherwise
            KSdir = input('Kilosort directory not found please specify manually the kilosort directory from which to construct spikestruct.mat:');
    end
    
    % create spikestruct variable
    sp = loadKSdir(KSdir);
    [spikestruct] = makespikestruct_MS(db, exp, sp, method);
    
    % save name according to KS directory
    [~, spikestruct.kilosort, ~] = fileparts(KSdir);
     
    if sum(ismember(KSdir(end-3: end), 'KS25')) == 4
        save([expdir '\spikestruct - KS25'], 'spikestruct', '-v7.3');
    else
        save([expdir '\spikestruct'], 'spikestruct', '-v7.3');
    end
    disp(['Exp: ' num2str(exp) ' spikestruct created!'])
    clear spikestruct
end

%% Find frame times if pupil camera signal is present
for exp = experiments
    topDir = 'E:\ms1121\Analysis Testing\';
    expdir = [topDir 'Exp_' num2str(exp) '_' db(exp).animal '_' db(exp).date];
    load([expdir '\spikestruct']);
    
    if isfield(spikestruct, 'frameTimes')
        spikestruct = rmfield(spikestruct, 'frameTimes');
    end
    
    if db(exp).nChans{end} > 0 % if camera signal present
        if exp < 48
            SR = 3e4;
            lfpdir = db(exp).dir;
        else
            SR = 2500;
            lfpdir = db(exp).LFPdir;
        end
        
        [LFP] = findlfpfile_MS(lfpdir);
        [frameTimes, ~] = detectFrames_MS(LFP, db(exp).nChans{1} + db(exp).nChans{end}, SR);
        hgsave(gcf, [expdir '\pupil_frames'])
        
        for cond = 1:numel(db(exp).injection)-1 % loop on conditions
            spikestruct.frameTimes{cond} = frameTimes(frameTimes > spikestruct.timepoints(cond)/1000 & frameTimes < spikestruct.timepoints(cond+1)/1000)*1000;
        end
        save([expdir '\spikestruct'], 'spikestruct', '-v7.3');
        disp(['Recording: ' expdir ' ' num2str(numel(frameTimes)) ' frames detected'])
    else
        disp(['Recording: ' expdir ' has no camera signal'])
    end
end
%% calculate PSD for all channels and frequencies
if calc_psd
  for exp = experiments
      LFPfile = dir([db(exp).LFPdir '\*.dat*']);
      if ~exist([LFPfile.folder '\chanPSDs.mat'], 'file')
          psdChannels_MS([LFPfile.folder '\' LFPfile.name], 385, 385, 2500);
          disp(['Exp: ' num2str(exp) ' PSD all channels complete!'])
      else
          disp(['Exp ' num2str(exp) ' - ' LFPfile.folder '\chanPSDs.mat already exists'])
      end
  end
end
%% find saturations on every channel
for exp = experiments
    expdir = [topDir '\Exp_' num2str(exp) '_' db(exp).animal '_' db(exp).date];
    load([expdir '\spikestruct']);
    if ~isfield(spikestruct, 'saturations')
        if isfile([expdir '\LFPallchannels'])
            load([expdir '\LFPallchannels']);
        else
            [LFPallchannels, LFPsaturations] = Probesaturations(db, exp, 385, 385, 2500);
            save([expdir '\LFPallchannels'], 'LFPallchannels', '-v7.3');
        end
        spikestruct.saturations = LFPsaturations;
        save([expdir '\LFPallchannels'], 'LFPallchannels', '-v7.3');
        save([expdir '\spikestruct'], 'spikestruct', '-v7.3');
    end
    disp(['Exp: ' num2str(exp) ' All channels saturations complete and saved to spikestruct!']);
end
%% Extracts waveform and calculates quality for each unit and saves to existing spikestruct variable
for exp = experiments
    expdir = [topDir '\Exp ' num2str(exp) ' ' db(exp).animal ' ' db(exp).injection{1:end-1}];
    qualityfilelist = dir([db(exp).dir '\*.qua.1.mat']);
    if ~isfile([db(exp).dir '\' qualityfilelist.name]) & ~isfile([db(exp).dir '\waveforms.mat']) % only run if files not present already
        % run post processing pipeline
        postprocessingPipeline([db(exp).dir '\stitchedData.dat'], db(exp).nChans{2}, 3e4);
    else
        disp(['Either or both ' db(exp).dir '\waveforms.mat or ' db(exp).dir '\' qualityfilelist.name...
            ' already exists - check if spikestruct already includes waveforms and quality fields']);
    end
    
    % sanity check - post processing run on most recent spikesorting
    load([expdir '\spikestruct'])
    if isfield(spikestruct, 'depth') & ~(numel(spikestruct.clusteridx) == numel(spikestruct.depth))
        disp(['Error in exp ' num2str(exp) ' - parameters from sp and waveforms.mat do not match, re-running code to update'])
    end
    
    qualityfilelist = dir([db(exp).dir '\*.qua.1.mat']);% find quality file
    
    % Unit waveforms, channels and depth
    load([db(exp).dir '\waveforms.mat'], 'maxWaveforms', 'chanMap');
    sp = loadKSdir(db(exp).dir);
    Channeldepth(:,1) = chanMap(:,2)-1; % channel on probe
    if ~isempty(strfind(db(exp).dir, 'Neuropix')) && exp >= 62 % if neuropixels and kilosort 2
        sp.ycoords = 20:10:3840;
        sp.ycoords(2:2:end+1) = 20:20:3840;
    end
    Channeldepth(:,2) = double(sp.ycoords(chanMap(:,2))); % depth
    
    noiseMUA = sum(spikestruct.clusteridx < 2);
    if noiseMUA > 0
        spikestruct.channel = circshift(Channeldepth(:,1), noiseMUA);
        spikestruct.depth = circshift(Channeldepth(:,2), noiseMUA);
        spikestruct.waveforms = circshift(maxWaveforms, noiseMUA, 1);
    else
        spikestruct.channel = Channeldepth(:,1);
        spikestruct.depth = Channeldepth(:,2);
        spikestruct.waveforms = maxWaveforms;
    end
    
    % waveform parameters
    for i = 1:size(maxWaveforms,1)
        spikeparameters(i,:) = SpikeParamsBD(maxWaveforms(i,:));
    end
    spikestruct.spikeparams = spikeparameters;
    
    % Unit quality and contamination
    load([db(exp).dir '\' qualityfilelist.name], 'unitQ');
    stableunitquality = unitQ(ismember(unitQ(:,1), spikestruct.clusteridx),[1 2 6]); % Isolation distance (column 2) and quality (column 3) of units with waveforms
    spikestruct.unitquality = stableunitquality(:,2:3); % Isolation distance (column 1) and quality (column 2)
    if ~(numel(spikestruct.clusteridx) == numel(spikestruct.depth))
        disp(['Error in exp ' num2str(exp) ' - parameters from sp and waveforms.mat still different - delete files, check phy and re run post processing'])
    else
        disp([db(exp).dir ' parameters from sp and waveforms.mat match - waveforms and quality fields added to spikestruct'])
        save([expdir '\spikestruct'], 'spikestruct', '-v7.3');
        clear maxWaveforms chanMap Channeldepth spikeparameters
    end
end
%% Analysis scripts to be used with output
if analysis
    for exp = experiments
        if region == 'PFC'
            spikeanalysisPFC(db, exp);
        elseif region == 'V1'
            spikeanalysisV1;
        else
            disp(['Region: ' region ' not recognised - should be PFC or V1'])
        end
        poweranalysis(db, exp);
        % singleunitanalysis(db, exp, db(exp).cond, 'PFC')
        % StLFP(db, experiments, conditions, LFPchannel, frequencyband, plotfig, varargin)
        % frequencybandStLFP(db, exp, conditions, 0, [60 90])
    end
end
end