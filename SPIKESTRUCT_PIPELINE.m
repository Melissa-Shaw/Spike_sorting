function  [spikestruct] = SPIKESTRUCT_PIPELINE(shared_drive,data_drive,db, exp, method, need_stitching) 
% method = 'manual' or 'auto' or 'unsorted'
% need_stitching = false or true (single recording or in patches) 
% shared_drive = 'S:' or 'X:'
% data_drive = 'F:' or 'E:'

addpath(genpath([shared_drive '\cortical_dynamics\Shared\Code\sortingPipeline']));
addpath(genpath([shared_drive '\cortical_dynamics\Shared\Code\github_cortex-lab_spikes']));
addpath([shared_drive '\cortical_dynamics\Shared\Code\github_kwikteam_npy-matlab']);
addpath([shared_drive '\cortical_dynamics\Shared\Code\matlib\IO']);

%% Make spikestruct
disp(['Making spikestruct for exp: ' num2str(exp)]);
topDir = [shared_drive '\cortical_dynamics\User\ms1121\Analysis Testing\'];
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
    case 'unsorted'
        KSdir = db(exp).dir;
    otherwise
        KSdir = input('Kilosort directory not found please specify manually the kilosort directory from which to construct spikestruct.mat:');
end

% create spikestruct variable
sp = loadKSdir(KSdir);
[spikestruct] = makespikestruct_MS(shared_drive,data_drive,db, exp, sp, method, need_stitching);

% save name according to KS directory
[~, spikestruct.kilosort, ~] = fileparts(KSdir);

if sum(ismember(KSdir(end-3: end), 'KS25')) == 4
    save([expdir '\spikestruct - KS25'], 'spikestruct', '-v7.3');
else
    save([expdir '\spikestruct'], 'spikestruct', '-v7.3');
    disp('Spikestruct saved.');
end
disp(['Exp: ' num2str(exp) ' spikestruct created!'])

%% Find frame times if pupil camera signal is present

%if isfield(spikestruct, 'frameTimes')
%    spikestruct = rmfield(spikestruct, 'frameTimes');
%end

%if db(exp).nChans{end} > 0 % if camera signal present
%    if exp < 48
%        SR = 3e4;
%        lfpdir = db(exp).dir;
%    else
%        SR = 2500;
%        lfpdir = db(exp).LFPdir;
%    end%

%    [LFP] = findlfpfile_MS(lfpdir);
%    [frameTimes, ~] = detectFrames_MS(LFP, db(exp).nChans{1} + db(exp).nChans{end}, SR);
%    hgsave(gcf, [expdir '\pupil_frames'])

%    for cond = 1:numel(db(exp).injection)-1 % loop on conditions
%        spikestruct.frameTimes{cond} = frameTimes(frameTimes > spikestruct.timepoints(cond)/1000 & frameTimes < spikestruct.timepoints(cond+1)/1000)*1000;
%    end
%    save([expdir '\spikestruct'], 'spikestruct', '-v7.3');
%    disp(['Recording: ' expdir ' ' num2str(numel(frameTimes)) ' frames detected'])
%else
%    disp(['Recording: ' expdir ' has no camera signal'])
%end
end