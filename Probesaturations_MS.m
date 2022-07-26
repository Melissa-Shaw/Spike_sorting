function [LFPallchannels, LFPsaturations] = Probesaturations_MS(db, exp, nChans, chans2ignore, SR)
% Probesaturations runs saturation detection function on whole probe to find best LFP channels to use for analyses
% inputs:
% LFPfile - name of LFP file
% nChans - number of channels in file
% chans2ignore - vhannels corresponding to analogue inputs - leave empty if none
% SR - sampling rate in Hz or file

% Outputs:
% structure LFPsaturations including fields for:
% time of saturations (milliseconds)
% number of saturations
% frequency of saturations (per minute)
% mean duration of saturations (seconds)
% total time of saturations (seconds)

d = dir([db(exp).LFPdir '\*.dat*']);

% find length of file
LFPfile = [d.folder '\' d.name];
nSampsTotal = d.bytes/nChans/2;

% set logcal for channels to include
chans = ones(size(1:nChans));
chans(chans2ignore) = 0;
chans = logical(chans);

fid = fopen(LFPfile, 'r');
tic
LFPallchannelsraw = fread(fid, [nChans nSampsTotal], '*int16');
LFPallchannelsraw = LFPallchannelsraw(chans,:);
disp(['All channels have been loaded'])
toc


tic
for chan = find(chans == 1) % channels of interest
    LFPallchannels(chan,:) = resample(double(LFPallchannelsraw(chan,:))', 1000, 2500)'; % from 2500Hz to 1KHz resolution
end
disp(['All channels have been resampled to 1kHz'])
% tic 
% LFPallchannels_khz = resample(double(LFPallchannels)', 100, 1000)';
toc


for chan = find(chans == 1)
    [LFPsaturationslogical{chan,1}, time{chan,1}, nSaturations{chan,1}, fSaturations{chan,1}, dSaturations{chan,1}, f] =...
        detectLFPsaturationsMD_BD(LFPallchannels(chan,:), 0.001, 'combined', 1, [0.2 1], ['Exp ' num2str(exp)], num2str(chan));
    tSaturations{chan} = single(time{chan}(logical(LFPsaturationslogical{chan}))*1000); % saturation times in ms
    close all
end
disp(['Saturations have been detected for all channels'])

for zerosats = find(cellfun(@isempty, tSaturations)) % change empty cells to 0
    tSaturations{zerosats} = 0;
    dSaturations{zerosats} = 0;
end

LFPsaturations.tSaturations = tSaturations';
LFPsaturations.nSaturations = cell2mat(nSaturations);
LFPsaturations.fSaturations = cell2mat(fSaturations);
LFPsaturations.dSaturations = cell2mat(dSaturations);
LFPsaturations.totaltimeSaturations = cell2mat(nSaturations).*cell2mat(dSaturations);

end
