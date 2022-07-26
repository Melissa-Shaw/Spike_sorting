function [frameTimes, frameTimesoffset, frameInd] = detectViStim(dataFilename, nChans, sr)
% A helper function for getvisualresponse.m script. It finds the times and indices
% of every consecutive peak in the voltage data recorded  on NIDQ 
% These peaks correspond to photodiode signals corresponding to image
% presentations
% Inputs: 
% dataFilename - a string with the name of the data file (in dat or bin formats).
% nChans - a total number of channels in the recording. 
% sr - a sampling rate.
% Outputs: frameTimes and frameInd.

chunkSize = 1000000;

fid = fopen(dataFilename, 'r');

d = dir(dataFilename);
nSampsTotal = d.bytes/nChans/2;
nChunksTotal = ceil(nSampsTotal/chunkSize);

% Load voltage data
chunkInd = 1;
NIDQ = [];
while 1
  fprintf(1, 'chunk %d/%d\n', chunkInd, nChunksTotal);
  dat = fread(fid, [nChans chunkSize], '*int16');
  if ~isempty(dat)
    NIDQ = [NIDQ dat(end-1,:)]; %#ok<AGROW>
  else
    break
  end
  
%   plot(NIDQ,'r')
%   hold on
%   displacement = 1;
%   threshold = (max(NIDQ) - min(NIDQ))/50;
%   NIDQthreshold = NIDQ;
%   NIDQthreshold(NIDQthreshold<=threshold) = 0;
%   NIDQthreshold(NIDQthreshold>threshold) = 1;
%   frameOnsets = NIDQthreshold - [NIDQthreshold(1:displacement) NIDQthreshold(1:end-displacement)];
%   frameOnsets(frameOnsets<0.5) = 0;
%   frameOnsets(frameOnsets>0.5) = 1;
%   frameInd = find(logical(int8(round(frameOnsets))));
%   plot(frameInd(2:end),(frameInd(2:end)-frameInd(1:end-1))*10, '.g', 'MarkerSize',20)
%   hold off
  
  chunkInd = chunkInd+1;
end

% Detect frame times
displacement = 1;
threshold = max(NIDQ) - 0.75*int16(std(single(NIDQ))); % threshold based on mode + 0.5SD - only works if presentation time is shorter than inter trial interval
% threshold = 3000; % arbitrary threshold 
NIDQthreshold = NIDQ;

% make all points greater than threshold = 1 and all less = 0
NIDQthreshold(NIDQthreshold<=threshold) = 0;
NIDQthreshold(NIDQthreshold>threshold) = 1;

frameOnsets = NIDQthreshold - [NIDQthreshold(1:displacement) NIDQthreshold(1:end-displacement)];
frameOffsets(frameOnsets<0) = 1;
frameOnsets(frameOnsets<=0.5) = 0;
frameOnsets(frameOnsets>0.5) = 1;

% onset
frameInd = find(logical(int8(round(frameOnsets))))';
frameTimes = frameInd*(1/sr)';

% offset
frameIndoff = find(logical(int8(round(frameOffsets))))';
frameTimesoffset = frameIndoff*(1/sr)';

disp([num2str(numel(frameTimes)) ' frames detected']);

plot(NIDQ,'r');
hold on
plot(ones(size(NIDQ))*double(threshold), 'b-')
plot(frameInd, (single(max(NIDQ))+5)*ones(size(frameInd)), '.g', 'MarkerSize', 10)
for i = 1:numel(frameTimes)
text(frameInd(i), (double(max(NIDQ))+5), num2str(i))
end
hold off
 % close
