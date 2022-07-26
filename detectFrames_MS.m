function [frameTimes, frameInd] = detectFrames_MS(dataFilename, nChans, sr)
% A helper function for eyeAnalysis script. It finds the times and indices
% of every consecutive peak in the voltage data recorded in an extra Open
% Ephys channel. These peaks correspond to the eye-tracker video frame
% refresh instances.
% Input: dataFilename - a string with the name of the data file (in dat or
%                       bin formats).
%        nChans - a total number of channels in the recording. The voltage
%                 channel is supposed to be the last one.
%        sr - a sampling rate.
% Output: frameTimes and frameInd.

chunkSize = 1000000;

fid = fopen(dataFilename, 'r');

d = dir(dataFilename);
nSampsTotal = d.bytes/nChans/2;
nChunksTotal = ceil(nSampsTotal/chunkSize);

% Load voltage data
chunkInd = 1;
eyeData = [];
while 1
  fprintf(1, 'chunk %d/%d\n', chunkInd, nChunksTotal);
  dat = fread(fid, [nChans chunkSize], '*int16');
  if ~isempty(dat)
    eyeData = [eyeData dat(end,:)]; %#ok<AGROW>
  else
    break
  end
  
%   plot(eyeData,'r')
%   hold on
%   displacement = 1;
%   threshold = (max(eyeData) - min(eyeData))/3;
%   eyeDataTh = eyeData;
%   eyeDataTh(eyeDataTh<=threshold) = 0;
%   eyeDataTh(eyeDataTh>threshold) = 1;
%   frameOnsets = eyeDataTh - [eyeDataTh(1:displacement) eyeDataTh(1:end-displacement)];
%   frameOnsets(frameOnsets<0.5) = 0;
%   frameOnsets(frameOnsets>0.5) = 1;
%   frameInd = find(logical(int8(round(frameOnsets))));
%   plot(frameInd(2:end),(frameInd(2:end)-frameInd(1:end-1))*10, '.g', 'MarkerSize',20)
%   hold off
  
  chunkInd = chunkInd+1;
end

% Detect frame times
displacement = 1;
threshold = (max(eyeData) - min(eyeData))/5;
eyeDataTh = eyeData;
eyeDataTh(eyeDataTh<=threshold) = 0;
eyeDataTh(eyeDataTh>threshold) = 1;
frameOnsets = eyeDataTh - [eyeDataTh(1:displacement) eyeDataTh(1:end-displacement)];
frameOnsets(frameOnsets<=0.5) = 0;
frameOnsets(frameOnsets>0.5) = 1;

frameInd = find(logical(int8(round(frameOnsets))))';
frameTimes = frameInd*(1/sr)';

plot(eyeData,'r')
hold on
plot(frameInd, double(max(eyeData)+1)*ones(size(frameInd)), '.g', 'MarkerSize', 20)


