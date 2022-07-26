 function psdChannels_MS(filename, nChans, chans2ignore, sr)
% The function reads binary files and displays psd average for a range of
% frequencies for all channels in the file.
%
% Input: filename should include the full path and the extension.
%        nChans is the total number of channels in the recording.
%        chans2ignore tells which channels should be excluded from analysis.
%          Typically these are extra channels that do not correspond to
%          probe LFP recording channels.
%        sr is sampling frequency. Default is 30000 Hz.
%
% Output: figs contains figure handles.

if nargin < 4
  sr = 3e4;
end

chunkSize = 10000000;

fid = [];

d = dir(filename);
nSampsTotal = d.bytes/nChans/2;
nChunksTotal = ceil(nSampsTotal/chunkSize);

chans = ones(size(1:nChans));
chans(chans2ignore) = 0;
chans = logical(chans);

opt.typespk1 = 'c';
opt.maxFreq = 1500;
opt.winfactor = 10;
opt.freqfactor = 1.333;
opt.tapers = 5; %#ok<STRNU>

filepath = fileparts(filename);

try
  
  fid = fopen(filename, 'r');
  
  chunkInd = 1;
  while 1
    fprintf(1, 'chunk %d/%d\n', chunkInd, nChunksTotal);
    dat = fread(fid, [nChans chunkSize], '*int16');
    if ~isempty(dat) && chunkInd~=nChunksTotal
      dat = dat(chans,:);
      
      % Calculate PSDs
      for ch = 1:sum(chans)
        fprintf(1, 'channel %d/%d\n', ch, sum(chans));
        %[freq, psd] = freqDependentWindowCoherence(double(dat(ch,:)-mean(dat(ch,:))), [], 1/sr, [], opt);
        [psd, freq] = ComputeAvPowerSpectrum(double(dat(ch,:)-mean(dat(ch,:))), sr);
        if ch == 1 && chunkInd == 1
          psdTotal = zeros(sum(chans),numel(psd));
        end
        psdTotal(ch,:) = psdTotal(ch,:)+psd;
      end
      
    else
      break
    end
    chunkInd = chunkInd+1;
  end
  
  psdAverage = psdTotal./(nChunksTotal-1);
  
  % Plot full PSD for all channels separately
  for ch = 1:sum(chans)
    loglog(freq,psdAverage(ch,:))
    ylabel('Power (\muV^2/Hz)')
    xlabel('Frequency (Hz)')
    title(['PSD_ch' num2str(ch)], 'Interpreter','none');
    hgsave(gcf, [filepath filesep 'PSD_ch' num2str(ch)]);
  end
  
  % Plot frequency range average PSDs for all channels combined
  f1 = freq-1; f1(f1>0) = 0; f1(f1<0) = 1; f1 = logical(f1);
  f4 = freq-4; f4(f4>0) = 0; f4(f4<0) = 1; f4 = f4-f1; f4 = logical(f4);
  f10 = freq-10; f10(f10>0) = 0; f10(f10<0) = 1; f10 = f10-f4-f1; f10 = logical(f10);
  f15 = freq-15; f15(f15>0) = 0; f15(f15<0) = 1; f15 = f15-f10-f4-f1; f15 = logical(f15);
  f30 = freq-30; f30(f30>0) = 0; f30(f30<0) = 1; f30 = f30-f15-f10-f4-f1; f30 = logical(f30);
  f60 = freq-60; f60(f60>0) = 0; f60(f60<0) = 1; f60 = f60-f30-f15-f10-f4-f1; f60 = logical(f60);
  f100 = freq-100; f100(f100>0) = 0; f100(f100<0) = 1; f100 = f100-f60-f30-f15-f10-f4-f1; f100 = logical(f100);
  f200 = freq-200; f200(f200>0) = 0; f200(f200<0) = 1; f200 = f200-f100-f60-f30-f15-f10-f4-f1; f200 = logical(f200);
  f500 = freq-500; f500(f500>0) = 0; f500(f500<0) = 1; f500 = f500-f200-100-f60-f30-f15-f10-f4-f1; f500 = logical(f500);
  f1500 = freq-500; f1500(f1500>0) = 0; f1500(f1500<0) = 1; f1500 = f1500-f500-f200-100-f60-f30-f15-f10-f4-f1; f1500 = logical(f1500);
  fBand = [f1; f4; f10; f15; f30; f60; f100; f200; f500; f1500];
  fBandTitle = {'0p1_1Hz'; '1_4Hz'; '4_10Hz'; '10_15Hz'; '15_30Hz'; '30_60Hz'; '60_100Hz'; '100_200Hz';'200_500Hz'; '500_1500Hz'};
  chBandPSDall = zeros(size(fBand,1),sum(chans));
  for b = 1:size(fBand,1)
    chBandPSD = zeros(1,sum(chans));
    for ch = 1:sum(chans)
      chBandPSD(ch) = sum(psdAverage(ch,fBand(b,:)))/sum(fBand(b,:));
    end
    chBandPSDall(b,:) = chBandPSD;
    plot(chBandPSD,1:numel(chBandPSD))
    xlabel('Power (\muV^2/Hz)')
    ylabel('Channel #')
    title(['Average PSD in band ' fBandTitle{b}], 'Interpreter','none');
    hgsave(gcf, [filepath filesep 'PSDaverage_band' fBandTitle{b}]);
  end
  
  save([filepath filesep 'chanPSDs.mat'], 'freq', 'psdAverage', 'chBandPSDall');
  
  close(gcf)
  fclose(fid);
  
catch me
  
  if ~isempty(fid)
    fclose(fid);
  end
  
  rethrow(me)
  
end
