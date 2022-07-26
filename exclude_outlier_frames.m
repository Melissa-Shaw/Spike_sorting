function [frameTimes, frames_removed] = exclude_outlier_frames(frameTimes,buffer)
frames_removed = [];
for i=1:numel(frameTimes)-1
  if i==1
    frame_before = -buffer;
  else
    frame_before = frameTimes(i-1);
  end
  frame_after = frameTimes(i+1);
  if frameTimes(i)>frame_before+buffer && frameTimes(i)<frame_after-buffer
    frames_removed = [frames_removed i];
  end
end
frameTimes(frames_removed) = [];
end