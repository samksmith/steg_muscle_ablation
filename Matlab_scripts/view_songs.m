function [songs] = view_songs(file_name, samp_freq)

fid = fopen(file_name, 'r');

if fid == -1
      error('Error: check that path has been set to find appropriate file.')
end

songs = fread(fid, 'float32');
[song_length,~] = size(songs);

% define default sample frequency
if nargin<2
    samp_freq = 195312.5;
end

%
% -------------
% Now plot songs in spectrogram form
% -------------

figure()
specgram(songs, 512, samp_freq);
caxis([-100 20])

%Uncomment these if you wish to save the figure. 
%filename = strcat(file_name,'1.fig')
%savefig(filename)

%Now plot calls in oscillogram form

%define X axis
time_between = 1/samp_freq;
%define Y axis
y = songs;

x = 0:time_between:(song_length-1)*time_between;
title(file_name);
%length(x) == length(y);

figure()
plot(x,y)