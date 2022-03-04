function y = bandpass_sks(song,samp_freq)

% define default sample frequency
if nargin<2
    samp_freq = 195312.5/2;
end

y = bandpass(song, [9000 44000], samp_freq);