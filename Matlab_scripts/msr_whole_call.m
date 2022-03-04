function [note_num, call_length, call_DF, entropy, pk_amp, rms, pk_power] = msr_whole_call(call, samp_rate, threshold, reset, note_starts, note_ends, INI_max);

% function [note_num, call_length, call_DF, entropy, pk_amp, rms, pk_power] = msr_whole_call(call, samp_rate, threshold, reset, note_starts, note_ends, INI_max);
%
% Defaults:
% INI_max = 200
% reset = 10
% threshold = 8
% samp_rate = 195312.5
%
% If needed, this routine will calculate onset and offset of notes by passing data to msr_note_times.


% Assign defaults
%----------------
if nargin < 7, INI_max = 200; end
if nargin < 4, reset = 10; end
if nargin < 3, threshold = 8; end
if nargin < 2, samp_rate = 195312.5; end

% If needed, calculate onset and offset of notes
if nargin<5
    [note_starts, note_ends, note_durs, INI] = msr_note_times(call, samp_rate, threshold, reset, INI_max);
end

if nargin <1, error('Call data required but not entered.'), end



% Calculate note number and call length (sec)
%--------------------------------------------
[note_num, c] = size(note_starts);
call_length = (note_ends(note_num,1)-note_starts(1,1))/1000; %Since note boundaries are in msec, diving by 1000 gives call duration in sec

a = ceil(note_starts(1,1)*samp_rate/1000); % convert from msec to samples
b = floor(note_ends(note_num,1)*samp_rate/1000);

% Calculate call dominant frequency
%----------------------------------
[Pxx, F] = pwelch(call(a:b),[],[],[], samp_rate);
[pk_power, pk_index] = max(Pxx);
call_DF = F(pk_index);

% Calculate call entropy
%-----------------------
entropy = mean(log(Pxx/mean(Pxx)));

% Calculate RMS, peak amplitude
%------------------------------
pk_amp = max(abs(call(a:b)));
rms = sqrt(mean(call(a:b).*call(a:b)));


