function [notes_max_Hz, notes_min_Hz, notes_FM, notes_resid] = msr_call_FM(call, sampling_rate, spec_resol, method, note_starts, note_ends, threshold, reset, INI_max);

% function [note_max_Hz, note_min_Hz, notes_FM, note_resid] =
% msr_call_FM(call, sampling_rate, spec_resol, method, note_starts, note_ends, threshold, reset, INI_max);
%
% INI_max = 200
% reset = 10
% threshold = 8
% method = 'q'
% spec_resol = 256
% sampling_rate = 195312.5
% 
% Can measure note starts and stops if omitted from input arguments.

if nargin <1, error('Call data required but not entered.'), end


% Check inputs and assign defaults
%---------------------------------
if nargin<9, INI_max = 200; end
if nargin<8, reset = 10; end
if nargin<7, threshold = 8; end
if nargin<4, method = 'q'; end
if nargin<3, spec_resol = 256; end
if nargin<2, sampling_rate = 195312.5; end
if nargin <1, error('Call data required but not entered.'), end

% If needed, calculate onset and offset of notes
if nargin<6
    [note_starts, note_ends, note_durs, INI] = msr_note_times(call, sampling_rate, threshold, reset, INI_max);
end

[note_num, c] = size(note_starts);

a = round(note_starts*sampling_rate/1000); % convert from msec to samples
b = round(note_ends*sampling_rate/1000);

a = a - spec_resol/2; % For purposes of measuring frequency parameter, need to expand window so that frequency analysis window can be centered at beginning and end of note. 
a(a<0) = 1; %a should never be negative! 
b = b + spec_resol/2;

for i=1:note_num
    [notes_max_Hz(i,1), notes_min_Hz(i,1), notes_FM(i,:), notes_resid(i,1)] = msr_note_Hz(call(a(i):b(i)), sampling_rate, spec_resol, method);
end