function [note_DFs] = msr_note_DFs(call, note_starts, note_ends, sampling_rate);

% function [note_DFs] = msr_note_DFs(call, note_starts, note_ends, sampling_rate);
% This works!

if nargin<4
    sampling_rate = 195312.5;
end

[rs, cs] = size(note_starts);
[re, ce] = size(note_ends);

if rs ~= re
    error('note starts and ends of different dimensions')
end

if cs ~= ce
    error('note starts and ends of different dimensions')
end

for i=1:rs
    note_start = round(note_starts(i,1)*sampling_rate/1000); %assuming note info stored as msec
    note_end = round(note_ends(i,1)*sampling_rate/1000);
    [Pxx, F] = pwelch(call(note_start:note_end),[],[],[], sampling_rate);
    
    [pk_power(i,1), pk_index] = max(Pxx);
    note_DFs(i,1) = F(pk_index);
end
