function [pk_amps, pk_times, qtr_pk_amps, qtr_pk_times, note_rmss] = msr_all_note_amps(call, note_starts, note_ends, samp_rate)
%function [pk_amps, pk_times, qtr_pk_amps, qtr_pk_times, note_rmss] = msr_all_note_amps(call, note_starts, note_ends, samp_rate)

if nargin<4
    samp_rate = 195312.5;
end

[note_num, c] = size(note_starts);

a = round(note_starts*samp_rate/1000);
b = round(note_ends*samp_rate/1000);

for i=1:note_num
    [pk_amps(i,1), pk_times(i,1), qtr_pk_amps(i,1), qtr_pk_times(i,1), note_rmss(i,1)] = msr_note_amps(call(a(i):b(i)), samp_rate);
end