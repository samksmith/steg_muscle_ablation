function [pk_amp, pk_time, qtr_pk_amp, qtr_pk_time, note_rms] = msr_note_amps(note, samp_rate)

% function [pk_amp, pk_time, qtr_pk_amp, qtr_pk_time, note_rms] = msr_note_amps(note, samp_rate)

if nargin < 2
    samp_rate = 195312.5;
end


[note_length, c] = size(note);

[pk_amp, pk_time] = max(note);

pk_time = pk_time/samp_rate;


qtr_note_length = round(note_length/4);

[qtr_pk_amp, qtr_pk_time] = max(note(1:qtr_note_length));

qtr_pk_time = qtr_pk_time/samp_rate;


note_sq = note .* note;

note_rms = sqrt(sum(note_sq)/note_length);