function [all_notes_matrix, labels] = msr_all_notes(call, samp_rate, note_starts, note_ends)

% [all_notes_matrix, labels] = msr_all_notes(call, samp_rate, note_starts, note_ends)

if nargin<2, samp_rate = 195312.5; end

[note_num, c] = size(note_starts);

% Begin with amplitude modulated measures of notes
%-------------------------------------------------
[pk_amps, pk_times, qpk_amps, qpk_times, note_rmss] = msr_all_note_amps(call, note_starts, note_ends, samp_rate);

note_durs = note_ends - note_starts;
INI = note_starts(2:note_num) - note_ends(1:note_num-1);

[rel_pk_amps, rel_qpk_amps, rel_note_rmss] = msr_norm_amps(pk_amps, qpk_amps, note_rmss);

INI(note_num,1) = 0; % Note that the inter-note interval needs a zero entered for its last position before writing to matrix to keep number of rows constant.
[rel_pk_tm, rel_qpk_tm, rel_INI] = msr_norm_times(pk_times, qpk_times, note_durs, INI);


% Next look at frequency modulation
%----------------------------------
[notes_max_Hz, notes_min_Hz, notes_FM, notes_resid] = msr_call_FM(call, samp_rate, 256, 'q', note_starts, note_ends); 
% dimensions of notes_FM are 1r x 3c if quadratic function -- which is default.

[notes_DF] = msr_note_DFs(call, note_starts, note_ends, samp_rate);

% Now compile various measures into a single matrix, in which the note number is the number of rows



all_notes_matrix = [note_starts, note_ends, note_durs, INI, pk_amps, pk_times, qpk_amps, qpk_times, note_rmss, rel_pk_amps, rel_qpk_amps, rel_note_rmss, rel_pk_tm, rel_qpk_tm, rel_INI, notes_max_Hz, notes_min_Hz, notes_FM, notes_resid, notes_DF]; %Note this assumes FM fit by quadratic, which is default.

labels = char('note.starts', 'note.ends', 'note.durs', 'INI', 'pk.amps', 'pk.times', 'qpk.amps', 'qpk.times', 'note.rmss', 'rel.pk.amps', 'rel.qpk.amps', 'rel.note.rmss', 'rel.pk.tm', 'rel.qpk.tm', 'rel.INI', 'max.Hz', 'min.Hz', 'FMa', 'FMb', 'FMc', 'FM.resid', 'note.DF');
% note labels has each label in a row. After saving data, or before writing to common file, will need to transpose to match matrices. For now this is convenient.
