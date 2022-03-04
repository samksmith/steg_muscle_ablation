function [rel_pk_amps, rel_qpk_amps, rel_note_rmss] = msr_norm_amps(pk_amps, qpk_amps, note_rmss)

% function [rel_pk_amps, rel_qpk_amps, rel_note_rmss] = msr_norm_amps(pk_amps, qpk_amps, note_rmss)
%
% S. Phelps, May 28, 2007

rel_pk_amps = pk_amps/max(pk_amps);

rel_qpk_amps = qpk_amps/max(pk_amps);

rel_note_rmss = note_rmss/max(note_rmss);