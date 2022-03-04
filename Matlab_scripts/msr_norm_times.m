function [rel_pk_tm, rel_qpk_tm, rel_INI] = msr_norm_times(pk_times, qpk_times, note_durs, INI)

% function [rel_pk_tm, rel_qpk_tm, rel_INI] = msr_norm_times(pk_times, qpk_times, note_durs, INI)
%
% S. Phelps, May 29, 2007

rel_pk_tm = pk_times./note_durs;

rel_qpk_tm = qpk_times./note_durs;

rel_INI = INI./note_durs;
