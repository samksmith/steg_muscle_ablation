# steg_muscle_ablation
Repository for data, scripts, and analysis of muscle ablation study in singing mice (Scotinomys teguina)

Matlab scripts:
auto_detect_view_sks.m - Detect singing mouse songs in a continuous recording (20 min long) and plot spectrograms for each song.
auto_extract_sks.m - Detect songs in a continuous recording (20 min long) and extract each song.
msr_multiple.m - main script that measures song variables and outputs a table of whole song measures and note matrices for each song for a list of songs.

Msr_multiple calls the following:
msr_note_times.m - Returns note starts, ends, duration, and internote interval. Calls find_long_call.m and remove_subnotes.m
msr_whole_call.m - Returns number of notes, call length, call dominant frequency, entropy, peak amplitude, RMS, and peak power.
msr_all_notes_arcL.m - returns matrix of attributes for each note including arc length and chord length. Calls msr_call_FM_arcL.m and msr_note_DFs.m.
msr_curve.m - fits curves to song notes.

