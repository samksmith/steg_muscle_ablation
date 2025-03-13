# steg_muscle_ablation
Repository for data, scripts, and analysis of muscle ablation study in singing mice (Scotinomys teguina). Histological slide images can be found on figshare project: Cricothyroid Muscle Ablation in Alston's singing mouse, (Scotinomys teguina).


Detecting songs in raw audio files, extracting songs, measuring song and notes variables are done using the following matlab scripts:

auto_detect_view_sks.m - Detect singing mouse songs in a continuous recording (20 min long) and plot spectrograms for each song.

auto_extract_sks.m - Detect songs in a continuous recording (20 min long) and extract each song.

msr_multiple.m - main script that measures song variables and outputs a table of whole song measures and note matrices for each song for a list of songs.

Msr_multiple calls the following:

msr_note_times.m - Returns note starts, ends, duration, and internote interval. Calls find_long_call.m and remove_subnotes.m

msr_whole_call.m - Returns number of notes, call length, call dominant frequency, entropy, peak amplitude, RMS, and peak power.

msr_all_notes_arcL.m - returns matrix of attributes for each note including arc length and chord length. Calls msr_call_FM_arcL.m and msr_note_DFs.m.

msr_curve.m - fits curves to song notes.

note_matrices/ - folder contains note measurement data.

muscle_ablation_muscle_measurements.csv - contains muscle area estimates.

Note matrices and muscle data are analyzed using the following R scripts:

prep_note_matrices.R - combines all note matrices for each individual, filters and QCs arc length measurements.

nmAL_analysis.R - calculates mean normalized arc length, calculates volume estimates for muscles, and generates models.

