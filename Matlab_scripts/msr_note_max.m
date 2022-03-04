function[note_max, place_max] = msr_note_max(file_name, samp_freq)

% defining my terms; assigning variables for future use
song = read_songs(file_name);
[note_starts, note_ends, note_durs] = msr_note_times(song);
num_notes = numel(note_starts);
note_max = 0;
place_max = 0;

% dealing with argument fuzziness
if nargin < 2
    samp_freq = 195312.5;
end

for i=1:num_notes
    this_start = note_starts(i,1); 
    this_end = note_ends(i,1); 
    note_sum = this_end - this_start;
    
    % write error message here in case note_sum =/= note_dur (i,1)
    if note_sum ~= note_durs(i,1)
        error('Error: check to make sure that starts and ends of notes are matching properly')
    end
    
    first_sample = round((this_start * samp_freq) / 1000);
    last_sample = round((this_end * samp_freq) / 1000);
    this_note = song(first_sample:last_sample,1);
    ab_this = abs(this_note); 
    [this_note_max, place_peak] = max(ab_this);
    this_place_max = place_peak / samp_freq;
    note_max(i,1) = this_note_max;
    place_max(i,1) = this_place_max;
end