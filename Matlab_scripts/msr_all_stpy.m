function [call_stpy, note_stpy, note_diff, note_labels, all_notes_matrix] = msr_all_stpy(call, samp_rate, threshold, reset, INI_max, fig_on, call_id)

% [call_stpy, note_stpy, note_diff, note_labels, all_notes_matrix] = msr_all_stpy(call, samp_rate, threshold, reset, INI_max, fig_on, call_id)
%
% Stereoptypy in this routine is measured by the (dis)similarity of adjacent notes. We take the difference between two notes in any parameter
% and divide this by the mean of the parameter across the two notes. To get a call-wide measure of stereotypy, we averagethese across note pairs.
%
% This routine takes 8 input arguments, the first being a call to analyze. If more than one call is present in the data, it will analyze the
% longest call. The latter arguments can be omitted for default values. Default values are
%
% if nargin < 7, call_id = 'unnamed call'
% if nargin < 6, fig_on = 0
% if nargin < 5, INI_max = 200
% if nargin < 4, reset = 10
% if nargin < 3, threshold = 8
% if nargin < 2, samp_rate = 195312.5
%
% If the fig_on input is true (=1), the routine will make a graph of note stereotypies for the respective variables, as well as a
% spectrogram. 
%
% The routine provides five output variables. CALL_STPY is a row vector describing the average stereotypy measures for each parameter.
% NOTE_STPY provides stereotypy measures for each pair of notes (rows), and each parameter (columns). NOTE_DIFFS provides analogous data
% that have not been weighted by note means. NOTE_LABELS provides a list of parameters measured across notes. Notice that the label list has
% the variables in rows and letters of the labels in columns. In contrast, the data outputs have variables in columns and notes in columns. 
% This is due to teh way Matlab treats text. The text will need to be transposed before being used as a labels in a spreadsheet or text file.
% The last output is the raw note measures the code uses to compute stereotypy.
% 
% The routine calls, directly or indirectly, nearly all of the routines in the "call analysis tools" folder. 
% The routines it calls directly are MSR_CALL_AMPS, MSR_ALL_NOTES, MSR_NOTE_TIMES.
% 
% S. Phelps 6-03-07.

% Set defaults
%---------
if nargin < 7, call_id = 'unnamed call'; end
if nargin < 6, fig_on = 0; end
if nargin < 5, INI_max = 200; end
if nargin < 4, reset = 10; end
if nargin < 3, threshold = 8; end
if nargin < 2, samp_rate = 195312.5; end



% Measure note number, onset and offset, trim "call" to just those samples within call.
%----------------------------------------------
[note_starts, note_ends, note_durs, INI] = msr_note_times(call, samp_rate, threshold, reset, INI_max);


% Measure attributes of notes. The resulting matrix has notes in a column, variable values in each column.
%----------------------------------------------
[all_notes_matrix, note_labels] = msr_all_notes(call, samp_rate, note_starts, note_ends);


% Measure stereotypy, assumes at least two notes
%-------------------------
[note_num, var_num] = size(all_notes_matrix); 

for i=2:note_num
    note_diff(i-1,:) = all_notes_matrix(i,:) - all_notes_matrix(i-1,:);
    note_stpy(i-1,:) = abs(2*note_diff(i-1,:)./(all_notes_matrix(i,:) + all_notes_matrix(i-1,:)));
end

call_stpy = mean(note_stpy);


% Plot curves.
%-------------------------
% Define dimensions of graph, plot spectrogram, if figure is being made
if (fig_on)
    figure, title(call_id);
    rows = ceil(var_num/4)+1;
    if rows>3, rows = 3; end
    columns = 4;
    subplot(rows,1,1), specgram(call, 256, samp_rate), caxis([-100, 35]); 
    
    % Show data for each column of note stats
    sp=5; %position of subplot skips first row for specgram.
    
    for i =1:var_num
        [note_num, c] = size(all_notes_matrix);     
        j=(i-1)*4 +1;
        var_name = (note_labels(i,:)); 
        
        if all_notes_matrix(note_num,i) == 0 % Redefine note number for INI measures to avoid zero filler
            note_num=note_num-1;
        end
        
        % Activate subplot window
        if (i == 9)|(i == 21) 
            figure, title(call_id); 
            sp=1;
        end    
        subplot(rows,columns,sp);
        
        %plot and label each variable
        T = 1:1:(note_num-1);
        T = T';  
        plot(T, note_stpy(1:(note_num-1),i), 'o'), xlabel(['note pair']), ylabel(['relative difference']), title(note_labels(i,:));
        
        sp=sp+1;
    end
end

% Concatenate data and labels for output
%---------------------------------------
%call_stats = [note_num, call_length, call_DF, entropy, pk_amp, rms, pk_power, curves];
%parameter_labels = char('note_num', 'call_length', 'call_DF', 'entropy', 'pk_amp', 'rms', 'pk_power');
%parameter_labels = char(call_stat_labels, curve_labels);
