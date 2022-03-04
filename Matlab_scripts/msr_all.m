function [call_stats, call_stat_labels, all_notes_matrix, note_labels] = msr_all(call, samp_rate, threshold, reset, time_code, curve_meth, INI_max, fig_on, call_id)

% [call_stats, call_stat_labels, all_notes_matrix, note_labels] = msr_all(call, samp_rate, threshold, reset, time_code, curve_meth, INI_max, fig_on, call_id)
%
% This routine takes 8 input arguments, the first being a call to analyze. If more than one call is present in the data, it will analyze the
% longest call. The latter arguments can be omitted for default values. Default values are
%
% if nargin < 8, fig_on = 1; end
% if nargin < 7, INI_max = 200; end
% if nargin < 6, curve_meth = 'q'; end
% if nargin < 5, time_code = 'note'; end
% if nargin < 4, reset = 10; end
% if nargin < 3, threshold = 8; end
% if nargin < 2, samp_rate = 195312.5; end
%
% If the fig_on input is true (=1), the routine will make a graph of all curves for the respective variables, as well as a
% spectrogram. 
%
% The routine provides three output variables. The first is a row vector describing all scalar descriptors of a call.
% This is a rather long list, including the parameters for curves that describe how note variables change over the course
% of a call. To clarify, it also lists the labels of each of these variables. Labels ending in "a", "b" and "c" correspond
% to the coefficients in a quadratic curve fit for the changes in that variable over the course of a curve. Labels ending in
% "res" are the average residuals from the corresponding curves. The third output variable is a matrix of variables measured for
% individual notes within a call (notes in rows, variables in columns) on which the curves were based. The final output is
% a list of labels for this matrix. Note that both label lists have the variables in rows and letters of the labels in columns.
% In contrast, the data outputs have variables in columns and notes in columns. This is due to the way Matlab treats text.
% The text will need to be transposed before being used as a labels in a spreadsheet or text file.
% 
% The routine calls, directly or indirectly, nearly all of the routines in the "call analysis tools" folder. 
% The routines it calls directly are MSR_CURVE, MSR_CALL_AMPS, MSR_WHOLE_CALL, MSR_ALL_NOTES, 
% 
% S. Phelps 6-03-07.

% Set defaults
%---------
if nargin < 9, call_id = 'unnamed call'; end
    if nargin < 8, fig_on = 1; end
if nargin < 7, INI_max = 200; end
if nargin < 6, curve_meth = 'q'; end
if nargin < 5, time_code = 'note'; end
if nargin < 4, reset = 10; end
if nargin < 3, threshold = 8; end
if nargin < 2, samp_rate = 195312.5; end



% Measure note number, onset and offset, trim "call" to just those samples within call.
%----------------------------------------------
[note_starts, note_ends, note_durs, INI] = msr_note_times(call, samp_rate, threshold, reset, INI_max);
%[note_num,c] = size(note_ends);

%a = ceil(note_starts(1,1)*samp_rate/1000);
%b = floor(note_ends(note_num,1)*samp_rate/1000);
%call = call(a:b);

%note_ends = note_ends - note_starts(1,1);
%note_starts = note_starts - note_starts(1,1) + 1000/samp_rate;


% Measure attributes of whole call. Each of these measures is a scalar.
%----------------------------------------------
[note_num, call_length, call_DF, entropy, pk_amp, rms, pk_power] = msr_whole_call(call, samp_rate, threshold, reset, note_starts, note_ends, INI_max);


% Measure attributes of notes. The resulting matrix has notes in a column, variable values in each column.
%----------------------------------------------
[all_notes_matrix, note_labels] = msr_all_notes(call, samp_rate, note_starts, note_ends);


% Measure and plot curves.
%-------------------------
[r, var_num] = size(all_notes_matrix);

% Define dimensions of graph, plot spectrogram, if figure is being made
if (fig_on)
    figure, title(call_id);
    rows = ceil(var_num/4)+1;
    if rows>3, rows = 3; end
    columns = 4;
    subplot(rows,1,1), specgram(call, 256, samp_rate), caxis([-100, 35]); 
end

% Calculate curves for each column of note stats
sp=5; %position of subplot skips first row for specgram.

for i =1:var_num
    
    j=(i-1)*4 +1;
    var_name = (note_labels(i,:)); 
    
    if all_notes_matrix(note_num,i) == 0 % Redefine note number for INI measures to avoid zero filler
        note_num=r-1;
    end
    
    if (fig_on) % Activate subplot window for msr_curve routine.
        if (i == 9)|(i == 21) 
            figure, title(call_id); 
            sp=1;
        end    
        subplot(rows,columns,sp);
        sp=sp+1;
    end
    
    [curves(1,j:j+2), curves(1,j+3)] = msr_curve(all_notes_matrix(1:note_num,i), note_starts(1:note_num), note_ends(1:note_num), time_code, curve_meth, fig_on, var_name);
    
    if i ==1    
        curve_labels = char(strcat(note_labels(i,:),'_a'), strcat(note_labels(i,:),'_b'), strcat(note_labels(i,:),'_c'), strcat(note_labels(i,:),'_res'));
    else
        curve_label_temp = char(strcat(note_labels(i,:),'_a'), strcat(note_labels(i,:),'_b'), strcat(note_labels(i,:),'_c'), strcat(note_labels(i,:),'_res'));
        curve_labels = char(curve_labels, curve_label_temp);
    end
    
    note_num = r;
    
    if fig_on % Label curve drawn by msr_curve
        if (sp<10), xlabel (''), end
        if (i == var_num)
            a = axis;
            x = 1.25*((a(1,2)-a(1,1)) + a(1,1));
            y = (a(1,4) - a(1,3))/2 + a(1,3);
            text(x,y,char(strcat('note number = ', num2str(note_num)), strcat('length = ', num2str(call_length)), strcat('DF = ', num2str(call_DF)), strcat('peak amp = ', num2str(pk_amp)), strcat('RMS amp = ', num2str(rms)))), xlabel('time(sec)'), ylabel('frequency(hz)'); 
        end 
    end
    if fig_on == 2 % Save each figure with a unique name
        if i == 8
            fig_file_name = strcat(call_id, '_A');
            saveas(gcf, fig_file_name, 'jpeg'); end
        if i == 20
            fig_file_name = strcat(call_id, '_B');
            saveas(gcf, fig_file_name, 'jpeg');
            close; end
        if i == 22
            fig_file_name = strcat(call_id, '_C');
            saveas(gcf, fig_file_name, 'jpeg');
            close; end
    end
end

% Concatenate data and labels for output
%---------------------------------------

call_stats = [note_num, call_length, call_DF, entropy, pk_amp, rms, pk_power, curves];

call_stat_labels = char('note_num', 'call_length', 'call_DF', 'entropy', 'pk_amp', 'rms', 'pk_power');
call_stat_labels = char(call_stat_labels, curve_labels);
