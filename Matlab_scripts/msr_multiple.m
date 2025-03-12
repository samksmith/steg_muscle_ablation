%The routine provides one output table describing all scalar descriptors of a call.
% This is a rather long list, including the parameters for curves that describe how note variables change over the course
% of a call. To clarify, the variable names are included in the table. Labels ending in "a", "b" and "c" correspond
% to the coefficients in a quadratic curve fit for the changes in that variable over the course of a curve. Labels ending in
% "res" are the average residuals from the corresponding curves. The table is written out as a comma-delimited tab file. 
% a matrix of measurements for each note is also generated for each call
% ("all_notes_matrix"), but not saved into a master table. There is an
% option for saving the note_matrix file that is commented out by default.
% 
% The routine calls, directly or indirectly, nearly all of the routines in the "call analysis tools" folder. 
% The routines it calls directly are MSR_CURVE, MSR_CALL_AMPS, MSR_WHOLE_CALL, MSR_ALL_NOTES, 
%
% Please note: if more than one call is present in the data, it will analyze the
% longest call.
% 
% Original script: S. Phelps 6-03-07. Adapted to automatically take
% measurements of multiple song files by S. Smith 03-03-22.

%% Initializations
% Change these default values
clear;clc;close all;
f_group = dir('41-13-1_1*.F32'); % list all the recordings from a particular individual
nFiles = length(f_group); % how many files are there

fig_on = 0; % should figures be generated and output?
INI_max = 200; % what is the max internote interval to consider still same song? (ms)
curve_meth = 'q';
time_code = 'note';
reset = 10; 
threshold = 8;
samp_rate = 195312.5/2;

% initialize a table that will hold the data for all the songs of a
% particular individual
% initialize size of table, the types of variables going in and the
% variable names
size_table = [nFiles 96];
varTypes = ["double","double","double","double","double","double","double","double","double","double",...
    "double","double","double","double","double","double","double","double","double","double","double","double",...
    "double","double","double","double","double","double","double","double","double","double","double","double",...
    "double","double","double","double","double","double","double","double","double","double","double","double",...
    "double","double","double","double","double","double","double","double","double","double","double","double",...
    "double","double","double","double","double","double","double","double","double","double","double","double",...
    "double","double","double","double","double","double","double","double","double","double","double","double",...
    "double","double","double","double","double","double","double","double","double","double","double","double","double","string"];
varNames = ["note_num","call_length","call_DF","entropy","pk_amp","rms","pk_power","note.starts_a",...
    "note.starts_b","note.starts_c","note.starts_res","note.ends_a","note.ends_b","note.ends_c","note.ends_res",...
    "note.durs_a","note.durs_b","note.durs_c","note.durs_res","INI_a","INI_b","INI_c","INI_res","pk.amps_a","pk.amps_b",...
    "pk.amps_c","pk.amps_res","pk.times_a","pk.times_b","pk.times_c","pk.times_res","qpk.amps_a","qpk.amps_b","qpk.amps_c",...
    "qpk.amps_res","qpk.times_a","qpk.times_b","qpk.times_c","qpk.times_res","note.rmss_a","note.rmss_b","note.rmss_c",...
    "note.rmss_res","rel.pk.amps_a","rel.pk.amps_b","rel.pk.amps_c","rel.pk.amps_res","rel.qpk.amps_a","rel.qpk.amps_b",...
    "rel.qpk.amps_c","rel.qpk.amps_res","rel.note.rmss_a","rel.note.rmss_b","rel.note.rmss_c","rel.note.rmss_res",...
    "rel.pk.tm_a","rel.pk.tm_b","rel.pk.tm_c","rel.pk.tm_res","rel.qpk.tm_a","rel.qpk.tm_b","rel.qpk.tm_c","rel.qpk.tm_res",...
    "rel.INI_a","rel.INI_b","rel.INI_c","rel.INI_res","max.Hz_a","max.Hz_b","max.Hz_c","max.Hz_res","min.Hz_a","min.Hz_b",...
    "min.Hz_c","min.Hz_res","FMa_a","FMa_b","FMa_c","FMa_res","FMb_a","FMb_b","FMb_c","FMb_res","FMc_a","FMc_b","FMc_c",...
    "FMc_res","FM.resid_a","FM.resid_b","FM.resid_c","FM.resid_res","note.DF_a","note.DF_b","note.DF_c","note.DF_res","call_ID"];
% actually make the empty table
calldata_table = table('Size',size_table,'VariableTypes',varTypes,'VariableNames',varNames);

for f = 1:nFiles % take each file and measure all relevant parameters
    f_group(f).name % chosing which file to open
    fid = fopen(f_group(f).name,'r'); % opening file
    [call, c1] = fread(fid,'float32'); % reading in data
    fclose('all');
    call_id = f_group(f).name; % record name of song
    
    % Measure note number, onset and offset, trim "call" to just those samples within call.
    [note_starts, note_ends, note_durs, INI] = msr_note_times(call, samp_rate, threshold, reset, INI_max);
    
    % Measure attributes of whole call. Each of these measures is a scalar.
    [note_num, call_length, call_DF, entropy, pk_amp, rms, pk_power] = msr_whole_call(call, samp_rate, threshold, reset, note_starts, note_ends, INI_max);
   
    % Measure attributes of notes. The resulting matrix has notes in a column, variable values in each column.
    [all_notes_matrix, note_labels] = msr_all_notes_arcL(call, samp_rate, note_starts, note_ends);

    % Measure and plot curves.
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
    
    % make an array of all the results from the song
    call_stats = [note_num, call_length, call_DF, entropy, pk_amp, rms, pk_power,curves];
    % double --> cell array so can add to the master table
    call_stats = num2cell(call_stats);
    call_stats{1,96} = call_id; % add in song ID to the data array
    % add the call stats to the master table of results for the individual
    calldata_table(f,:) = call_stats;
    
    % if you want to save the matrix of each note measurement, 
    % make sure the following lines are uncommented 
    % name for file writing out
    nameRoot = call_id(1:end-4); % take the time stamp and song number from original file
    filename_notesmatr = strcat(nameRoot,'_notematrix.txt'); % make name for notesmatrix file
    % variable names
    names = ["note.starts","note.ends","note.durs","INI","pk.amps","pk.times",...
        "qpk.amps","qpk.times","note.rmss","rel.pk.amps","rel.qpk.amps","rel.note.rmss",...
        "rel.pk.tm","rel.qpk.tm","rel.INI","max.Hz","min.Hz","FMa","FMb","FMc","FM.resid","note.DF"];
    % make a table for notes matrix (currently is saved as an array)
    notes_matrix = array2table(all_notes_matrix,'VariableNames',names);
    % actually write the table out
    writetable(notes_matrix,filename_notesmatr,'WriteRowNames',true);
end
% filename for master data out table
filename = '41-13-1_presurgery.txt';
% write the master table of results out
writetable(calldata_table,filename,'WriteRowNames',true);

