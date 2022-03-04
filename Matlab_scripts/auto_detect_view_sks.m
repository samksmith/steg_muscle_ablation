%% Script to automatically detect S.teg Songs from continuous recordings
% and find the starts and stops of individual songs. Show each song chunk.

%% Initializations
% Change these default values if needed
clear;clc;close all;
f_group = dir('53-2-1*f32'); % list all the recordings from a particular chamber
indName = "53-2-1"; % name of the individual recorded
% Sorting files by order recorded
[~,idx] = sort([f_group.datenum]);

% how many files are there
nFiles = length(f_group); % how many files are there
Fs = 195312.5/2; % sampling rate of the recordings
FsN = 1000; % new sampling  rate for downsampling recordings
lowF = 12e3;highF = 43e3; %range of note frequency
syllSmooth = .01; %smoothing of syllables, in seconds
[a,b] = rat(Fs/FsN); % ratio of original to downsampled sampling rates
gap_thr = 0.5 * FsN; % sec (2 vocalizations separated by gap_thr is treated as 2 separate songs)
% 
bpFilt = designfilt('bandpassfir','FilterOrder',20, ...
          'CutoffFrequency1',lowF,'CutoffFrequency2',highF, ...
          'SampleRate',Fs); % defining the filter we designed to detect songs within files
%      
signal_thr_ch1 = 5e-4; % To be adjusted empirically depending upon background noise of recordings

Times = struct('Name',[],'Starts',[],'Stops',[]); % initializing a structure where we record the 
% start and stop of each song detected in a particular file

% this for loop finds songs in each file, records starts and stops to Times
% variable, and writes each song out
for f = 1:nFiles % walk through each file and find songs
    tic % starts recording time
    f_group(idx(f)).name % chosing which file to open
    fid_ch1 = fopen(f_group(idx(f)).name,'r'); % opening file

    [ch1Raw, c1] = fread(fid_ch1,'float32'); % reading in data
    
    fclose('all');
    
    dataOutCh1 = filtfilt(bpFilt,ch1Raw); % zero phase digital filtering on signal

    dsCh1= resample(dataOutCh1.^2,b,a); % downsample but include FIR filter (better than smooth-->decimate fcn)
    
    tMouse = dsCh1>signal_thr_ch1; % any signal higher than noise assigned to tMouse variable
   
    Times(f).Name = f_group(idx(f)).name; % add name of file to Times structure
    
    t = find(tMouse == 1); % where in the file is there signal?
    if (~isempty(t)) % if there is signal somewhere in the file
        
        [temp_starts, temp_stops] =  findStartStops_sks(t,gap_thr); % use findStartsStops function to figure out the start and stop
        % of the song - see other function for details
        song_idx = (temp_stops - temp_starts) > 2 * FsN; % only record the starts and stops of a song if it is longer than 2 seconds
        
        fin_starts = round(temp_starts*(Fs/FsN)); % temp_stops and temp_starts are recorded at downsampled rate. Find the real
        fin_stops = round(temp_stops*(Fs/FsN)); % starts and stops in the original file (regular sampling rate). Note you have to round
        % the value so you have a whole number.
        
        Times(f).Starts = fin_starts(song_idx); % assign starts and stops to the Times structure
        Times(f).Stops = fin_stops(song_idx); 
    else % if there was no signal above noise in the file
        Times(f).Starts = NaN; % record that no songs there
        Times(f).Stops = NaN; % record that no songs there
    end
    if ~isnan(Times(f).Starts) % If there are songs in the file
        for thisSong = 1:length(Times(f).Starts) % create song files for each instance of a song
            thisStart = Times(f).Starts(thisSong); % start time assigned from Times structure
            thisStop = Times(f).Stops(thisSong); % end time assigned from Times structure
            begIndex = round(thisStart - 97656.25); % start file 1 sec before detected 1st note
            endIndex = round(thisStop + 97656.25); % end file 1 sec after last detected note
            if endIndex > 120000000 % if adding 1 sec takes you past end of file, make
                endIndex = 120000000; % end of file the correct length - note my files are 20 min long - need to be changed if hour long
            elseif begIndex < 1 % if song starts within 1 second of original file start
                begIndex = 1; % make song file begin at the beginning of the original file
            end
            chunk = ch1Raw(begIndex:endIndex); % define the song
            filtChunk = bandpass_sks(chunk); % filter song
            figure()
            specgram(filtChunk,512,97656.25); % make a spectrogram
            caxis([-100 20])
        end
    else % if there are no songs in the file
        strcat("No songs detected in ",Times(f).Name) % print out, no songs detected in file-name.
    end

    toc % stop recording time -- will print out how long it took to do the above process.
end