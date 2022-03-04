function [note_starts, note_ends, note_durs, INI] = msr_note_times_entropy(song, samp_rate, threshold, reset, INI_max)

% Check input arguments and assign defaults
if nargin<5
    INI_max = 200;           %time in msec that defines transition to new call within recording
end

if (nargin<4)
    reset=10;                %reset is time in msec that need to be above threshold to be out of note
end

if (nargin<3)
    threshold = 4;          %number of standard deviations above background entropy level needed to detect call
end

if (nargin<2)
    samp_rate = 195312.5;
end

%get entropy measure to work off of--note, this takes quite a bit of the
%processing time
[entropy] = read_entropy(song, samp_rate);

%clean out the inevitable NaNs
entropy = entropy(isfinite(entropy));

%get a baseline for the top of the entropy reading
peak = max(entropy);

variable = mean(entropy);
nastd = nanstd(entropy(1:100));

%set some basic values
threshold =  variable - (threshold*nastd);
%reset = round(reset*samp_rate/1000); %define criteria (# subthreshold samples) for resetting note counter to being outside a note
[song_dur, c] = size(entropy); %get the length of the song defined
B=1;                        %counts beginning of notes
H=1;                        %suprathreshold counter
E=1;                        %end note counter
L=1;                        %subthreshold counter
note = 0;                   %logical flag to keep track of whether in note or out of note
last_lo=1;                  %running index of last subthreshold sample
last_hi=1;                  %running index of last suprathreshold sample
note_starts=1;
note_ends=1;
INI(1,1)=1;
Lo_crit = 3;                %number of consecutive above threshold samples to indiciate in note

spec_resol = 256;           %this should agree with value in msr_note_Hz. will be used to define minimum note length at end.

%Extract call information
for i=1:song_dur
    if (note~=1)
        if(entropy(i,1)<threshold)   %if currently not in a note, and current input below threshold **** THIS LINE TAKES UP 50% OF COMPUTATION
            
            if ((i-last_lo)==1)             % if not in note, and entropy(i) below threshold, will ask if last sample below threshold as well
                L=L+1;                      % if last_lo was also below threshold, will increment lows counter to L+1
            else
                L=1;                        % if last_lo was not below threshold, will reset L to 1, because need repeated samples above threshold to detect onset of note when note is off
            end
            
            if L>Lo_crit                        %if more than Lo_crit samples are below threshold, will score this as beginning of note
                note_starts(B,1)=i-Lo_crit;     %score this as beginning of note
                B=B+1;                          %advance beginning of note counter
                H=1;                            %reset suprathreshold counter
                note=1;                         %turn on note flag
            end
            last_lo=i;                      % set last_hi to current value of i to store for next pass
        end    
        
    elseif (entropy(i,1)>threshold)    %if in note but input is above threshold **** THIS LINE TAKES UP 40% OF COMPUTATION
                                       %count # consecutive entries below threshold; note considered to end if followed by 10ms samples below thresh
        if ((i-last_hi)==1)                 %if the difference is one, they are consecutive
            H=H+1;                          %increment subthreshold counter
        else
            H=1;                            %This is key -- prevents counter from adding up events around zero-crossings within call
        end
        if H >= reset                   %if exceeds #msecs expected of silence (usually 10)
            note_ends(E,1)=i-reset;
            E=E+1;                      %increments end-of-note counter
            L=1;                        %resets L
            note=0;                     %resets note to 'off'
        end
        last_hi=i;
    end
end


note_num_S = size(note_starts);
note_num_E = size(note_ends);
if note_num_S ~= note_num_E
    error('note numbers are different based on onset and offset!')
end

note_starts = 256000*note_starts/samp_rate; %express times in msec
note_ends = 256000*note_ends/samp_rate;


%note_ends = note_ends - note_starts(1,1);
%note_starts = note_starts - note_starts(1,1);   %set beginning of call as time zero

note_durs = note_ends - note_starts;

INI = note_starts(2:note_num_S) - note_ends(1:note_num_E-1);

% Make sure there aren't multiple calls being included
if max(INI) > INI_max 
    [note_starts, note_ends, note_durs, INI] = find_long_call(note_starts, note_ends, INI, INI_max);
end

if max(INI) == INI_max
    warning('Warning: The inter-note interval reaches maximum allowed -- check to see if end notes lost from call. May need to use a higher value for "INI_max".')
end

if min(INI) == reset
    warning('Warning: The inter-note interval reaches the minimum allowed -- check spectrogram to see if gaps within notes are interpreted as INIs. May need to use a higher value for "reset".')
end

% Make sure no sub-notes or extraneous noised are being included
%min_note_dur = 5*spec_resol/2 * 1000/samp_rate; % Under default settings, this corresponds to 3.3msec. More importantly, it gives four points on a spectrograph, which is the minimum needed to fit quadratic curve uniquely.
%if min(note_durs) <= min_note_dur
%    [note_starts, note_ends, note_durs, INI] = remove_subnotes(note_starts, note_ends, min_note_dur);
%end

