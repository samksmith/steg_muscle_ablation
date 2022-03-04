function [note_starts, note_ends, note_durs, INI] = msr_note_times(call, samp_rate, threshold, reset, INI_max);
% function [note_starts, note_ends, note_durs, INI] = msr_note_times(call, samp_rate, threshold, reset, INI_max);


% Check input arguments and assign defaults
if nargin<5
    INI_max = 200;           %time in msec that defines transition to new call within recording
end

if (nargin<4)
    reset=1;                %reset is time in msec that need to be subthreshold to be out of note
end

if (nargin<3)
    threshold = 8;          %number of standard deviations above background noise needed to detect call
end

if (nargin<2)
    samp_rate = 195312.5;
end


%Define variables
[r,c] = size(call);
threshold = threshold*std(call(r-1000:r,1));
reset = round(reset*samp_rate/1000); %define criteria (# subthreshold samples) for resetting note counter to being outside a note
[call_dur, c] = size(call); %call_dur determines number of samples to cycle through
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
Hi_crit = 3;                %number of consecutive above threshold samples to indiciate in note
call = abs(call);           %make call amplitude only

spec_resol = 256;           %this should agree with value in msr_note_Hz. will be used to define minimum note length at end.

%Extract call information
for i=1:call_dur
    if (note~=1)
        if(call(i,1)>threshold)   %if currently not in a note, and current input above threshold **** THIS LINE TAKES UP 50% OF COMPUTATION
            
            if ((i-last_hi)==1)             % if not in note, and call(i) above threshold, will ask if last sample above threshold as well
                H=H+1;                      % if last_hi was also above threshold, will increment high counter to H+1
            else
                H=1;                        % if last_hi was not above threshold, will reset H to 1, because need repeated samples above threshold to detect onset of note when note is off
            end
            
            if H>Hi_crit                        %if more than Hi_crit samples are above threshold, will score this as beginning of note
                note_starts(B,1)=i-Hi_crit;     %score this as beginning of note
                B=B+1;                          %advance beginning of note counter
                H=1;                            %reset suprathreshold counter
                note=1;                         %turn on note flag
            end
            last_hi=i;                      % set last_hi to current value of i to store for next pass
        end    
        
    elseif (call(i,1)<threshold)    %if in note but input is beneath threshold **** THIS LINE TAKES UP 40% OF COMPUTATION
                                    %count # consecutive entries below threshold; note considered to end if followed by 10ms samples below thresh
        if ((i-last_lo)==1)                 %if the difference is one, they are consecutive
            L=L+1;                          %increment subthreshold counter
        else
            L=1;                            %This is key -- prevents counter from adding up events around zero-crossings within call
        end
        if L >= reset                   %if exceeds #msecs expected of silence (usually 10)
            note_ends(E,1)=i-reset;
            E=E+1;                      %increments end-of-note counter
            L=1;                        %resets L
            note=0;                     %resets note to 'off'
        end
        last_lo=i;
    end
end


note_num_S = size(note_starts);
note_num_E = size(note_ends);
if note_num_S ~= note_num_E
    error('note numbers are different based on onset and offset!')
end

note_starts = 1000*note_starts/samp_rate; %express times in msec
note_ends = 1000*note_ends/samp_rate;


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
min_note_dur = 5*spec_resol/2 * 1000/samp_rate; % Under default settings, this corresponds to 3.3msec. More importantly, it gives four points on a spectrograph, which is the minimum needed to fit quadratic curve uniquely.
if min(note_durs) <= min_note_dur
    [note_starts, note_ends, note_durs, INI] = remove_subnotes(note_starts, note_ends, min_note_dur);
end