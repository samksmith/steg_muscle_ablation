function [note_starts, note_ends, note_starts_T, note_ends_T] = msr_call_amps(call, samp_rate, threshold, reset); 
% function [note_starts, note_ends, note_starts_T, note_ends_T] = msr_call_amps(call, samp_rate,threshold, reset); 


% Check input arguments and assign defaults 
if (nargin<4) 
   reset=5;                %reset is time in msec that need to be subthreshold to be out of note 
end 

if (nargin<3) 
   threshold = 3*std(call(1:1000)); 
end 

if (nargin<2) 
   samp_rate = 195312.5; 
end 


%Define variables 
reset = round(reset*samp_rate/1000); 
[call_dur, c] = size(call); %call_dur determines number of samples to cycle through 
B=1;                        %counts beginning of notes 
H=1;                        %suprathreshold counter 
E=1;                        %end note counter 
L=1;                        %subthreshold counter 
note = 0;                   %logical flag to keep track of whether in note or out of note 
last2lo=[1,1];              %running index of last two subthreshold samples 
last2hi=[1,1];              %running index of last two suprathreshold samples 
note_starts=1; 
note_ends=1; 
Hi_crit = 3;                %number of consecutive above threshold samples to indiciate in note 
call = abs(call);           %make call amplitude only 


%Extract call information 
for i=1:call_dur 
   if (note~=1)&&(call(i)>threshold)   %if currently not in a note, and current input above threshold 

       last2hi(1,1)=last2hi(1,2); 
       last2hi(1,2) = i; 
       if ((last2hi(1,2)-last2hi(1,1))==1) 
           H=H+1; 
       else 
           H=1; 
       end 

       if H>Hi_crit 
           note_starts(B,1)=i-Hi_crit;     %score this as beginning of note 
           B=B+1;                          %advance note counter 
           H=1;                            %reset suprathreshold counter 
           note=1;                         %turn on note flag 
       end 
   end 

   if ((note==1)&&(call(i)<threshold))     %if in note but input is beneath threshold 
                                           %count # consecutive entries below threshold 
                                           %note considered to end if followed by 5ms samples below thresh 
       last2lo(1,1)=last2lo(1,2);          %assign last subthreshold index to LAST2(1,1) 
       last2lo(1,2)=i;                     %assign current index to LAST2(1,2) 
       if ((last2lo(1,2)-last2lo(1,1))==1) %if the difference is one, they are consecutive 
           L=L+1;                          %increment subthreshold counter 
       else 
           L=1;                            %This is key -- prevents counter from adding up events around zero-crossings within call 
       end 
       if L == reset                   %if exceeds #msecs expected of silence (usually 5) 
           note_ends(E,1)=i-reset; 
           E=E+1;                      %increments end-of-note counter 
           L=1;                        %resets L 
           note=0;                     %resets note to 'off' 
       end 
   end 
end 

note_num_S = size(note_starts); 
note_num_E = size(note_ends); 
if note_num_S ~= note_num_E 
   error('note numbers are different based on onset and offset!') 
end 

note_starts_T = 1000*note_starts/samp_rate; %express times in msec 
note_ends_T = 1000*note_ends/samp_rate; 


%note_ends = note_ends - note_starts(1,1); 
%note_starts = note_starts - note_starts(1,1);   %set beginning of call as time zero 

%note_durs = note_ends - note_starts; 

%for i = 1:(note_num_S - 1) 
%   INI(i)=note_starts(i+1) - note_ends(i); 
%end 
