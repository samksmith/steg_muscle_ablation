function [data_curve, data_resid] = msr_curve(data, note_starts, note_ends, time, method, fig_on, var_name, call_id)

% [data_curve, data_resid] = msr_curve(data, note_starts, note_ends, time, method, fig_on, var_name, call_id)
%
% This function takes an input of data, and several optional variables, and fits a quadratic curve
% to the changes in note structure over time. The NOTE_STARTS and NOTE_ENDS variables are used only if the
% time axis for the curve is defined in msec from onset, rather than in note number. The TIME axis is
% defined as either 'note' or 'msec'. The default is 'note'. The only method allowed for curve fitting
% currently is a quadratic equation, defined by the input 'q', and also the default. I may add others if
% they seem useful. The FIG_ON variable is a logic switch. When it's true (=1), the routine will plot the
% curve and data in a currently active window. If a variable name is given, this will show up as the
% plot's title. If this routine is called by MSR_ALL, the plot is embedded within the active figure.
%
% S. Phelps, 6-3-2007

if nargin<8, call_id = 'call curves plot'; end
if nargin<7, var_name = ''; end
if nargin<6, fig_on=0; end
if nargin<5, method = 'q'; end
if nargin<4, time = 'note'; end
if nargin<3
    note_starts =0;
    note_ends = 0;
end

[note_num, c] = size(data);

switch(time)
case('note')
    T = 0:1:(note_num-1);                     % This will define "time" as note number, with first note as note "zero".
    T = T';                                   % Defining first note as zero means intercepts are estimates of first note.
case('msec')
    T = (note_starts+note_ends)/2;        % Defines time as center of note in msec, from onset of call.
end


switch(method)                                % may want to subtract min value from data, add it back after curve fitting
case('q')
    
    [data_curve,s] = polyfit(T, data, 2); % Fits data to second order (quadratic) polynomial, where s is used by polyval to compute size of residuals
    [y, data_resid] = polyval(data_curve, T, s);
    data_resid = mean(data_resid);
    
    if(fig_on)
        plot(T, y, '-', T, data, 'o'), xlabel(['time (', time, ')']);
        if nargin == 7
            title(var_name)
        end
    end
end

