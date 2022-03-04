function [note_max_Hz, note_min_Hz, Hz_curve, note_resid] = msr_note_Hz(note, sampling_rate, spec_resol, method, fig_on);
% function [note_max_Hz, note_min_Hz, Hz_curve, note_resid] =
% msr_note_Hz(note, sampling_rate, spec_resol, method, fig_on);

%check inputs and assign defaults
if nargin<5
    fig_on = 0;
end
if nargin<4
    method = 'q';
end
if nargin<3
    spec_resol = 256;
end
if nargin<2
    sampling_rate = 195312.5;
end

[note_dur, c] = size(note); %check whether note is a column or row vector
if note_dur <= 1
    error('note too short')
end

%---------

[spec, F, T] = specgram(note, spec_resol, sampling_rate); %compute spectrogram
[power_pks, pk_index] = max(spec);

[end_time, c] = size(T); % in call to note hz, expand definition of note to include +/- NFFT/2 on either side of note boundaries.  
for i=1:end_time
    pk_Hz(i,1) = F(pk_index(i));
end

note_max_Hz = max(pk_Hz);
note_min_Hz = min(pk_Hz);


switch(method) %haven't looked at methods 'e' and 'e2'
case('q')
    [Hz_curve,s] = polyfit(T, pk_Hz, 2); % fits pk frequency data to second order (quadratic) polynomial, where p is the model and s is used by polyval to compute size of standard dev in residuals
    [y, note_SD] = polyval(Hz_curve, T, s);
    note_resid = mean(note_SD);
    if(fig_on)
        plot(T, y, '-', T, pk_Hz, 'o');
    end
    
case('e')
    [Hz_curve,s] = polyfit(T, log10(pk_Hz), 1); % fits log of pk frequency data to straight line
    [y, note_SD] = polyval(T, Hz_curve, s);
    if(fig_on)
        plot(T, y, '-', T, log10(pk_Hz), 'o');
    end
    
case('e2')
    [Hz_curve,s] = polyfit(T, log10(pk_Hz), 2); % fits log of pk frequency data to quadratic
    [y, note_SD] = polyval(T, Hz_curve, s);
    if(fig_on)
        plot(T, y, '-', T, log10(pk_Hz), 'o');
    end
end