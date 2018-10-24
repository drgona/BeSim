function sine_waves = time_transform(samples,Ts)

days = samples*Ts/86400;  %24h = 86400 sec
weeks = days/7;  %24h = 86400 sec
months = days/(365/12); 
% amplitude = 1;

x1 = 0:(months*2*pi/samples):months*2*pi-(months*2*pi/samples);
x2 = 0:(weeks*2*pi/samples):weeks*2*pi-(weeks*2*pi/samples);
x3 = 0:(days*2*pi/samples):days*2*pi-(days*2*pi/samples);
time_months = sin(x1)';
time_weeks = sin(x2)';
time_days = sin(x3)';

sine_waves = [time_days time_weeks time_months];

end