sig_len = 1024;          % number of points in signal
N = 1024;                % FFT length
fs = 80e6;              % sampling frequency
t = (1:sig_len)/fs;      % time vector
f = (250)*1e6;           % signal frequency
a = 1/N;                 % signal amplitude
snr1 = 30;               % signal-to-noise ratio
sim_len = 4096;          % how long the simulation must run for
design = 'wideband_slim_functest.slx'; % design we're simulating

real_sig = a*cos(2*pi*f*t);
an = 10^((20*log10(a/sqrt(2)) - snr1)/10);
noise_sig = sqrt(an)*randn(1,sig_len);

%input split up for bram's
d0 = noise_sig(1:2:end);
d1 = real_sig(1:2:end);
d2 = noise_sig(2:2:end);
d3 = real_sig(2:2:end);

%simulate the design
tic;
sim(design, sim_len);
T = toc;
