sig_len = 1024;          % number of points in signal
N = 1024;                % FFT length
fs = 80e6;              % sampling frequency
t = (1:sig_len)/fs;      % time vector
f = (250)*1e6;           % signal frequency
a = 1/N;                 % signal amplitude
snr1 = 30;               % signal-to-noise ratio
sim_len = 4096;          % how long the simulation must run for
design = 'wideband_functest.slx'; % design we're simulating

cos_sig = a*cos(2*pi*f*t);

an = 10^((20*log10(a/sqrt(2)) - snr1)/10);
zero_arr = zeros(1,sig_len);
imp_train = zero_arr;
slc = 128:128:sig_len;
imp_train(slc) = 0.4;
noise_sig = sqrt(an)*randn(1,sig_len);

%input split up for bram's
d0 = zero_arr(1:2:end);
d1 = imp_train(1:2:end);
d2 = zero_arr(2:2:end);
d3 = imp_train(2:2:end);

%simulate the design
tic;
sim(design, sim_len);
T = toc;
