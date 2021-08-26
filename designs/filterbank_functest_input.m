sig_len = 64;           % number of points in signal
fs = 80e6;               % sampling frequency
t = (1:sig_len)/fs;      % time vector
f = (250)*1e6;           % signal frequency
a = 1/sig_len;           % signal amplitude
snr1 = 30;               % signal-to-noise ratio
sim_len = 4096;          % how long the simulation must run for
design = 'filterbank_functest.slx'; % design we're simulating

an = 10^((20*log10(a/sqrt(2)) - snr1)/10);
noise_sig0 = zeros(1,sig_len);
noise_sig0(1:1:sig_len) = sqrt(an)*randn(1,sig_len);
square_sig1 = zeros(1,sig_len);
square_sig1(sig_len/4:1:3*sig_len/4) = 1.0/1024.0;

%input split up for bram's
d0 = noise_sig0(1:sig_len);
d1 = square_sig1(1:sig_len);

%simulate the design
tic;
simout = sim(design, sim_len);
T = toc;

val_index = find(simout.val_output.data(:)>0,1,'first')
