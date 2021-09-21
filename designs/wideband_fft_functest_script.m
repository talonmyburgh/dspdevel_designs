%FFT parameters
nof_points = 1024;                % FFT length
wideband_factor = 2;
pipeline_nof_points = nof_points/wideband_factor;
dual_processing = 1;
reorder_freq = 0;
shift_schedule = 0;
input_data_width = 18;
stage_data_width = 18;
output_data_width = 18;

%Signal generation parameters
sig_len = 1024;          % number of points in signal
fs = 0.1e9;               % sampling frequency
t = (1:sig_len)/fs;      % time vector
f = 32*1e6;           % signal frequency
a = (2^input_data_width)/8; % signal amplitude
snr1 = 30;               % signal-to-noise ratio
bram_depth = 512;
bram_addr = log2(bram_depth);

k=1:nof_points;
sig_1 = a*sin((441*2*pi/nof_points)*k); % expect delta in the second bin
sig_2 = a*sin((521*2*pi/nof_points)*k); % expect delta in the fourth bin
an = 10^((20*log10(a/sqrt(2)) - snr1)/10);
delta_sig = zeros(sig_len,1);
delta_sig(7) = 750;
real_sig = delta_sig;

% zero_sig = zeros(sig_len)

% noise_sig = sqrt(an)*randn(1,sig_len);
noise_sig = delta_sig;

%Simulation parame7ters
sim_len = 8192;                   % how long the simulation must run for
design = 'wideband_functest.slx'; % design we're simulating

%BRAM configuration
d0 = real_sig(1:wideband_factor:end);
d1 = noise_sig(1:wideband_factor:end);
d2 = real_sig(2:wideband_factor:end);
d3 = noise_sig(2:wideband_factor:end);
d4 = real_sig(3:wideband_factor:end);
d5 = noise_sig(3:wideband_factor:end);
d6 = real_sig(4:wideband_factor:end);
d7 = noise_sig(4:wideband_factor:end);

%Simulate
tic;
simout=sim(design, sim_len);
T = toc;

%Collect Inputs from the design
in_re  = reshape(simout.re_input.data(1:sim_len,:)', [], 1);
in_im  = reshape(simout.im_input.data(1:sim_len,:)', [], 1);

%Collect Outputs from the design and process them in wideband_fft_process_output()
dv_index = find(simout.dv_out.data(:)>0,1,'first');
out_re = simout.re_out.data(dv_index:end,:);
out_im = simout.im_out.data(dv_index:end,:);

[output_a, output_b, output_x] = wideband_fft_process_output(out_re, out_im, wideband_factor, nof_points, dual_processing);
[theoretical_output_a, theoretical_output_b] = wideband_fft_model(in_re, in_im, nof_points, dual_processing, reorder_freq);

% plot(x,abs(output_b(:,end)),x,abs(theoretical_output_b(:,end)))