%FFT parameters
nof_points = 1024;                % FFT length
wideband_factor = 2;
pipeline_nof_points = nof_points/wideband_factor;
dual_processing = 1;
reorder_freq = 0;
shift_schedule = nof_points-1;
input_data_width = 18;
stage_data_width = 18;
output_data_width = 18;

%Signal generation parameters
sig_len = 1024;          % number of points in signal
fs = 0.1e9;               % sampling frequency
t = (1:sig_len)/fs;      % time vector
f = 32*1e6;           % signal frequency
a = (2^input_data_width)/20; % signal amplitude
snr1 = 30;               % signal-to-noise ratio
bram_depth = 512;
bram_addr = log2(bram_depth);

k=1:nof_points;
sig_1 = a*sin((11*2*pi/nof_points)*k); % expect delta in the second bin
sig_2 = a*sin((521*2*pi/nof_points)*k); % expect delta in the fourth bin
an = 10^((20*log10(a/sqrt(2)) - snr1)/10);
delta_sig = zeros(sig_len,1);
delta_sig(7) = 750;
real_sig = sig_1;

% zero_sig = zeros(sig_len)

noise_sig = sqrt(an)*randn(1,sig_len);
% noise_sig = delta_sig;

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
[theoretical_output_a, theoretical_output_b, theoretical_output_x] = wideband_fft_model(in_re, in_im, nof_points, dual_processing, reorder_freq);

if dual_processing
	x  = 1:size(output_a, 1);
	out_a = abs(output_a(:, end))
	out_b = abs(output_b(:, end))
	theory_a = abs(theoretical_output_a(:, end))
	theory_b = abs(theoretical_output_b(:, end))

	% normalise the output signals
	out_a = out_a/max(out_a)
	out_b = out_b/max(out_b)
	theory_a = theory_a/max(theory_a)
	theory_b = theory_b/max(theory_b)

	subplot(3,2,1)
	plot(in_re(1:length(x)))
	title('Input A (Real)')
	subplot(3,2,2)
	plot(in_im(1:length(x)))
	title('Input B (Real)')
	
	subplot(3,2,3)
	plot(x,out_a,x,theory_a)
	title('Output A (Abs)')
	legend('HDL','Ideal');
	subplot(3,2,4)
	plot(x,out_b,x,theory_b)
	title('Output B (Abs)')
	legend('HDL','Ideal');
	
	subplot(3,2,5)
	plot(abs(theory_a - out_a))
	title('Output A Difference')
	subplot(3,2,6)
	plot(abs(theory_b - out_b))
	title('Output B Difference')	 
else
	
	x  = 1:size(output_x, 1);
	out_x_re = real(output_x(:, end))
	out_x_im = imag(output_x(:, end))
	theory_re = real(theoretical_output_x(:, end))
	theory_im = imag(theoretical_output_x(:, end))

	% normalise the output signals
	out_re = out_re/max(out_re)
	out_im = out_im/max(out_im)
	theory_re = theory_re/max(theory_re)
	theory_im = theory_im/max(theory_im)

	subplot(3,2,1)
	plot(in_re)
	title('Real Input')
	subplot(3,2,2)
	plot(in_im)
	title('Imag Input')
	
	subplot(3,2,3)
	plot(x,out_re,x,theory_re)
	title('Real Output')
	legend('HDL','Ideal');
	subplot(3,2,4)
	plot(x,out_im,x,theory_im)
	title('Imag Output')
	legend('HDL','Ideal');
	
	subplot(3,2,5)
	plot(abs(theory_re - out_re))
	title('Real Output Difference')
	subplot(3,2,6)
	plot(abs(theory_im - out_im))
	title('Imag Output Difference')	 
end
