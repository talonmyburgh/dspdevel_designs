% PFB Mask Values
wideband_factor = 2;     % parallel ingress/egress rate
nof_streams   = 1;       % number of streams
nof_points = 1024;       % FFT size but also number of phases in PFB
nof_taps   = 2;          % signal_length/nof_points
fwidth = 0.8;
pipeline_point_size = round(wideband_factor/nof_points);
dual_processing = 1;
reorder_freq = 0;
shift_schedule = 1023;
input_data_width = 18;
prefilter_out_data_width = 18;
fft_out_data_width = 18;
enable_prefilter = true;

% Simulation values
design = 'pfb_wideband_functest.slx'; % design we're simulating 
sim_len = 16000;          % how long the simulation must run for

%Signal generation parameters
sig_len = 1024;          % number of points in signal
fs = 0.1e9;               % sampling frequency
t = (1:sig_len)/fs;      % time vector
f = 32*1e6;           % signal frequency
a = (2^input_data_width)/10; % signal amplitude
snr1 = 30;               % signal-to-noise ratio
bram_depth = 512;
bram_addr = log2(bram_depth);

k=1:nof_points;
sig_1 = a*sin((32*2*pi/nof_points)*k); % expect delta in the second bin
sig_1 = sig_1;
sig_2 = zeros(sig_len);

%%fill BRAMs%%
d0 = sig_2(1:wideband_factor*nof_streams:nof_points);
d1 = sig_1(1:wideband_factor*nof_streams:nof_points);
d2 = sig_2(2:wideband_factor*nof_streams:nof_points);
d3 = sig_1(2:wideband_factor*nof_streams:nof_points);

%%%Simulate the design%%%
tic;
simout = sim(design, sim_len);
T = toc;

%Inputs
% %take 1 real spectra out - neg and pos frequencies if ~dual, else A_re and
% %B_re pos frequencies
numspectra = 2;
len = numspectra*nof_points;

%%Interleave recorded inputs - real and imaginary%%
%Collect Inputs from the design
in_re  = reshape(simout.re_input.data(1:sim_len,:)', [], 1);
in_im  = reshape(simout.im_input.data(1:sim_len,:)', [], 1);

%%Determine index at which filter-data is valid%%
fil_valid_index = (find(simout.dv_fil.data(:)>0,1,'first')); % interleaved wb_factor signals

%%Determine index at which output-data is valid%%
out_valid_index = (find(simout.dv_out.data(:)>0,1,'first')); % interleaved wb_factor signals

fil_output_array_re = simout.re_fil.data(fil_valid_index:end,:);
fil_output_array_im = simout.im_fil.data(fil_valid_index:end,:);
pfb_output_re = simout.re_out.data(out_valid_index:end,:);
pfb_output_im = simout.im_out.data(out_valid_index:end,:);

prefilter_output_re = prefilter_process_output(fil_output_array_re, wideband_factor, nof_points);
prefilter_output_im = prefilter_process_output(fil_output_array_im, wideband_factor, nof_points);
[output_a, output_b, output_x] = wideband_fft_process_output(pfb_output_re, pfb_output_im, wideband_factor, nof_points, dual_processing);
[theoretical_prefilter_output,theoretical_pfb_output_a,theoretical_pfb_output_b,theoretical_pfb_output_c] = pfb_model(in_re, ...
    in_im, nof_points, dual_processing, reorder_freq, nof_taps, fwidth,enable_prefilter);
%%Plotting%%
%Select an output slice to display
slice = 10;
scale = 1024;
subplot(6,1,1);
plot(in_re(1:nof_points));
title('Real Input');
subplot(6,1,2);
plot(in_im(1:nof_points));
title('Imaginary Input');
subplot(6,1,3);
x_prefilter = 0:1:nof_points-1;
plot(x_prefilter,prefilter_output_re(:,slice),x_prefilter,real(theoretical_prefilter_output(:,slice)));
title('Real prefilter output');
legend('HDL', 'Ideal')
subplot(6,1,4);
plot(x_prefilter,prefilter_output_re(:,slice),x_prefilter,real(theoretical_prefilter_output(:,slice)));
title('Imaginary prefilter output');
legend('HDL', 'Ideal')
subplot(6,1,5);
x_pfb = 0:1:nof_points/2 -1;
plot(x_pfb,abs(output_a(:,slice)),x_pfb,abs(theoretical_pfb_output_a(:,slice))/scale);
title("PFB output A, actual and theoretical");
legend('HDL', 'Ideal')
subplot(6,1,6);
plot(x_pfb,abs(output_b(:,slice)),x_pfb,abs(theoretical_pfb_output_b(:,slice))/scale);
title("PFB output B, actual and theoretical");
legend('HDL', 'Ideal')


