% Filterbank Mask values
wideband_factor = 2;     % parallel ingress/egress rate
nof_points = 128;        % FFT size but also number of phases in PFB
nof_taps   = 4;          % signal_length/nof_points
fwidth = 0.8;

% Simulation values
sim_len = 4096;          % how long the simulation must run for
design = 'filterbank_functest.slx'; % design we're simulating

% Inferred values
filter_length = nof_points*nof_taps;

%%%Setup input_signal and populate simulink BRAMs%%%

%%noise input_sig%%
a = 1.0;                 % signal amplitude
snr1 = 30;               % signal-to-noise ratio
an = 10^((20*log10(a/sqrt(2)) - snr1)/10);
input_sig(1:1:nof_points) = sqrt(an)*randn(1,nof_points);
%%square input_sig%%
% input_sig= zeros(nof_points,1);
% input_sig(nof_points/4:1:3*nof_points/4) = 1.0/1024.0;
%%pulse input_sig%%
% input_sig = zeros(nof_points,1);
% input_sig(32) = 1;
%%ramp input_sig%%
% input_sig(:, 1) = (1:nof_points)/(nof_points*2);

%%fill BRAMs%%
d0 = input_sig(1:wideband_factor:nof_points);
d1 = input_sig(2:wideband_factor:nof_points);

%%%Simulate the design%%%
tic;
simout = sim(design, sim_len);
T = toc;

%%Interleave recorded inputs%%
interleaved_input = zeros(sim_len*wideband_factor, 1);
for i=1:1:wideband_factor
    interleaved_input(i:wideband_factor:end) = simout.data_input.data(1:sim_len,i);
end

%%Slice interleaved input into filter_lengths%%
input_slice_count = ((size(interleaved_input, 1)-filter_length)/nof_points);
input_slices = zeros(filter_length, input_slice_count);
for slice_index=0:input_slice_count
    input_slices(:, slice_index+1) = interleaved_input(slice_index*nof_points+1:slice_index*nof_points+filter_length, 1);
end

%%Determine index at which output-data is valid%%
dv_index = find(simout.val_output.data(:)>0,1,'first')-1;
%%Interleave recorded output%%
output_array = simout.data_output.data(dv_index:end,:);
prefilter_output = prefilter_process_output(output_array, wideband_factor, nof_points);

%%Calculate Theoretical Output%%
theoretical_prefilter_output = prefilter_model(interleaved_input, nof_points, nof_taps, fwidth);

%%Plotting%%
%Select an output slice after all taps have been filled
slice_index = max(nof_taps-1, 1);
slice_count = 4;

slice_end_index = slice_index+slice_count-1;
filter_output_x = slice_index*nof_points:(slice_end_index+1)*nof_points-1;

subplot(4,1,1)
% plot(reshape(input_slices(:, slice_index:slice_end_index), [], 1))
plot(input_slices(:, slice_index:slice_end_index))
title('Input')
subplot(4,1,2)
theoretical_section = reshape(theoretical_prefilter_output(:, slice_index:slice_end_index), [], 1)
plot(filter_output_x,theoretical_section)
title('Theoretical')
subplot(4,1,3)
output_section = reshape(prefilter_output(:, slice_index:slice_end_index), [], 1)
plot(filter_output_x,output_section)
title('Output')
subplot(4,1,4)
plot(filter_output_x,abs(output_section - theoretical_section))
title('Absolute Difference')