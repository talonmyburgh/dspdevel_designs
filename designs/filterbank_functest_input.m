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

%%%Gather output, calculate theoretical output and compare%%%

%%Generate theoretical window%%
window_coeffs = hanning(filter_length);
filter_coeffs = window_coeffs .* transpose(sinc(fwidth*((1:filter_length)/nof_points - (nof_taps/2))));

%%Interleave recorded inputs%%
interleaved_input = zeros(sim_len, 1);
for i=1:1:wideband_factor
    interleaved_input(i:wideband_factor:end) = simout.data_input.data(1:sim_len/wideband_factor,i);
end
%%Slice interleaved input into filter_lengths%%
input_slices = reshape(interleaved_input, filter_length, []);

%%Interleave recorded output%%
interleaved_output = zeros(size(simout.data_output.data,1)*wideband_factor, 1);
for i=1:1:wideband_factor
    interleaved_output(i:wideband_factor:end) = simout.data_output.data(:,i);
end
%%Determine index at which output-data is valid%%
valid_index = (find(simout.val_output.data(:)>0,1,'first') - 1)*wideband_factor -1; % interleaved wb_factor signals

sim_len_less_valid_index = sim_len-valid_index;
whole_output_last_index = (sim_len - mod(sim_len_less_valid_index, nof_points));
%%Slice interleaved valid-output-data into as many nof_points slices%%
output_slices = reshape(interleaved_output(valid_index:whole_output_last_index-1), nof_points, []);

%%Calculate Theoretical Output%%
slice_count = sim_len/filter_length;
theoretical_output_slices = zeros(nof_points, slice_count); 
%Commutate filter across P rows (resulting in 'nof_taps' cols)
phased_filter = reshape(filter_coeffs, nof_points, nof_taps);
for slice_index=1:slice_count
    %Commutate input-slice across P rows (resulting in 'nof_taps' cols)%%
    phased_input = reshape(input_slices(:,slice_index), nof_points, nof_taps);
    theoretical_output_slices(:, slice_index) = sum(phased_input .* phased_filter, 2); % sum across taps (dim=2)
    %sum(reshape(1:6, 2, 3) .* [1, 1, 1; 2, 2, 2], 2) % proof of the above computation
end

%%Plotting%%
%Select an output slice after all taps have been filled
slice_index = nof_taps;

subplot(1,4,1)
plot(input_slices(:, slice_index))
title('Input')
subplot(1,4,2)
theoretical_section = theoretical_output_slices(:, slice_index)
plot(theoretical_section)
title('Theoretical')
subplot(1,4,3)
output_section = output_slices(:, slice_index+nof_taps-1)
plot(output_section)
title('Output')
subplot(1,4,4)
plot(abs(output_section - theoretical_section))
title('Absolute Difference')