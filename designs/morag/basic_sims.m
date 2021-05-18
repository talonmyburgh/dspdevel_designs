shifting = 1023;
shiftreg = timeseries(shifting);
sim_len = 2048;

%% CASPER
design = 'casper_wb_1024pt_8in.slx';

T = 0;

input_file = '1024_sine_fs-800mhz_fsig-250mhz'
file = fopen(input_file);
re_sig = fscanf(file, '%f');

% parse parameters from file name
h = char(input_file);
ix = strfind(h,'_');  % get underscore locations
N = str2num(h(1:ix(1)-1));
fs = str2num(h(ix(end)-6:ix(end)-4))*1e6;

% demux input signal
d0 = re_sig(1:8:end);
d1 = re_sig(2:8:end);
d2 = re_sig(3:8:end);
d3 = re_sig(4:8:end);
d4 = re_sig(5:8:end);
d5 = re_sig(6:8:end);
d6 = re_sig(7:8:end);
d7 = re_sig(8:8:end);

tic;
% simulate design
sim(design, sim_len);
T = toc;

% interleave output samples and plot
dmux_out = 4;

% find valid data index
val_id = find(sync_out);
val_len = (sim_len + 1)*dmux_out - val_id*dmux_out;

% interleave output
fft_re(1:4:val_len) = out_re(val_id+1:end);
fft_re(2:4:val_len) = out_re1(val_id+1:end);
fft_re(3:4:val_len) = out_re2(val_id+1:end);
fft_re(4:4:val_len) = out_re3(val_id+1:end);

fft_im(1:4:val_len) = out_im(val_id+1:end);
fft_im(2:4:val_len) = out_im1(val_id+1:end);
fft_im(3:4:val_len) = out_im2(val_id+1:end);
fft_im(4:4:val_len) = out_im3(val_id+1:end);

fft_complex = fft_re(1:N/2) + fft_im(1:N/2)*1j;  % discards neg freq components, (N/2 - N) is next spectra
output_fft = abs(fft_complex);

freq_ax = (1:N/2)*(fs/(N));

figure;
subplot(1,1,1);
plot(freq_ax/1e6, output_fft)
xlabel('freq (MHz)')

%% ASTRON simulation
% % still confused about output
% 
% sim_len = 5000;
% 
% % demux input signal
% d0 = real_sig(1:2:end);
% d1 = real_sig(2:2:end);
% 
% % simulate design
% sim('wideband_sine.slx', sim_len)
% 
% % one input only ----------------------------------------------------------
% val_id = find(ans.out_sync1);
% val_len = sim_len + 1 - val_id;
% 
% fft_re_one(1:round(val_len/2)) = ans.re_2(val_id+1:2:end);
% fft_im_one(1:round(val_len/2)) = ans.im_2(val_id+1:2:end);
% 
% fft_complex_one = fft_re_one(1:N/2) + fft_im_one(1:N/2)*1j;
% output_fft_one = abs(fft_complex_one);
% 
% freq_ax_one = (1:(N/2))*(fs/(N));
% 
% figure;
% subplot(1,1,1)
% plot(freq_ax_one/1e6, output_fft_one)
% xlabel('freq (MHz)')
% title('wb\_factor = 1')
% 
% % multi input -------------------------------------------------------------
% 
% % interleave output samples and plot
% dmux_out = 2;
% % 
% val_id = find(ans.out_sync);
% % range is (1:sim_len+1) NOT (0:sim_len)
% val_len = (5000 + 1)*dmux_out - val_id*dmux_out;
% 
% % range is (val_id+1:sim_len+1) NOT (val_id:sim_len)
% 
% % demux two real outputs
% fft_re_1(1:2:round(val_len/2)-1) = ans.re_0(val_id+2:2:end);
% fft_re_1(2:2:round(val_len/2)) = ans.re_1(val_id+2:2:end);
% 
% % fft_re_2(1:2:round(val_len/2)-1) = ans.re_0(val_id+2:2:end);
% % fft_re_2(2:2:round(val_len/2)) = ans.re_1(val_id+2:2:end);
% 
% fft_im_1(1:2:round(val_len/2)-1) = ans.im_0(val_id+2:2:end);
% fft_im_1(2:2:round(val_len/2)) = ans.im_1(val_id+2:2:end);
% 
% % fft_im_2(1:2:round(val_len/2)-1) = ans.im_0(val_id+2:2:end);
% % fft_im_2(2:2:round(val_len/2)) = ans.im_1(val_id+2:2:end);
% 
% % fft_re(1:2:val_len) = ans.re_0(val_id+1:end);
% % fft_re(2:2:val_len) = ans.re_1(val_id+1:end);
% % 
% % fft_im(1:2:val_len) = ans.im_0(val_id+1:end);
% % fft_im(2:2:val_len) = ans.im_1(val_id+1:end);
% 
% 
% fft_complex = fft_re_1(1:N/2) + fft_im_1(1:N/2)*1j;
% output_fft = abs(fft_complex);
% 
% freq_ax = (1:(N/2))*(fs/(N));
% 
% figure;
% subplot(1,1,1)
% plot(freq_ax/1e6, output_fft)
% xlabel('freq (MHz)')
% title('wb\_factor = 2')