%--------------------------------------------------------------------------
% Author: Morag Brown
% generates input signals for FFT testing 
%--------------------------------------------------------------------------

% create input signal
sig_len = 2048;          % number of points in signal
N = 1024;                 % FFT length
fs = 800e6;              % sampling frequency
t = (1:sig_len)/fs;      % time vector
f = (250)*1e6;           % signal frequency
a = 1/N;                 % signal amplitude
snr1 = 30;               % signal-to-noise ratio

% generate signals
real_sig = a*cos(2*pi*f*t);

an = 10^((20*log10(a/sqrt(2)) - snr1)/10);
noise_sig = sqrt(an)*randn(1,sig_len);

sig = real_sig + noise_sig;

% write signal values to file
f1_name = sprintf('/media/morag/linux_storage/storage_home/university/test_system/fft_testing/simulation/inputs/%d_sine_noise_snr-%ddb_fs-%dmhz_fsig-%dmhz',N,snr1,fs/1e6,f/1e6);
file = fopen(f1_name, 'w');
fprintf(file, '%e\n', sig);
fclose(file);

f2_name = sprintf('/media/morag/linux_storage/storage_home/university/test_system/fft_testing/simulation/inputs/%d_sine_fs-%dmhz_fsig-%dmhz',N,fs/1e6,f/1e6);
file = fopen(f2_name, 'w');
fprintf(file, '%e\n', real_sig);
fclose(file);

f3_name = sprintf('/media/morag/linux_storage/storage_home/university/test_system/fft_testing/simulation/inputs/%d_noise_fs-%dmhz_fsig-%dmhz',N,fs/1e6,f/1e6);
file = fopen(f3_name, 'w');
fprintf(file, '%e\n', noise_sig);
fclose(file);

% % plots
% l = 80;
% 
% figure;
% subplot(3,1,1)
% plot(t(1:l),real_sig(1:l));
% title('real sinusoid')
% subplot(3,1,2)
% plot(t(1:l),noise_sig(1:l));
% title('noise')
% subplot(3,1,3)
% plot(t(1:l),sig(1:l));
% title('real sinusoid + noise')
% xlabel('time (s)')