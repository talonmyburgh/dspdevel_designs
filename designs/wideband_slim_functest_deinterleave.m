sync_index = find(out.sync_out.data(:)>0,1,'first')
re_out = out.re_out.data(sync_index:end)
im_out = out.im_out.data(sync_index:end)
len = length(re_out)

even = 1:2:len 
odd = 2:2:len 

re_real = re_out(even)
re_imag = im_out(odd)
re_complex = complex(re_real, re_imag)
re_mag = abs(re_complex)
% re_pha = angle(re_complex)

im_real = re_out(odd)
im_imag = im_out(even)
im_complex = complex(im_real, im_imag)
im_mag = abs(im_complex)
% im_pha = angle(im_complex)

subplot(1,2,1)
plot(re_mag)
title('Real Magnitude')
% subplot(2,2,3)
% plot(re_pha)
% title('Real Phase')

subplot(1,2,2)
plot(im_mag)
title('Imag Magnitude')
% subplot(2,2,4)
% plot(im_pha)
% title('Imag Phase')
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(2,2,1)
% plot(re_real)
% title('Real Real')
% subplot(2,2,3)
% plot(re_imag)
% title('Real Imag')
% 
% subplot(2,2,2)
% plot(im_real)
% title('Imag Real')
% subplot(2,2,4)
% plot(im_imag)
% title('Imag Imag')
