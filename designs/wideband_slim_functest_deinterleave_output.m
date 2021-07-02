%FFT deinterleaving side of things
N = 1024;                % FFT length
P = 2;
M = N/P;
dual = 1;
reorder_freq = 0;

%Inputs
% %take 1 real spectra out - neg and pos frequencies if ~dual, else A_re and
% %B_re pos frequencies
numspectra = 2;
len = numspectra*N;
dv_index = find(ans.dv_out.data(:)>0,1,'first');
re_in  = reshape(ans.re_input.data(1:1024,:)',[2048,1]);
im_in  = reshape(ans.im_input.data(1:1024,:)',[2048,1]);

%outputs
re_out = ans.re_out.data(dv_index:dv_index+len,:)'; 
im_out = ans.im_out.data(dv_index:dv_index+len,:)';

%complex out index
if ~dual
    Xre = zeros(P,M);
    Xim = zeros(P,M);
    f_tmp = 0:P:(M-1)*P;
    f = zeros(P,M);
    for p=1:1:P
        f(p,:) = p + f_tmp;
        Xre(p,:) = re_out(f(p,:));
        Xim(p,:) = im_out(f(p,:));
    end
    X=Xre + 1j*Xim;

    ideal = fft(re_in + 1j * im_in);
    half_fft = N/2;
    if reorder_freq                                       %frequencies are needed to be [neg, 0, pos]
        ideal = fftshift(ideal);
    else                                                  %frequencies are [0, pos, neg]
    end

    %Scale the MATLAB FFT to be the same size as the ASTRON one
    ideal = normaliseandraise(X, ideal);
    
    %Plotting
    x = 1:1:N;
    subplot(3,2,1)
    plot(re_in)
    title('Real Input')
    subplot(3,2,2)
    plot(im_in)
    title('Imag Input')
    subplot(3,2,3)
    plot(x,Xre,x,real(ideal))
    title('X - real')
    legend('Astron','Matlab');
    subplot(3,2,4)
    plot(x,Xim,x,imag(ideal))
    title('X - imag')
    legend('Astron','Matlab');
    subplot(3,2,5)
    plot(x,abs(X),x,abs(ideal));
    title('X - Magnitude');
    legend('Astron','Matlab');

%2 real out index
else
    even = 2:2:M;
    odd = 1:2:M;
    Are = zeros(P,M/2);
    Aim = zeros(P,M/2);
    Bre = zeros(P,M/2);
    Bim = zeros(P,M/2);
    
    for p=1:1:P
        Are(p,:) = re_out(p,even);
        Bre(p,:) = re_out(p,odd);
        Aim(p,:) = im_out(p,even);
        Bim(p,:) = im_out(p,odd);
    end
    Are = reshape(Are,[1,N/2]);
    Aim = reshape(Aim,[1,N/2]);
    Bre = reshape(Bre,[1,N/2]);
    Bim = reshape(Bim,[1,N/2]);
    
    A = Are+1j*Aim;
    B = Bre+1j*Bim;

    %ideal fft... we need to split the output to get what we see from HDL.
    ideal = fft(re_in + 1j * im_in);
    R_ideal = real(ideal);
    R_flip = R_ideal(:);
    R_flip(2:end) = R_flip(end:-1:2);
    
    I_ideal = imag(ideal);
    I_flip = I_ideal(:);
    I_flip(2:end) = I_flip(end:-1:2);
    
    B_ideal = (1/2) .* (R_ideal + 1j * I_ideal + R_flip - 1j * I_flip);
    B_ideal = B_ideal(1:N/2);
    A_ideal = (1/(2*1j)) * (R_ideal + 1j * I_ideal - R_flip + 1j * I_flip);
    A_ideal = A_ideal(1:N/2);

    % %Scale the MATLAB FFT to be the same size as the ASTRON one
    A_ideal = normaliseandraise(A,A_ideal);
    B_ideal = normaliseandraise(B,B_ideal);

    %Plotting
    x = 1:1:N/2;
    subplot(4,2,1)
    plot(im_in)
    title('A Input')
    subplot(4,2,2)
    plot(re_in)
    title('B Input')
    subplot(4,2,3)
    plot(x,Are,x,real(A_ideal))
    title('A - real')
    legend('Astron','Matlab');
    subplot(4,2,4)
    plot(x,Aim,x,imag(A_ideal))
    title('A - imag')
    legend('Astron','Matlab');
    subplot(4,2,5)
    plot(x,Bre,x,real(B_ideal))
    title('B - real')
    legend('Astron','Matlab');
    subplot(4,2,6)
    plot(x,Bim,x,imag(B_ideal))
    title('B - imag')
    legend('Astron','Matlab');
    subplot(4,2,7)
    plot(x,abs(A),x,abs(A_ideal))
    title('A - magnitude')
    legend('Astron','Matlab');
    subplot(4,2,8)
    plot(x,abs(B),x,abs(B_ideal))
    title('B - magnitude')
    legend('Astron','Matlab');
end

%assumes complex array inputs
function aff_arr = normaliseandraise(arr_to_match, arr_to_alter)
    arr_to_match_real = real(arr_to_match);
    arr_to_match_imag = imag(arr_to_match);
    arr_to_alter_real = real(arr_to_alter);
    arr_to_alter_imag = imag(arr_to_alter);

    %DC of both:
    arr_to_match_re_dc = mean(abs(arr_to_match_real));
    arr_to_match_im_dc = mean(abs(arr_to_match_imag));
    arr_to_alter_re_dc = mean(abs(arr_to_alter_real));
    arr_to_alter_im_dc = mean(abs(arr_to_alter_imag));

    %zero DC:
    arr_to_match_real_0dc = arr_to_match_real-arr_to_match_re_dc;
    arr_to_match_imag_0dc = arr_to_match_imag-arr_to_match_im_dc;
    arr_to_alter_real = arr_to_alter_real-arr_to_alter_re_dc;
    arr_to_alter_imag = arr_to_alter_imag-arr_to_alter_im_dc;
    
    %scale:
    arr_to_match_re_scale = max(abs(arr_to_match_real_0dc));
    arr_to_match_im_scale = max(abs(arr_to_match_imag_0dc));
    arr_to_alter_re_scale = max(abs(arr_to_alter_real));
    arr_to_alter_im_scale = max(abs(arr_to_alter_imag));

    if(arr_to_alter_re_scale == 0)
        arr_to_alter_re_scale = 1;
    end
    if(arr_to_alter_im_scale == 0)
        arr_to_alter_im_scale = 1;
    end
    
    arr_to_alter_real = arr_to_match_re_dc+(arr_to_alter_real/arr_to_alter_re_scale)*arr_to_match_re_scale;
    arr_to_alter_imag = arr_to_match_im_dc+(arr_to_alter_imag/arr_to_alter_im_scale)*arr_to_match_im_scale;

    aff_arr = arr_to_alter_real + arr_to_alter_imag*1i;
end