N = 1024;
P = 1;
M = N/P;
dual = 1;
reorder_freq = 0;

%Inputs
% %take 1 real spectra out - neg and pos frequencies if ~dual, else A_re and
% %B_re pos frequencies
numspectra = 1;
len = N * numspectra;
dv_index = find(out.dv_out.data(:)>0,1,'first');
re_in  = out.re_input.data(1:len);
im_in  = out.im_input.data(1:len);

% tvals = 1:1:N;
% in_im_vec = zeros(N,2);
% in_im_vec(:,1) = tvals;
% in_im_vec(17,2) = 0.3; 
% 
% in_re_vec = zeros(N,2);
% in_re_vec(:,1) = tvals;


%input index
t_tmp = 0:P:(M-1)*P;
t = zeros(P,M);
for p=1:1:P
    t(p,:) = p + t_tmp;
end
% in_im_wrk = zeros(N,2);
% in_re_wrk = zeros(N,2);
% in_im_wrk(:,2) = in_im_vec(t(1,:),2);
% in_im_wrk(:,1) = in_im_vec(:,1);
% in_re_wrk(:,1) = in_re_vec(:,1);
% in_re_wkr(:,2) = in_re_vec(t(1,:),2);

%outputs
re_out = out.re_out.data(dv_index:dv_index+len,:)'; 
im_out = out.im_out.data(dv_index:dv_index+len,:)';

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
%real out index
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
end

if dual
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
    subplot(4,2,4)
    plot(x,Aim,x,imag(A_ideal))
    title('A - imag')
    subplot(4,2,5)
    plot(x,Bre,x,real(B_ideal))
    title('B - real')
    subplot(4,2,6)
    plot(x,Bim,x,imag(B_ideal))
    title('B - imag')
    subplot(4,2,7)
    plot(x,abs(A),x,abs(A_ideal))
    title('A - magnitude')
    subplot(4,2,8)
    plot(x,abs(B),x,abs(B_ideal))
    title('B - magnitude')
else
    ideal = fft(re_in + 1j * im_in);
    half_fft = N/2;
    if reorder_freq                                       %frequencies are needed to be [neg, 0, pos]
        ideal = fftshift(ideal);
    else                                                  %frequencies are [0, pos, neg]
    end
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
    subplot(3,2,4)
    plot(x,Xim,x,imag(ideal))
    title('X - imag')
    subplot(3,2,5)
    plot(x,abs(X),x,abs(ideal))
    title('X - Magnitude')
end
