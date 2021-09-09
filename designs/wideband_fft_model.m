function [output_a, output_b] = wideband_fft_model(in_re, in_im, nof_points, dual_processing, reorder_freq)
    %complex out index
    if dual_processing
        ideal = fft(in_re + 1j * in_im);
        R_ideal = real(ideal);
        R_flip = R_ideal(:);
        R_flip(2:end) = R_flip(end:-1:2);
        
        I_ideal = imag(ideal);
        I_flip = I_ideal(:);
        I_flip(2:end) = I_flip(end:-1:2);
        
        output_b = (1/2) .* (R_ideal + 1j * I_ideal + R_flip - 1j * I_flip);
        output_b = output_b(1:nof_points/2);
        output_a = (1/(2*1j)) * (R_ideal + 1j * I_ideal - R_flip + 1j * I_flip);
        output_a = output_a(1:nof_points/2);

    else
        ideal = fft(in_re + 1j * in_im);
        if reorder_freq                                       %frequencies are needed to be [neg, 0, pos]
            ideal = fftshift(ideal);
        else                                                  %frequencies are [0, pos, neg]
        end
        output_a = ideal;
        output_b = [];
    end
end