function [output_a, output_b, output_x] = wideband_fft_model(in_re, in_im, nof_points, dual_processing, reorder_freq)
    input_signal_length = size(in_re)-mod(size(in_re),nof_points);
    re_input_sliced = reshape(in_re(1:input_signal_length),nof_points,[]);
    im_input_sliced = reshape(in_im(1:input_signal_length),nof_points,[]);
    num_input_slices = size(re_input_sliced,2);
    output_a = zeros(nof_points/2,num_input_slices);
    output_b = zeros(nof_points/2,num_input_slices);
    output_x = zeros(nof_points,num_input_slices);
    %complex out index
    for slice_index = 1 : num_input_slices
        if dual_processing
            ideal = fft(re_input_sliced(:,slice_index) + 1j * im_input_sliced(:,slice_index));
            R_ideal = real(ideal);
            R_flip = R_ideal(:);
            R_flip(2:end) = R_flip(end:-1:2);
            
            I_ideal = imag(ideal);
            I_flip = I_ideal(:);
            I_flip(2:end) = I_flip(end:-1:2);
            
            tmp_output_a = (1/2) .* (R_ideal + 1j * I_ideal + R_flip - 1j * I_flip);
            output_a(:,slice_index) = tmp_output_a(1:nof_points/2);
            tmp_output_b = (1/(2*1j)) * (R_ideal + 1j * I_ideal - R_flip + 1j * I_flip);
            output_b(:,slice_index) = tmp_output_b(1:nof_points/2);
        else
            ideal = fft(re_input_sliced(:,slice_index) + 1j * im_input_sliced(:,slice_index));
            if reorder_freq                                       %frequencies are needed to be [neg, 0, pos]
                ideal = fftshift(ideal);
            else                                                  %frequencies are [0, pos, neg]
            end
            output_x(:,slice_index) = ideal;
        end

    end
end