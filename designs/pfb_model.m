function[theoretical_prefilter_output, theoretical_fft_output_a,...
    theoretical_fft_output_b, theoretical_fft_output_c] = pfb_model(in_re, in_im, nof_points, dual_processing, reorder_freq, nof_taps, fwidth, enable_prefilter)
%pfb_model - Description
%
% This function approximately models the behaviour of the ASTRON pfb block in the CASPER
% toolflow. It uses the MATLAB prefilter and wideband FFT functions with floating point precision and is to be used as an 
% ideal.
    if enable_prefilter
    theoretical_prefilter_output_re = prefilter_model(in_re, nof_points, nof_taps, fwidth);
    theoretical_prefilter_output_im = prefilter_model(in_im, nof_points, nof_taps, fwidth);
    else
        theoretical_prefilter_output_re = NaN;
        theoretical_prefilter_output_im = NaN;
    end
    [output_a,output_b,output_c] = wideband_fft_model(in_re, in_im ...
                                ,nof_points, dual_processing, reorder_freq);

    theoretical_fft_output_a = output_a;
    theoretical_fft_output_b = output_b; 
    theoretical_fft_output_c = output_c;
    theoretical_prefilter_output = theoretical_prefilter_output_re + 1j .* theoretical_prefilter_output_im;
end