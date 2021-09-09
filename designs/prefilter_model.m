function theoretical_prefilter_output = prefilter_model(input_signal, nof_points, nof_taps, fwidth)
        %prefilter_model - Description
        %
        % This function approximately models the behaviour of the ASTRON prefilter block in the CASPER
        % toolflow. It uses MATLAB functions as well as floating point precision and is to be used as an 
        % ideal.

        signal_length = length(input_signal);
        filter_length = nof_points*nof_taps;
        %%Generate theoretical window%%
        window_coeffs = hanning(filter_length);
        filter_coeffs = window_coeffs .* transpose(sinc(fwidth*((1:filter_length)/nof_points - (nof_taps/2))));

        %Figure out how many sets of input phases/output phases there will be.
        slice_count = signal_length/nof_points;
        %reshape the input signal into slices that are nof_points long
        input_signal = reshape(input_signal,nof_points,slice_count);
        zero_padding = zeros(nof_points,nof_taps);
        input_signal = horzcat(zero_padding,input_signal);
        %reshape the coefficients into nof_points x nof_taps
        phased_filter = reshape(filter_coeffs, nof_points, nof_taps);
        %outputs
        theoretical_prefilter_output = zeros(nof_points, slice_count);
        for slice_index = 1:slice_count-1
                %Commutate input-slice across P rows (resulting in 'nof_taps' cols)%%
                phased_input = input_signal(:,slice_index:slice_index+nof_taps-1);
                theoretical_prefilter_output(:, slice_index) = sum(phased_input .* phased_filter, 2); % sum across taps (dim=2)
                %sum(reshape(1:6, 2, 3) .* [1, 1, 1; 2, 2, 2], 2) % proof of the above computation    
        end
end