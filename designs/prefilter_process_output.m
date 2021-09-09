function prefilter_output = prefilter_process_output(output_array, wideband_factor, nof_points)
    %%Interleave recorded output%%
    interleaved_output = zeros(size(output_array,1)*wideband_factor, 1);
    for i=1:1:wideband_factor
        interleaved_output(i:wideband_factor:end) = output_array(:,i);
    end
    %%Determine index at which output-data is valid%%
    valid_output_length = size(interleaved_output);
    whole_output_last_index = (valid_output_length - mod(valid_output_length,nof_points));
    prefilter_output = reshape(interleaved_output(1:whole_output_last_index-1), nof_points, []);
end