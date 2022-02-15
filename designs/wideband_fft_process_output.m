function [output_a, output_b, output_x] = wideband_fft_process_output(re_out, im_out, wideband_factor, nof_points, dual_processing)
    % re_out, im_out has wb_factor columns 
    pipelined_nof_points = nof_points/wideband_factor;
    num_output_slices = floor(size(re_out,1)/nof_points);
    output_a = zeros(nof_points/2, num_output_slices);
    output_b = zeros(nof_points/2, num_output_slices);
    output_x = zeros(nof_points, num_output_slices);
    for slice_index = 1:num_output_slices
        % Process pipelined_nof_points at a time where each is [pipelined_nof_points x wideband_factor].
        % Things aren't actually interleaved across wideband_factor, so collect across rows first
        slice_start = 1+pipelined_nof_points*(slice_index-1);
        re_output_slice = re_out(slice_start:slice_start+pipelined_nof_points-1,:); 
        im_output_slice = im_out(slice_start:slice_start+pipelined_nof_points-1,:); 

        % Real output
        if dual_processing
            odd = 1:2:pipelined_nof_points;
            even = 2:2:pipelined_nof_points;
            Are = zeros(pipelined_nof_points/2, wideband_factor);
            Aim = zeros(pipelined_nof_points/2, wideband_factor);
            Bre = zeros(pipelined_nof_points/2, wideband_factor);
            Bim = zeros(pipelined_nof_points/2, wideband_factor);
            
            for p=1:1:wideband_factor
                Are(:, p) = re_output_slice(odd, p);
                Bre(:, p) = re_output_slice(even, p);
                Aim(:, p) = im_output_slice(odd, p);
                Bim(:, p) = im_output_slice(even, p);
            end
            Are = reshape(Are,[nof_points/2,1]);
            Aim = reshape(Aim,[nof_points/2,1]);
            Bre = reshape(Bre,[nof_points/2,1]);
            Bim = reshape(Bim,[nof_points/2,1]);
            output_a(:,slice_index) = Are+1j*Aim;
            output_b(:,slice_index) = Bre+1j*Bim;
        % Complex output
        else
            Xre = zeros(wideband_factor,pipelined_nof_points);
            Xim = zeros(wideband_factor,pipelined_nof_points);
            f_tmp = 0:wideband_factor:(pipelined_nof_points-1)*wideband_factor;
            f = zeros(wideband_factor,pipelined_nof_points);
            for p=1:1:wideband_factor
                f(p,:) = p + f_tmp;
                Xre(p,:) = re_output_slice(f(p,:));
                Xim(p,:) = im_output_slice(f(p,:));
            end
            output_x(:,slice_index) = Xre + 1j*Xim;

        end
    end
end