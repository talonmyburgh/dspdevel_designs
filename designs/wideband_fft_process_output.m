function [output_a, output_b] = wideband_fft_process_output(wideband_factor, nof_points, dual_processing)
    pipelined_nof_points = nof_points/wideband_factor;

        %complex out index
    if ~dual_processing
        Xre = zeros(wideband_factor,pipelined_nof_points);
        Xim = zeros(wideband_factor,pipelined_nof_points);
        f_tmp = 0:wideband_factor:(pipelined_nof_points-1)*wideband_factor;
        f = zeros(wideband_factor,pipelined_nof_points);
        for p=1:1:wideband_factor
            f(p,:) = p + f_tmp;
            Xre(p,:) = re_out(f(p,:));
            Xim(p,:) = im_out(f(p,:));
        end
        output_a=Xre + 1j*Xim;
    %real out index
    else
        even = 2:2:pipelined_nof_points;
        odd = 1:2:pipelined_nof_points;
        Are = zeros(wideband_factor,pipelined_nof_points/2);
        Aim = zeros(wideband_factor,pipelined_nof_points/2);
        Bre = zeros(wideband_factor,pipelined_nof_points/2);
        Bim = zeros(wideband_factor,pipelined_nof_points/2);
        
        for p=1:1:wideband_factor
            Are(p,:) = re_out(p,even);
            Bre(p,:) = re_out(p,odd);
            Aim(p,:) = im_out(p,even);
            Bim(p,:) = im_out(p,odd);
        end
        Are = reshape(Are,[1,nof_points/2]);
        Aim = reshape(Aim,[1,nof_points/2]);
        Bre = reshape(Bre,[1,nof_points/2]);
        Bim = reshape(Bim,[1,nof_points/2]);
        output_a = Are+1j*Aim;
        output_b = Bre+1j*Bim;
    end
end