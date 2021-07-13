% pinc for dds compiler sim
fout = 15.625;
fsys = 156.25;
pwidth = 29;

pinc_d = (fout*2^pwidth)/(fsys);
pinc = pinc_d/2^pwidth

% system reset array
arr = zeros(1, 5000);
arr(2:1:end) = 1;
%arr(50) = 0;
a = timeseries(arr);