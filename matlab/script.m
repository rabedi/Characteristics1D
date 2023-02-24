dwh = datawheader;
fid = fopen('sd_0_tPos_20000.txt', 'r');
dwh = dwh.read(fid);
fclose(fid);

dwh.plotRaw(3);
fclose('all');
close('all');
%pause;

fid = fopen('sd_0__Summary.txt', 'r');
dwh = dwh.read(fid);
fclose(fid);

dwh.plotRaw(2);
fclose('all');
close('all');
a = 12;

