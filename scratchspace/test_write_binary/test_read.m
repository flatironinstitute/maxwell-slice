fid = fopen('dd.bin','rb');
x = fread(fid,[4],'*double'); % read everything in
x
fclose(fid);
