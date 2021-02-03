fid = fopen('dd.bin','rb');
discard = fread(fid,[1 1],'*int32'); % discard the beginning bytes
x = fread(fid,[4],'*double'); % read everything in
x
fclose(fid);