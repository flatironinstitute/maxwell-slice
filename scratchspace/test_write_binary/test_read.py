import numpy as np

with open('dd.bin', 'rb') as f:
    # skip the beginning bytes
    # discard = np.fromfile(f, count=1, dtype=np.int32)
    x = np.fromfile(f, dtype=np.float64)

# fid = fopen('dd.bin','rb');
# discard = fread(fid,[1 1],'*int32'); % discard the beginning bytes
# x = fread(fid,[4],'*double'); % read everything in
# x
# fclose(fid);
print(x)