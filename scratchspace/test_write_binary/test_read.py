import numpy as np

with open('dd.bin', 'rb') as f:
    x = np.fromfile(f, dtype=np.float64)

print(x)
