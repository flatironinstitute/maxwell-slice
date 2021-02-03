import numpy as np

with open('data.bin', 'rb') as f:
    x = np.fromfile(f, dtype=np.float64)

n = int(len(x) / 7)
x = np.reshape(x, (n, 7))
y = np.zeros((n, 4), dtype=np.complex)
y[:, 0] = x[:, 0] + 1j * x[:, 1]
y[:, 1] = x[:, 2] + 1j * x[:, 3]
y[:, 2] = x[:, 4] + 1j * x[:, 5]
y[:, 3] = x[:, 6]

print(y)
print(y.shape)