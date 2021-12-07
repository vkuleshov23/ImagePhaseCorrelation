import matplotlib.pyplot as plt
from matplotlib import image
from matplotlib import pyplot
from numpy import asarray
import numpy as np


file1 = "tank_first.bmp"
file2 = "tank_second.bmp"
# file2 = "tank_first.bmp"
# file1 = "tank_second.bmp"

first = plt.imread(file1);
second = plt.imread(file2);

r1 = first[:, :, 0]
r2 = second[:, :, 0]

U1 = np.fft.fft2(r1)
U2 = np.fft.fft2(r2)

tmp = U1 * U2.conjugate()

arr = abs(np.fft.ifft2(tmp.T)).T

plt.matshow(abs(arr))
plt.colorbar()
plt.show()

indices = np.where(arr == arr.max())
print(indices[0], indices[1])
# plt.plot(spectr1)
