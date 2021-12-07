from scipy.io import wavfile
import matplotlib.pyplot as plt
from matplotlib import image
from matplotlib import pyplot
from PIL import Image
from numpy import asarray
import sys
import numpy as np

# data = asarray(Image.open(sys.argv[1]))
im = plt.imread(sys.argv[1])
data = im[:, :, 0]
with open(sys.argv[2], 'w') as f:
    print(len(data), file=f)
    print(len(data[0]), file=f)
    for i in data:
        for j in i:
            print(j, file=f, end=' ')
        print('', file=f)
