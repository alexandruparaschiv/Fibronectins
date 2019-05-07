#!/usr/bin/env python3


import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

img = Image.open("./fibrils.png")
plt.imshow(img)
plt.show()


w=10
h=10
fig=plt.figure(figsize=(8, 8))
columns = 4
rows = 5
for i in range(1, columns*rows +1):
    img = Image.open("./fibrils.png")
    fig.add_subplot(rows, columns, i)
    plt.imshow(img)
plt.show()
