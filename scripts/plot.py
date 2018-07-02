import numpy as np

from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

import matplotlib.pyplot as plt

lines = open("log.out", "r").readlines()
tags = lines[0].rstrip().split()
print (tags)
mat = []
for i in range(1, len(lines)):
  mat.append(map(float, lines[i].rstrip().split()))

print (mat)

mat = np.array(mat)
dists = squareform(mat)
linkage_matrix = linkage(dists, "single")
dendrogram(linkage_matrix, labels=tags, orientation='right')
plt.title("test")
plt.show()
