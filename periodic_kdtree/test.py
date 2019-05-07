#!/usr/bin/env python3

from periodic_kdtree import PeriodicCKDTree
import numpy as np

# Boundaries (0 or negative means open boundaries in that dimension)
bounds = np.array([30.0, 30.0, -1])   # xy periodic, open along z

# Points
n = 10000
data = 30.0 * np.random.randn(n, 3)

# Build kd-tree
T = PeriodicCKDTree(bounds, data)

# Find 4 closest neighbors to a random point
# (d[j], i[j]) = distance and index of jth closest point
d, i = T.query([45.0, 10.0, 10.0], k=4)

# Find neighbors within a fixed distance of a point
neighbors = T.query_ball_point([45.0, 10.0, 10.0], r=3.0)

print(neighbors)
