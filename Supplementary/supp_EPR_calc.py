
import matplotlib.pyplot as plt
import numpy as np

cov = np.array([[0, -75.0319, -84.6724, 80.883],
               [-81.7206, 0, 87.801, -94.4123],
               [-77.063, 83.5677, 0, -73.0699],
               [82.8215, -75.061, -78.0624, 0]])

cov = cov + 94.4123
cov[0, 0] = 0
cov[1, 1] = 0
cov[2, 2] = 0
cov[3, 3] = 0

# cov = cov.reshape((2, 2, 2, 2))


def sum_coordinates(cov, H=4, W=4):
    sum_cov = np.zeros((H * 2, W * 2))
    for y1 in range(H):
        for x1 in range(W):
                for y2 in range(H):
                    for x2 in range(W):
                            if x1 == x2 and y1 == y2:
                                continue
                            sum_cov[y1 + y2, x1 + x2] += (cov[y1, x1] + cov[y2, x2])
    return sum_cov




