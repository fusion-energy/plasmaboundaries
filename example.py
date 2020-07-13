# Magnetic flux fonction
# This script is an implementation of the method described in
# “One size fits all” analytic solutions to the Grad–Shafranov equation
# A. J. Cerfon and J. P. Freidberg, Physics of Plamas 17 032502 (2010)
# https://doi.org/10.1063/1.3328818

import matplotlib.pyplot as plt
import numpy as np
from plasma_boundaries import parameters, model

import numpy as np

# plasma parameters
params = parameters.ITER


# compute psi
psi = model.compute_psi(params, config="double-null")

# plot the results
xmin, xmax = 0, 2
ymin, ymax = -2, 2
step = 0.01
x = np.arange(xmin, xmax, step=step)
y = np.arange(ymin, ymax, step=step)

X, Y = np.meshgrid(x, y)
Z = psi(X, Y)  # compute magnetic flux

plt.figure()

CS = plt.contourf(X, Y, Z, levels=15)
plt.contour(X, Y, Z, levels=[0], colors="white")
plt.colorbar(CS, label="Magnetic flux $\Psi$")
plt.xlabel('Radius $R/R_0$')
plt.ylabel('Height $Z/R_0$')
plt.axis('square')
plt.show()
