# plasma-boundaries
# This script is an implementation of the method described in
# “One size fits all” analytic solutions to the Grad–Shafranov equation
# A. J. Cerfon and J. P. Freidberg, Physics of Plamas 17 032502 (2010)
# https://doi.org/10.1063/1.3328818

import matplotlib.pyplot as plt
import numpy as np
from plasma_boundaries import parameters, compute_psi

# plasma parameters
params = parameters.ITER

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True, figsize=(10, 4.8))
for ax, config in zip([ax1, ax2, ax3], ["non-null", "single-null", "double-null"]):
    # compute psi
    psi = compute_psi(params, config=config)

    # plot the results
    xmin, xmax = 0.6, 1.35
    ymin, ymax = -0.8, 0.7
    x = np.arange(xmin, xmax, step=0.01)
    y = np.arange(ymin, ymax, step=0.01)

    X, Y = np.meshgrid(x, y)
    Z = psi(X, Y)  # compute magnetic flux

    # ax.title("ITER single-null")
    levels = np.linspace(Z.min(), 0, num=25)
    CS = ax.contourf(X, Y, Z, levels=levels, vmax=0)
    ax.contour(X, Y, Z, levels=[0], colors="black")
    ax.set_title('ITER ' + config)
    ax.set_xlabel('Radius $R/R_0$')
    ax.set_aspect("equal")
ax1.set_ylabel('Height $Z/R_0$')
fig.colorbar(CS, label="Magnetic flux $\Psi$", format="%.3f")
plt.show()
