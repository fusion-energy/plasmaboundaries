# plasma-boundaries
# This script is an implementation of the method described in
# “One size fits all” analytic solutions to the Grad–Shafranov equation
# A. J. Cerfon and J. P. Freidberg, Physics of Plamas 17 032502 (2010)
# https://doi.org/10.1063/1.3328818

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


import numpy as np
import plasma_boundaries

# plasma parameters
params = [plasma_boundaries.NSTX_single_null, plasma_boundaries.NSTX_double_null]

fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(10, 4.8))


for ax, config, param in zip(
        [ax1, ax2],
        ["single-null", "double-null"],
        params):

    # compute psi
    psi = plasma_boundaries.compute_psi(param, config=config)

    # plot the results
    xmin, xmax = 0.1, 2
    ymin, ymax = -2, 2
    x = np.arange(xmin, xmax, step=0.01)
    y = np.arange(ymin, ymax, step=0.01)

    X, Y = np.meshgrid(x, y)
    Z = psi(X, Y)  # compute magnetic flux

    # add filled contours
    levels2 = np.unique(np.linspace(Z.min(), -Z.min(), num=100))
    norm = mcolors.TwoSlopeNorm(vmin=Z.min(), vcenter=0., vmax=-Z.min())
    CSF = ax.contourf(
        X, Y, Z, levels=levels2, cmap="coolwarm", norm=norm,
        vmax=-Z.min(), extend="max")

    # add contours
    levels = np.unique(
        np.append(
            np.linspace(Z.min(), 0, num=10), np.linspace(0, Z.max(), num=20)))
    CS = ax.contour(
        X, Y, Z, levels=levels[levels != 0], colors="black",
        linestyles="solid")
    separatrix = ax.contour(
        X, Y, Z, levels=[0], colors="black", linestyles="dashed")
    ax.clabel(separatrix, inline=True, fmt=r"$\Psi = $%.0f")

    ax.set_title('NSTX ' + config)
    ax.set_xlabel('Radius $R/R_0$')
    ax.set_aspect("equal")
ax1.set_ylabel('Height $Z/R_0$')
plt.colorbar(CSF, label="Magnetic flux $\Psi$", format="%.3f")
plt.show()
