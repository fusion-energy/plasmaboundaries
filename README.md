# plasma-boundaries

This code computes and plots analytical solutions of the Grad-Shafranov (GS) equation for studying plasma equilibrium, stability and transport in fusion reactors based on the work of A. Cerfon and J. Freidberg.
Details on the method can be found in [*"One size fits all" analytical solutions to the Grad-Shafranov equation*, Physics of Plasmas 17 (2010)](https://doi.org/10.1063/1.3328818)

## Installation

First clone this repository
```git clone https://github.com/RemiTheWarrior/plasma-boundaries```

Install the dependencies
```pip install -r requirements.txt```

## Usage

First compute the magnetic flux <img src="https://render.githubusercontent.com/render/math?math=\Psi"> from plasma-boundaries based on a specific set of parameters.
In this example, the built-in ITER plasma parameters will be used:
```python
import plasma_boundaries

# plasma parameters
params = plasma_boundaries.ITER

# compute magnetic flux psi(x, y)
psi = plasma_boundaries.compute_psi(params, config='double-null')
```

The magnetic flux can now be calculated for any coordinates and ploted with matplotlib:
```python
print(psi(1.0, 0))

# plot the results
import matplotlib.pyplot as plt
import numpy as np

xmin, xmax = 0.6, 1.4
ymin, ymax = -0.6, 0.6
x = np.arange(xmin, xmax, step=0.01)
y = np.arange(ymin, ymax, step=0.01)
X, Y = np.meshgrid(x, y)
Z = psi(X, Y)  # compute magnetic flux

levels = np.linspace(Z.min(), 0, num=25)
CS = plt.contourf(X, Y, Z, levels=levels, vmax=0)
plt.contour(X, Y, Z, levels=[0], colors="black") # display the separatrix

plt.colorbar(CS, label="Magnetic flux $\Psi$")
plt.xlabel('Radius $R/R_0$')
plt.ylabel('Height $Z/R_0$')
plt.gca().set_aspect("equal")
plt.show()
```
In `compute_psi`, the argument `config` can also be set to `'single-null'` or `'non-null'` for other plasma shapes.

### Custom plasma parameters
Parameters can also be defined by creating the parameters dictionary:
```python
params = {
    "A": -0.155,
    "epsilon": 0.32,
    "elongation": 1.7,
    "triangularity": 0.33,
}
```
<img src="https://user-images.githubusercontent.com/40028739/87342222-31fb6700-c54b-11ea-87ab-ac2f57f7b113.png">
