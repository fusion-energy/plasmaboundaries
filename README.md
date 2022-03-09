[![CI with install](https://github.com/fusion-energy/plasmaboundaries/actions/workflows/ci_with_install.yml/badge.svg)](https://github.com/fusion-energy/plasmaboundaries/actions/workflows/ci_with_install.yml)

[![Upload Python Package](https://github.com/fusion-energy/plasmaboundaries/actions/workflows/python-publish.yml/badge.svg)](https://github.com/fusion-energy/plasmaboundaries/actions/workflows/python-publish.yml)
[![conda-publish](https://anaconda.org/conda-forge/plasmaboundaries/badges/version.svg)](https://anaconda.org/conda-forge/plasmaboundaries)
[![PyPI version](https://badge.fury.io/py/plasmaboundaries.svg)](https://badge.fury.io/py/plasmaboundaries)

# plasma-boundaries

This code computes and plots analytical solutions of the Grad-Shafranov (GS)
equation for studying plasma equilibrium, stability and transport in fusion
reactors based on the work of A. Cerfon and J. Freidberg.
Details on the method can be found in [*"One size fits all" analytical solutions to the Grad-Shafranov equation*, Physics of Plasmas 17 (2010)](https://doi.org/10.1063/1.3328818)

Documentation can be found [here](https://plasma-boundaries.readthedocs.io/en/latest/).

## Installation (Conda)
You can install plasma-boundaries using [Conda](https://docs.conda.io/en/latest/) by running:
```bash
conda install -c fusion-energy -c conda-forge plasmaboundaries
```

## Installation (pip)
You can install plasma-boundaries using [Pip](https://pip.pypa.io/en/stable/) by running:
```bash
pip install plasmaboundaries
```

Alternatively you can pip install from the repository:
```bash
pip install git+https://github.com/fusion-energy/plasmaboundaries.git
```


## Usage

First compute the magnetic flux <img src="https://render.githubusercontent.com/render/math?math=\Psi"> from plasma-boundaries based on a specific set of parameters.
In this example, the built-in ITER plasma parameters will be used:
```python
import plasmaboundaries

# plasma parameters
params = plasmaboundaries.ITER

# compute magnetic flux psi(R, z)
psi = plasmaboundaries.compute_psi(params, config='double-null')
```

The magnetic flux can now be calculated for any coordinates and plotted with matplotlib:
```python
print(psi(1.0, 0))

# plot the results
import matplotlib.pyplot as plt
import numpy as np

rmin, rmax = 0.6, 1.4
zmin, zmax = -0.6, 0.6
r = np.arange(rmin, rmax, step=0.01)
z = np.arange(zmin, zmax, step=0.01)
R, Z = np.meshgrid(r, z)
PSI = psi(R, Z)  # compute magnetic flux

levels = np.linspace(PSI.min(), 0, num=25)
CS = plt.contourf(R, Z, PSI, levels=levels, vmax=0)
plt.contour(R, Z, PSI, levels=[0], colors="black") # display the separatrix

plt.colorbar(CS, label="Magnetic flux $\Psi$")
plt.xlabel('Radius $R/R_0$')
plt.ylabel('Height $z/R_0$')
plt.gca().set_aspect("equal")
plt.show()
```
In `compute_psi`, the argument `config` can also be set to `'single-null'` or `'non-null'` for other plasma shapes.

<img src="https://user-images.githubusercontent.com/40028739/87403291-f8fbda80-c5bc-11ea-971e-7856043855de.png">
<img src="https://user-images.githubusercontent.com/40028739/87404184-1c735500-c5be-11ea-93a3-16ed588bf3c6.png">

### Custom plasma parameters
Parameters can also be defined by creating the parameters dictionary:
```python
params = {
    "A": -0.155,
    "aspect_ratio": 0.32,
    "elongation": 1.7,
    "triangularity": 0.33,
}
```

## Run the tests

You can run the tests with:
```pytest tests```
