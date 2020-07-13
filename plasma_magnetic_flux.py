# Magnetic flux fonction
# This script is an implementation of the method described in
# “One size fits all” analytic solutions to the Grad–Shafranov equation
# A. J. Cerfon and J. P. Freidberg, Physics of Plamas 17 032502 (2010)
# https://doi.org/10.1063/1.3328818

import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from scipy.optimize import fsolve


def psi(X, Y, coefficients_c, pkg=np):
    """returns the value of magnetic flux at point (X, Y)
     according to coefficients ci

    Args:
        X (float or numpy.array): x coordinate
        Y (float or numpy.array): y coordinate
        coefficients_c (list): list of floats, the ci coefficients
        pkg (callable, optional): if set to np (resp. sp), numpy (resp. sympy)
         objects will be used. Defaults to np.

    Returns:
        float or numpy.array: value(s) of magnetic flux
    """
    psi_1 = 1
    psi_2 = X**2
    psi_3 = Y**2 - X**2*pkg.log(X)
    psi_4 = X**4 - 4*X**2*Y**2
    psi_5 = 2*Y**4 - 9*Y**2*X**2 + 3*X**4*pkg.log(X) - 12*X**2*Y**2*pkg.log(X)
    psi_6 = X**6 - 12*X**4*Y**2 + 8*X**2*Y**4
    psi_7 = 8*Y**6 - 140*Y**4*X**2 + 75*Y**2*X**4 - 15*X**6*pkg.log(X) + \
        180*X**4*Y**2*pkg.log(X) - 120*X**2*Y**4*pkg.log(X)

    psis = [psi_1, psi_2, psi_3, psi_4, psi_5, psi_6, psi_7]
    val = X**4/8 + A*(1/2*X**2*pkg.log(X) - X**4/8) + \
        sum([coefficients_c[i]*psis[i] for i in range(len(coefficients_c))])
    return val


def psi_up_down_asymetric(X, Y, coefficients_c, pkg=np):
    """returns the value of magnetic flux at point (X, Y)
     according to coefficients ci

    Args:
        X (float or numpy.array): x coordinate
        Y (float or numpy.array): y coordinate
        coefficients_c (list): list of floats, the ci coefficients
        pkg (callable, optional): if set to np (resp. sp), numpy (resp. sympy)
         objects will be used. Defaults to np.

    Returns:
        float or numpy.array: value(s) of magnetic flux
    """
    
    psi_1 = 1
    psi_2 = X**2
    psi_3 = Y**2 - X**2*pkg.log(X)
    psi_4 = X**4 - 4*X**2*Y**2
    psi_5 = 2*Y**4 - 9*Y**2*X**2 + 3*X**4*pkg.log(X) - 12*X**2*Y**2*pkg.log(X)
    psi_6 = X**6 - 12*X**4*Y**2 + 8*X**2*Y**4
    psi_7 = 8*Y**6 - 140*Y**4*X**2 + 75*Y**2*X**4 - 15*X**6*pkg.log(X) + \
        180*X**4*Y**2*pkg.log(X) - 120*X**2*Y**4*pkg.log(X)
    psi_8 = Y
    psi_9 = Y*X**2
    psi_10 = Y**3 - 3*Y*X**2*pkg.log(X)
    psi_11 = 3*Y*X**4 - 4*Y**3*X**2
    psi_12 = 8*Y**5 - 45*Y*X**4 - 80*Y**3*X**2*pkg.log(X) + \
        60*Y*X**4*pkg.log(X)

    psis = [
        psi_1, psi_2, psi_3, psi_4, psi_5, psi_6, psi_7, psi_8, psi_9,
        psi_10, psi_11, psi_12]

    val = X**4/8 + A*(1/2*X**2*pkg.log(X) - X**4/8) + \
        sum([coefficients_c[i]*psis[i] for i in range(len(coefficients_c))])
    return val


def derivatives(f, coefficients_c, order):
    """Computes the derivatives of magnetic flux.
    Does not computes xy or yx derivatives.

    Args:
        f (callable f(x, y, coefficients_c, pkg)): function of magnetic flux
        coefficients_c (list): coefficients ci
        order (int): order of differenciation

    Returns:
        list: [x_derivative, y_derivative]
    """
    x, y = sp.symbols("x y")
    f_ = f(X=x, Y=y, coefficients_c=coefficients_c, pkg=sp)
    f_x = sp.diff(f_, *[x for i in range(order)])
    f_y = sp.diff(f_, *[y for i in range(order)])
    return [f_x, f_y]


# Plasma parameters

# ITER
epsilon = 0.32  # a/R_0
minor_radius = 1  # a
major_radius = minor_radius/epsilon  # R_0
A = -0.155  # A, arbitrary ?
elongation = 1.7  # kappa
triangularity = 0.33  # delta
# # Spheromak
# epsilon = 0.95
# elongation = 1
# triangularity = 0.2
# A = 1


# NSTX
epsilon = 0.78
elongation = 2
triangularity = 0.35
A = 0  # use A = -0.05 for single null configuration ; A = 0 else

alpha = np.arcsin(triangularity)  # alpha

N_1 = -(1 + alpha)/(epsilon*elongation**2)
N_2 = (1 - alpha)/(epsilon*elongation**2)
N_3 = -elongation/(epsilon*elongation*np.cos(alpha)**2)


def constraints(p):
    """Creates all the constraints for numerical solving

    Args:
        p (list): coefficients c_i

    Returns:
        list: list of constraints
    """
    psi_x_sp, psi_y_sp = derivatives(psi, p, 1)
    psi_xx_sp, psi_yy_sp = derivatives(psi, p, 2)

    def psi_x(x, y):
        return psi_x_sp.subs('x', x).subs('y', y)

    def psi_y(x, y):
        return psi_y_sp.subs('x', x).subs('y', y)

    def psi_xx(x, y):
        return psi_xx_sp.subs('x', x).subs('y', y)

    def psi_yy(x, y):
        return psi_yy_sp.subs('x', x).subs('y', y)

    list_of_equations = [
        psi(1 + epsilon, 0, p),
        psi(1 - epsilon, 0, p),
        psi(1 - triangularity*epsilon, elongation*epsilon, p),
        psi_x(1 - triangularity*epsilon, elongation*epsilon),
        psi_yy(1 + epsilon, 0) + N_1*psi_x(1 + epsilon, 0),
        psi_yy(1 - epsilon, 0) + N_2*psi_x(1 - epsilon, 0),
        psi_xx(1 - triangularity*epsilon, elongation*epsilon) +
        N_3*psi_y(1 - triangularity*epsilon, elongation*epsilon)
        ]
    return list_of_equations


def constraints_single_null(p):
    """Creates all the constraints for numerical solving in
     the case of a double null divertor configuration

    Args:
        p (list): coefficients c_i

    Returns:
        list: list of constraints
    """
    psi_x_sp, psi_y_sp = derivatives(psi_up_down_asymetric, p, 1)
    psi_xx_sp, psi_yy_sp = derivatives(psi_up_down_asymetric, p, 2)

    def psi_x(x, y):
        return psi_x_sp.subs('x', x).subs('y', y)

    def psi_y(x, y):
        return psi_y_sp.subs('x', x).subs('y', y)

    def psi_xx(x, y):
        return psi_xx_sp.subs('x', x).subs('y', y)

    def psi_yy(x, y):
        return psi_yy_sp.subs('x', x).subs('y', y)
    x_sep, y_sep = 1-1.1*triangularity*epsilon, -1.1*elongation*epsilon
    list_of_equations = [
        psi_up_down_asymetric(1 + epsilon, 0, p),
        psi_up_down_asymetric(1 - epsilon, 0, p),
        psi_up_down_asymetric(1 - triangularity*epsilon, elongation*epsilon, p),
        psi_up_down_asymetric(x_sep, y_sep, p),
        psi_y(1 + epsilon, 0),
        psi_y(1 - epsilon, 0),
        psi_x(1 - triangularity*epsilon, elongation*epsilon),
        psi_x(x_sep, y_sep),
        psi_y(x_sep, y_sep),

        psi_yy(1 + epsilon, 0) + N_1*psi_x(1 + epsilon, 0),
        psi_yy(1 - epsilon, 0) + N_2*psi_x(1 - epsilon, 0),
        psi_xx(1 - triangularity*epsilon, elongation*epsilon) +
        N_3*psi_y(1 - triangularity*epsilon, elongation*epsilon),
        ]
    return list_of_equations


def constraints_double_null(p):
    """Creates all the constraints for numerical solving in
     the case of a double null divertor configuration

    Args:
        p (list): coefficients c_i

    Returns:
        list: list of constraints
    """
    psi_x_sp, psi_y_sp = derivatives(psi, p, 1)
    psi_xx_sp, psi_yy_sp = derivatives(psi, p, 2)

    def psi_x(x, y):
        return psi_x_sp.subs('x', x).subs('y', y)

    def psi_y(x, y):
        return psi_y_sp.subs('x', x).subs('y', y)

    def psi_xx(x, y):
        return psi_xx_sp.subs('x', x).subs('y', y)

    def psi_yy(x, y):
        return psi_yy_sp.subs('x', x).subs('y', y)

    x_sep, y_sep = 1-1.1*triangularity*epsilon, 1.1*elongation*epsilon
    list_of_equations = [
        psi(1 + epsilon, 0, p),
        psi(1 - epsilon, 0, p),
        psi(x_sep, y_sep, p),
        psi_x(x_sep, y_sep),
        psi_y(x_sep, y_sep),

        psi_yy(1 + epsilon, 0) + N_1*psi_x(1 + epsilon, 0),
        psi_yy(1 - epsilon, 0) + N_2*psi_x(1 - epsilon, 0),
        ]
    return list_of_equations


# compute coefficients
x_0 = np.zeros(7) + 1

constraints_ = constraints_double_null
psi_ = psi

coefficients = fsolve(constraints_, x_0)

# plot the results
xmin, xmax = 0, 2
ymin, ymax = -2, 2
step = 0.01
x = np.arange(xmin, xmax, step=step)
y = np.arange(ymin, ymax, step=step)

X, Y = np.meshgrid(x, y)
Z = psi_(X, Y, coefficients)  # compute magnetic flux

plt.figure()

# levels = np.linspace(-0.15, 0.15, num=15)
CS = plt.contourf(X, Y, Z, levels=15)
plt.contour(X, Y, Z, levels=[0], colors="white")
plt.colorbar(CS, label="Magnetic flux $\Psi$")
plt.xlabel('Radius $R/R_0$')
plt.ylabel('Height $Z/R_0$')
plt.axis('square')
plt.show()
