import plasmaboundaries
import numpy as np
from scipy.optimize import fsolve
import sympy as sp
import matplotlib.pyplot as plt


SHIFT = 0.1  # a 10% shift is assumed for computing the X point coordinates


def test_points(aspect_ratio, elongation, triangularity):
    """Compute the coordinates of inner and outer equatorial points and high
    point based on plasma geometrical parameters.

    Args:
        aspect_ratio (float): minor radius / major radius
        elongation (float): plasma elongation
        triangularity (float): plasma triangularity

    Returns:
        ((float, float), (float, float), (float, float)): points (x, y)
            coordinates
    """
    outer_equatorial_point = (1 + aspect_ratio, 0)
    inner_equatorial_point = (1 - aspect_ratio, 0)
    high_point = (1 - triangularity*aspect_ratio, elongation*aspect_ratio)
    return outer_equatorial_point, inner_equatorial_point, high_point


def val_from_sp(expression):
    """Transforms a sympy expression to a callable function
    f(x, y)

    Args:
        expression (sympy.Add): sympy expression to be converted
            which has symbols 'x' and 'y' in it.
    """
    def val(x, y):
        return expression.subs('x', x).subs('y', y)
    return val


def constraints(p, params, config):
    """Creates set of constraints for parametric GS solution for any
    plasma configuration.

    Args:
        p (list): c_i coefficients (floats)
        params (dict): contains the plasma parameters
            (aspect_ratio, elongation, triangularity, A)
        config (str): shape of the plasma 'non-null', 'single-null',
            'double-null'.

    Returns:
        list: set of constraints
    """

    # create sympy expressions for derivatives
    def psi(x, y):
        return plasmaboundaries.psi(x, y, p, params["A"], config, pkg='sp')
    psi_x_sp, psi_y_sp = plasmaboundaries.derivatives(psi, 1)
    psi_xx_sp, psi_yy_sp = plasmaboundaries.derivatives(psi, 2)

    psi_x, psi_y = val_from_sp(psi_x_sp), val_from_sp(psi_y_sp)
    psi_xx, psi_yy = val_from_sp(psi_xx_sp), val_from_sp(psi_yy_sp)

    # create test points
    aspect_ratio = params["aspect_ratio"]
    triangularity = params["aspect_ratio"]
    elongation = params["elongation"]
    outer_equatorial_point, inner_equatorial_point, high_point = \
        test_points(aspect_ratio, elongation, triangularity)

    # constraints common to all configurations
    N_1, N_2, N_3 = params["N_1"], params["N_2"], params["N_3"]

    list_of_constraints = [
        psi(*outer_equatorial_point),
        psi(*inner_equatorial_point),
        psi_yy(*outer_equatorial_point) + N_1*psi_x(*outer_equatorial_point),
        psi_yy(*inner_equatorial_point) + N_2*psi_x(*inner_equatorial_point),
        ]

    # add constraints common to non-null and single-null
    if config in ["non-null", "single-null"]:
        list_of_constraints += [
            psi_xx(*high_point) + N_3*psi_y(*high_point),
            psi(*high_point),
            psi_x(*high_point),
            ]

    # add constraints for single-null only
    if config in ["single-null", "double-null"]:
        x_sep, y_sep = 1-(1+SHIFT)*triangularity*aspect_ratio, \
            -(1+SHIFT)*elongation*aspect_ratio
        list_of_constraints += [
            psi(x_sep, y_sep),
            psi_x(x_sep, y_sep),
            psi_y(x_sep, y_sep),
        ]
    if config == "single-null":
        list_of_constraints += [
            psi_y(*outer_equatorial_point),
            psi_y(*inner_equatorial_point),
        ]
    return list_of_constraints


def compute_N_i(params):
    """Computes N_1, N_2 and N_3 coefficients based on plasma parameters

    Args:
        params (dict): contains the plasma parameters
            (aspect_ratio, elongation, triangularity, A)

    Returns:
        (float, float, float): (N_1, N_2, N_3)
    """
    triangularity = params["triangularity"]
    aspect_ratio = params["aspect_ratio"]
    elongation = params["elongation"]

    alpha = np.arcsin(triangularity)  # alpha
    N_1 = -(1 + alpha)/(aspect_ratio*elongation**2)
    N_2 = (1 - alpha)/(aspect_ratio*elongation**2)
    N_3 = -elongation/(aspect_ratio*elongation*np.cos(alpha)**2)

    return N_1, N_2, N_3


def compute_psi(params, config="non-null", return_coeffs=False):
    """Computes the magnetic flux fonction

    Args:
        params (dict): contains the plasma parameters
            (aspect_ratio, elongation, triangularity, A)
        config (str, optional): shape of the plasma
            "non-null", "single-null", "double-null". Defaults to "non-null".
        return_coeffs (bool, optional): If True, will also return the
            coefficients c_i. Defaults to False.

    Returns:
        (callable) or (callable, list): Magnetic flux fonction and
            coefficients c_i (only if return_coeffs is True)
    """
    params["N_1"], params["N_2"], params["N_3"] = compute_N_i(params)

    if config in ["non-null", "double-null"]:
        coefficient_number = 7
    elif config == "single-null":
        coefficient_number = 12
    x_0 = np.zeros(coefficient_number)
    coefficients = fsolve(constraints, x_0, args=(params, config), xtol=10)

    def new_psi(X, Y, pkg='np'):
        return plasmaboundaries.psi(
            X, Y, c_i=coefficients, A=params["A"], config=config, pkg=pkg)

    if return_coeffs:
        return new_psi, coefficients
    else:
        return new_psi


def get_separatrix_coordinates(params, config, step=0.01):
    """Creates a list of points describing the separatrix

    Args:
        params (dict): contains the plasma parameters
            (aspect_ratio, elongation, triangularity, A)
        config (str): shape of the plasma
            "non-null", "single-null", "double-null".
        step (float, optional): Resolution of the domain. Defaults to 0.01.

    Raises:
        ValueError: If no separatrix is found within the points of interest

    Returns:
        numpy.array: list of points coordinates
    """

    aspect_ratio = params["aspect_ratio"]
    elongation = params["elongation"]
    triangularity = params["triangularity"]
    inner_equatorial_point = (1 - aspect_ratio, 0)
    outer_equatorial_point = (1 + aspect_ratio, 0)
    if config == "non-null" or config == "double-null":
        if config == "non-null":
            low_point = (1 - aspect_ratio*triangularity, -elongation*aspect_ratio)
        else:
            low_point = (1 - 1.1*triangularity*aspect_ratio, -1.1*elongation*aspect_ratio)
        high_point = (low_point[0], -low_point[1])
    elif config == "single-null":
        low_point = (1 - 1.1*triangularity*aspect_ratio, -1.1*elongation*aspect_ratio)
        high_point = (1 - aspect_ratio*triangularity, elongation*aspect_ratio)

    xmin, xmax = inner_equatorial_point[0], outer_equatorial_point[0]
    ymin, ymax = low_point[1], high_point[1]

    security_factor = 1.1
    x = np.arange(xmin/security_factor, xmax*security_factor, step=step)
    y = np.arange(ymin*security_factor, ymax*security_factor, step=step)

    X, Y = np.meshgrid(x, y)

    # compute psi
    psi = compute_psi(params, config=config)
    Z = psi(X, Y)  # compute magnetic flux

    # add contour
    separatrix = plt.contour(X, Y, Z, levels=[0], colors="white", linestyles="dashed")

    # extract points coordinates
    points = np.empty((0, 2))
    for p in separatrix.collections[0].get_paths():
        v = p.vertices
        points = np.append(points, v, 0)

    # remove contours
    for cont in separatrix.collections:
        cont.remove()
    return points
