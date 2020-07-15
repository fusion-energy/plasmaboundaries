from plasmaboundaries.magnetic_flux import psi_up_down_symmetric, \
    psi_up_down_asymmetric, derivatives

import numpy as np
from scipy.optimize import fsolve
import sympy as sp


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
    A = params["A"]
    aspect_ratio = params["aspect_ratio"]
    triangularity = params["triangularity"]
    elongation = params["elongation"]
    N_1, N_2, N_3 = params["N_1"], params["N_2"], params["N_3"]

    if config == "non-null" or config == "double-null":
        psi_ = psi_up_down_symmetric
    else:
        psi_ = psi_up_down_asymmetric

    def psi(x, y):
        return psi_(x, y, p, A, pkg='sp')
    psi_x_sp, psi_y_sp = derivatives(psi, 1)
    psi_xx_sp, psi_yy_sp = derivatives(psi, 2)

    psi_x = val_from_sp(psi_x_sp)
    psi_y = val_from_sp(psi_y_sp)
    psi_xx = val_from_sp(psi_xx_sp)
    psi_yy = val_from_sp(psi_yy_sp)

    arguments = [psi, (psi_x, psi_y), (psi_xx, psi_yy), A, aspect_ratio, triangularity, elongation, (N_1, N_2, N_3)]
    if config == "non-null":
        list_of_equations = constraints_non_null(*arguments)
    elif config == "single-null":
        list_of_equations = constraints_single_null(*arguments)
    elif config == "double-null":
        list_of_equations = constraints_double_null(*arguments)
    return list_of_equations


def constraints_non_null(
        f, first_order_d, second_order_d, A, aspect_ratio,
        triangularity, elongation, N_coeffs):
    """Creates set of constraints for parametric GS solution non-null
    plasma configuration.

    Args:
        f (callable): function f(x, y)
        first_order_d ((callable, callable)): first order derivatives of the
        function f
        second_order_d ((callable, callable)): second order derivatives of the
        function f: (d^2/dy^2(f), d^2/dx^2(f))
        A (float): Plasma parameter
        aspect_ratio (float): Plasma parameter
        triangularity (float): Plasma parameter
        elongation (float): Plasma parameter
        N_coeffs ((float, float, float)): Coefficients N_1, N_2, N_3 based on
        plasma
        parameters

    Returns:
        list: set of constraints
    """
    fx, fy = first_order_d
    fxx, fyy = second_order_d
    N_1, N_2, N_3 = N_coeffs

    list_of_equations = [
        f(1 + aspect_ratio, 0),
        f(1 - aspect_ratio, 0),
        f(1 - triangularity*aspect_ratio, elongation*aspect_ratio),
        fx(1 - triangularity*aspect_ratio, elongation*aspect_ratio),
        fyy(1 + aspect_ratio, 0) + N_1*fx(1 + aspect_ratio, 0),
        fyy(1 - aspect_ratio, 0) + N_2*fx(1 - aspect_ratio, 0),
        fxx(1 - triangularity*aspect_ratio, elongation*aspect_ratio) +
        N_3*fy(1 - triangularity*aspect_ratio, elongation*aspect_ratio)
        ]
    return list_of_equations


def constraints_single_null(
        f, first_order_d, second_order_d, A, aspect_ratio,
        triangularity, elongation, N_coeffs):
    """Creates set of constraints for parametric GS solution single-null
    plasma configuration.

    Args:
        f (callable): function f(x, y)
        first_order_d ((callable, callable)): first order derivatives of the
        function f
        second_order_d ((callable, callable)): second order derivatives of the
        function f: (d^2/dy^2(f), d^2/dx^2(f))
        A (float): Plasma parameter
        aspect_ratio (float): Plasma parameter
        triangularity (float): Plasma parameter
        elongation (float): Plasma parameter
        N_coeffs ((float, float, float)): Coefficients N_1, N_2, N_3 based on
        parameters

    Returns:
        list: set of constraints
    """
    fx, fy = first_order_d
    fxx, fyy = second_order_d
    N_1, N_2, N_3 = N_coeffs

    x_sep, y_sep = 1-1.1*triangularity*aspect_ratio, -1.1*elongation*aspect_ratio
    list_of_equations = [
        f(1 + aspect_ratio, 0),
        f(1 - aspect_ratio, 0),
        f(1 - triangularity*aspect_ratio, elongation*aspect_ratio),
        f(x_sep, y_sep),
        fy(1 + aspect_ratio, 0),
        fy(1 - aspect_ratio, 0),
        fx(1 - triangularity*aspect_ratio, elongation*aspect_ratio),
        fx(x_sep, y_sep),
        fy(x_sep, y_sep),
        fyy(1 + aspect_ratio, 0) + N_1*fx(1 + aspect_ratio, 0),
        fyy(1 - aspect_ratio, 0) + N_2*fx(1 - aspect_ratio, 0),
        fxx(1 - triangularity*aspect_ratio, elongation*aspect_ratio) +
        N_3*fy(1 - triangularity*aspect_ratio, elongation*aspect_ratio),
        ]
    return list_of_equations


def constraints_double_null(
        f, first_order_d, second_order_d, A, aspect_ratio,
        triangularity, elongation, N_coeffs):
    """Creates set of constraints for parametric GS solution double-null
    plasma configuration.

    Args:
        f (callable): function f(x, y)
        first_order_d ((callable, callable)): first order derivatives of the
        function f
        second_order_d ((callable, callable)): second order derivatives of the
        function f: (d^2/dy^2(f), d^2/dx^2(f))
        A (float): Plasma parameter
        aspect_ratio (float): Plasma parameter
        triangularity (float): Plasma parameter
        elongation (float): Plasma parameter
        N_coeffs ((float, float, float)): Coefficients N_1, N_2, N_3 based on
        plasma
        parameters

    Returns:
        list: set of constraints
    """
    fx, fy = first_order_d
    fxx, fyy = second_order_d
    N_1, N_2, N_3 = N_coeffs

    x_sep, y_sep = 1-1.1*triangularity*aspect_ratio, 1.1*elongation*aspect_ratio
    list_of_equations = [
        f(1 + aspect_ratio, 0),
        f(1 - aspect_ratio, 0),
        f(x_sep, y_sep),
        fx(x_sep, y_sep),
        fy(x_sep, y_sep),
        fyy(1 + aspect_ratio, 0) + N_1*fx(1 + aspect_ratio, 0),
        fyy(1 - aspect_ratio, 0) + N_2*fx(1 - aspect_ratio, 0),
        ]
    return list_of_equations


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


def compute_coefficients_c_i(params, constraints, config):
    """Calculates the coefficients c_i based on a set of constraints

    Args:
        params (dict): contains the plasma parameters
        (aspect_ratio, elongation, triangularity, A)
        constraints (callable): list of equations
        coefficient_number (int): number of constraints/coefficients
         (7 if up-down symetrical, 12 if up-down asymmetrical)

    Returns:
        list: coefficients c_i (floats)
    """
    N_1, N_2, N_3 = compute_N_i(params)
    params["N_1"], params["N_2"], params["N_3"] = N_1, N_2, N_3

    if config == "non-null" or config == "double-null":
        coefficient_number = 7
    elif config == "single-null":
        coefficient_number = 12
    x_0 = np.zeros(coefficient_number) + 1

    coefficients = fsolve(constraints, x_0, args=(params, config))
    return coefficients


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
    if config == "non-null" or config == "double-null":
        psi = psi_up_down_symmetric
    elif config == "single-null":
        psi = psi_up_down_asymmetric

    coefficients = compute_coefficients_c_i(
        params, constraints=constraints,
        config=config)

    def new_psi(X, Y, pkg='np'):
        return psi(
            X, Y, c_i=coefficients, A=params["A"], pkg=pkg)

    if return_coeffs:
        return new_psi, coefficients
    else:
        return new_psi
