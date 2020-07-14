from plasmaboundaries.magnetic_flux import psi_up_down_symmetric, \
    psi_up_down_asymmetric, derivatives

import numpy as np
from scipy.optimize import fsolve
import sympy as sp


def val_from_sp(expression):
    def val(x, y):
        return expression.subs('x', x).subs('y', y)
    return val


def constraints(p, params, config):

    A = params["A"]
    epsilon = params["epsilon"]
    triangularity = params["triangularity"]
    elongation = params["elongation"]
    N_1, N_2, N_3 = params["N_1"], params["N_2"], params["N_3"]

    if config == "non-null" or config == "double-null":
        psi_ = psi_up_down_symmetric
    else:
        psi_ = psi_up_down_asymmetric

    def psi(x, y):
        return psi_(x, y, p, A, pkg=sp)
    psi_x_sp, psi_y_sp = derivatives(psi, 1)
    psi_xx_sp, psi_yy_sp = derivatives(psi, 2)

    psi_x = val_from_sp(psi_x_sp)
    psi_y = val_from_sp(psi_y_sp)
    psi_xx = val_from_sp(psi_xx_sp)
    psi_yy = val_from_sp(psi_yy_sp)

    arguments = [psi, (psi_x, psi_y), (psi_xx, psi_yy), A, epsilon, triangularity, elongation, (N_1, N_2, N_3)]
    if config == "non-null":
        list_of_equations = constraints_non_null(*arguments)
    elif config == "single-null":
        list_of_equations = constraints_single_null(*arguments)
    elif config == "double-null":
        list_of_equations = constraints_double_null(*arguments)
    return list_of_equations


def constraints_non_null(
        f, first_order_d, second_order_d, A, epsilon,
        triangularity, elongation, N_coeffs):

    fx, fy = first_order_d
    fxx, fyy = second_order_d
    N_1, N_2, N_3 = N_coeffs

    list_of_equations = [
        f(1 + epsilon, 0),
        f(1 - epsilon, 0),
        f(1 - triangularity*epsilon, elongation*epsilon),
        fx(1 - triangularity*epsilon, elongation*epsilon),
        fyy(1 + epsilon, 0) + N_1*fx(1 + epsilon, 0),
        fyy(1 - epsilon, 0) + N_2*fx(1 - epsilon, 0),
        fxx(1 - triangularity*epsilon, elongation*epsilon) +
        N_3*fy(1 - triangularity*epsilon, elongation*epsilon)
        ]
    return list_of_equations


def constraints_single_null(
        f, first_order_d, second_order_d, A, epsilon,
        triangularity, elongation, N_coeffs):

    fx, fy = first_order_d
    fxx, fyy = second_order_d
    N_1, N_2, N_3 = N_coeffs

    x_sep, y_sep = 1-1.1*triangularity*epsilon, -1.1*elongation*epsilon
    list_of_equations = [
        f(1 + epsilon, 0),
        f(1 - epsilon, 0),
        f(1 - triangularity*epsilon, elongation*epsilon),
        f(x_sep, y_sep),
        fy(1 + epsilon, 0),
        fy(1 - epsilon, 0),
        fx(1 - triangularity*epsilon, elongation*epsilon),
        fx(x_sep, y_sep),
        fy(x_sep, y_sep),
        fyy(1 + epsilon, 0) + N_1*fx(1 + epsilon, 0),
        fyy(1 - epsilon, 0) + N_2*fx(1 - epsilon, 0),
        fxx(1 - triangularity*epsilon, elongation*epsilon) +
        N_3*fy(1 - triangularity*epsilon, elongation*epsilon),
        ]
    return list_of_equations


def constraints_double_null(
        f, first_order_d, second_order_d, A, epsilon,
        triangularity, elongation, N_coeffs):

    fx, fy = first_order_d
    fxx, fyy = second_order_d
    N_1, N_2, N_3 = N_coeffs

    x_sep, y_sep = 1-1.1*triangularity*epsilon, 1.1*elongation*epsilon
    list_of_equations = [
        f(1 + epsilon, 0),
        f(1 - epsilon, 0),
        f(x_sep, y_sep),
        fx(x_sep, y_sep),
        fy(x_sep, y_sep),
        fyy(1 + epsilon, 0) + N_1*fx(1 + epsilon, 0),
        fyy(1 - epsilon, 0) + N_2*fx(1 - epsilon, 0),
        ]
    return list_of_equations


def compute_N_i(params):
    """Computes the N_1 N_2 and N_3 coefficients based on plasma parameters

    Args:
        params (dict): contains the plasma params
        (epsilon, elongation, triangularity, A)

    Returns:
        (float, float, float): (N_1, N_2, N_3)
    """
    triangularity = params["triangularity"]
    epsilon = params["epsilon"]
    elongation = params["elongation"]

    alpha = np.arcsin(triangularity)  # alpha
    N_1 = -(1 + alpha)/(epsilon*elongation**2)
    N_2 = (1 - alpha)/(epsilon*elongation**2)
    N_3 = -elongation/(epsilon*elongation*np.cos(alpha)**2)

    return N_1, N_2, N_3


def compute_coefficients_c_i(params, constraints, config):
    """calculates the coefficients c_i based on a set of constraints

    Args:
        params (dict): contains the plasma params
        (epsilon, elongation, triangularity, A)
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
    """Compute the magnetic flux fonction

    Args:
        params (dict): contains the plasma parameters
        config (str, optional): shape of the plasma
         "non-null", "single-null", "double-null". Defaults to "non-null".
        return_coeffs (bool, optional): If True, will also return the
         coefficients c_i. Defaults to False.

    Returns:
        (callable) or (callable, list): Magnetic flux fonction and
         coefficients c_i (if return_coeffs is True)
    """
    if config == "non-null" or config == "double-null":
        psi = psi_up_down_symmetric
    elif config == "single-null":
        psi = psi_up_down_asymmetric

    coefficients = compute_coefficients_c_i(
        params, constraints=constraints,
        config=config)

    def new_psi(X, Y, pkg=np):
        return psi(
            X, Y, coefficients_c=coefficients, A=params["A"], pkg=pkg)

    if return_coeffs:
        return new_psi, coefficients
    else:
        return new_psi
