from . import magnetic_flux
from . import parameters

import numpy as np
from scipy.optimize import fsolve


def constraints(p, parameters):
    """Creates all the constraints for numerical solving

    Args:
        p (list): coefficients c_i

    Returns:
        list: list of constraints
    """
    A = parameters["A"]
    epsilon = parameters["epsilon"]
    triangularity = parameters["triangularity"]
    elongation = parameters["elongation"]
    N_1, N_2, N_3 = parameters["N_1"], parameters["N_2"], parameters["N_3"]
    psi = magnetic_flux.psi_up_down_symmetric
    psi_x_sp, psi_y_sp = magnetic_flux.derivatives(psi, p, A, 1)
    psi_xx_sp, psi_yy_sp = magnetic_flux.derivatives(psi, p, A, 2)

    def psi_x(x, y):
        return psi_x_sp.subs('x', x).subs('y', y)

    def psi_y(x, y):
        return psi_y_sp.subs('x', x).subs('y', y)

    def psi_xx(x, y):
        return psi_xx_sp.subs('x', x).subs('y', y)

    def psi_yy(x, y):
        return psi_yy_sp.subs('x', x).subs('y', y)

    list_of_equations = [
        psi(1 + epsilon, 0, p, A),
        psi(1 - epsilon, 0, p, A),
        psi(1 - triangularity*epsilon, elongation*epsilon, p, A),
        psi_x(1 - triangularity*epsilon, elongation*epsilon),
        psi_yy(1 + epsilon, 0) + N_1*psi_x(1 + epsilon, 0),
        psi_yy(1 - epsilon, 0) + N_2*psi_x(1 - epsilon, 0),
        psi_xx(1 - triangularity*epsilon, elongation*epsilon) +
        N_3*psi_y(1 - triangularity*epsilon, elongation*epsilon)
        ]
    return list_of_equations


def constraints_single_null(p, parameters):
    """Creates all the constraints for numerical solving in
     the case of a double null divertor configuration

    Args:
        p (list): coefficients c_i

    Returns:
        list: list of constraints
    """
    A = parameters["A"]
    epsilon = parameters["epsilon"]
    triangularity = parameters["triangularity"]
    elongation = parameters["elongation"]
    N_1, N_2, N_3 = parameters["N_1"], parameters["N_2"], parameters["N_3"]
    psi = magnetic_flux.psi_up_down_asymetric
    psi_x_sp, psi_y_sp = magnetic_flux.derivatives(psi, p, A, 1)
    psi_xx_sp, psi_yy_sp = magnetic_flux.derivatives(psi, p, A, 2)

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
        psi(1 + epsilon, 0, p, A),
        psi(1 - epsilon, 0, p, A),
        psi(1 - triangularity*epsilon, elongation*epsilon, p, A),
        psi(x_sep, y_sep, p, A),
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


def constraints_double_null(p, parameters):
    """Creates all the constraints for numerical solving in
     the case of a double null divertor configuration

    Args:
        p (list): coefficients c_i
        psi (callable): function phi to use
        parameters (dict): contains parameters

    Returns:
        list: list of constraints
    """
    A = parameters["A"]
    epsilon = parameters["epsilon"]
    triangularity = parameters["triangularity"]
    elongation = parameters["elongation"]
    N_1, N_2, N_3 = parameters["N_1"], parameters["N_2"], parameters["N_3"]
    psi = magnetic_flux.psi_up_down_symmetric
    psi_x_sp, psi_y_sp = magnetic_flux.derivatives(psi, p, A, 1)
    psi_xx_sp, psi_yy_sp = magnetic_flux.derivatives(psi, p, A, 2)

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
        psi(1 + epsilon, 0, p, A),
        psi(1 - epsilon, 0, p, A),
        psi(x_sep, y_sep, p, A),
        psi_x(x_sep, y_sep),
        psi_y(x_sep, y_sep),
        psi_yy(1 + epsilon, 0) + N_1*psi_x(1 + epsilon, 0),
        psi_yy(1 - epsilon, 0) + N_2*psi_x(1 - epsilon, 0),
        ]
    return list_of_equations


def compute_coefficients_c_i(params, constraints, coefficient_number):
    """calculates the coefficients c_i based on a set of constraints

    Args:
        params (dict): contains the plasma parameters
        (epsilon, elongation, triangularity, A)
        constraints (callable): list of equations
        coefficient_number (int): number of constraints/coefficients
         (7 if up-down symetrical, 12 if up-down asymetrical)

    Returns:
        list: coefficients c_i (floats)
    """
    N_1, N_2, N_3 = parameters.compute_N_i(params)
    params["N_1"], params["N_2"], params["N_3"] = N_1, N_2, N_3
    x_0 = np.zeros(coefficient_number) + 1

    coefficients = fsolve(constraints, x_0, args=(params))
    return coefficients


def compute_psi(params, config="non-null", return_coeffs=False):

    if config == "non-null" or config == "double-null":
        psi = magnetic_flux.psi_up_down_symmetric
        if config == "non-null":
            constraints_ = constraints
        else:
            constraints_ = constraints_double_null
        coefficients = compute_coefficients_c_i(
            params, constraints=constraints_, coefficient_number=7)
    else:
        constraints_ = constraints_single_null
        psi = magnetic_flux.psi_up_down_asymetric
        coefficients = compute_coefficients_c_i(
            params, constraints=constraints_, coefficient_number=12)

    def new_psi(X, Y, pkg=np):
        return psi(X, Y, coefficients_c=coefficients, A=params["A"], pkg=pkg)

    if return_coeffs:
        return new_psi, coefficients
    else:
        return new_psi
