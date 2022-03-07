import plasmaboundaries

import sympy as sp
import numpy as np


def GS_equation(psi, A, x, y):
    """Grad-Shafranov equation x*d(psi_x/x)/dx + psi_yy - (1_A)*x**2 - A = 0

    Args:
        psi (sympy.Add): Function psi(x, y) (surface flux)
        A (float): plasma parameter such as A = R_0**4/Psi_0**2 * F * dF/dpsi
        x (sp.symbol): x-coordinate
        y (sp.symbol): y-coordinate

    Returns:
        sympy.Add: left hand-side of GS equation
    """
    psi_x = sp.diff(psi, x)  # x first order derivative
    psi_yy = sp.diff(psi, y, y)  # y second order derivative
    eq = x*sp.diff((1/x)*psi_x, x) + psi_yy
    eq += - (1 - A)*x**2 - A
    return eq


def test_non_null():
    """Test the function compute_psi with non-null configuration
    """
    A = -0.05
    params = {
        "A": A,
        "aspect_ratio": 0.3,
        "triangularity": 0.7,
        "elongation": 2
    }
    config = "non-null"
    psi = plasmaboundaries.compute_psi(params, config=config)

    x, y = sp.symbols("x y")

    psi_sp = psi(x, y, pkg='sympy')
    x_min, x_max = 0.001, 2
    y_min, y_max = -2, 2

    x_test, y_test = \
        np.linspace(x_min, x_max, num=50), \
        np.linspace(y_min, y_max, num=50)
    for point_x, point_y in zip(x_test, y_test):
        val = GS_equation(psi_sp, A, x, y).subs(x, point_x).subs(y, point_y)
        assert np.isclose(float(val), 0)


def test_single_null():
    """Test the function compute_psi with single-null configuration
    """
    A = 0
    params = {
        "A": A,
        "aspect_ratio": 0.3,
        "triangularity": 0.7,
        "elongation": 2
    }
    config = "single-null"
    psi = plasmaboundaries.compute_psi(params, config=config)

    x, y = sp.symbols("x y")

    psi_sp = psi(x, y, pkg='sympy')

    x_min, x_max = 0.001, 2
    y_min, y_max = -2, 2

    x_test, y_test = \
        np.linspace(x_min, x_max, num=50), \
        np.linspace(y_min, y_max, num=50)
    for point_x, point_y in zip(x_test, y_test):

        val = GS_equation(psi_sp, A, x, y).subs(x, point_x).subs(y, point_y)
        assert np.isclose(float(val), 0)


def test_double_null():
    """Test the function compute_psi with double-null configuration
    """
    A = 0
    params = {
        "A": A,
        "aspect_ratio": 0.3,
        "triangularity": 0.7,
        "elongation": 2
    }
    config = "double-null"
    psi = plasmaboundaries.compute_psi(params, config=config)

    x, y = sp.symbols("x y")

    psi_sp = psi(x, y, pkg='sympy')

    x_min, x_max = 0.001, 2
    y_min, y_max = -2, 2

    x_test, y_test = \
        np.linspace(x_min, x_max, num=50), \
        np.linspace(y_min, y_max, num=50)
    for point_x, point_y in zip(x_test, y_test):
        val = GS_equation(psi_sp, A, x, y).subs(x, point_x).subs(y, point_y)
        assert np.isclose(float(val), 0)
