import numpy as np

def compute_N_i(parameters):
    triangularity = parameters["triangularity"]
    epsilon = parameters["epsilon"]
    elongation = parameters["elongation"]

    alpha = np.arcsin(triangularity)  # alpha
    N_1 = -(1 + alpha)/(epsilon*elongation**2)
    N_2 = (1 - alpha)/(epsilon*elongation**2)
    N_3 = -elongation/(epsilon*elongation*np.cos(alpha)**2)

    return N_1, N_2, N_3


ITER = {
    "epsilon": 0.32,  # a/R_0
    "A": -0.155,  # A, arbitrary ?
    "elongation": 1.7,  # kappa
    "triangularity": 0.33,  # delta
}

spheromak_equilibrium = {
    "epsilon": 0.95,
    "A": 1,
    "elongation": 1,
    "triangularity": 0.2
}

NSTX_double_null = {
    "epsilon": 0.78,
    "A": 0,
    "elongation": 2,
    "triangularity": 0.35
}

NSTX_single_null = {
    "epsilon": 0.78,
    "A": -0.05,
    "elongation": 2,
    "triangularity": 0.35
}
