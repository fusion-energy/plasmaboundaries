import numpy as np


ITER = {
    "aspect_ratio": 0.32,  # a/R_0
    "A": -0.155,  # A, arbitrary ?
    "elongation": 1.7,  # kappa
    "triangularity": 0.33,  # delta
}

spheromak_equilibrium = {
    "aspect_ratio": 0.95,
    "A": 1,
    "elongation": 1,
    "triangularity": 0.2
}

NSTX_double_null = {
    "aspect_ratio": 0.78,
    "A": 0,
    "elongation": 2,
    "triangularity": 0.35
}

NSTX_single_null = {
    "aspect_ratio": 0.78,
    "A": -0.05,
    "elongation": 2,
    "triangularity": 0.35
}
