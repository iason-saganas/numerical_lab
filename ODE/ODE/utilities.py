import numpy as np


def derivative_of_exp_example_1(x, y):
    """
    Represents a derivative of an example differential equation:

        y' = -y

    The function signature has to have the independent variable as first argument, even if not used,
    for the solvers to work properly.

    :param x: np.array,  The independent variable
    :param y: np.array,  The dependent variable
    """
    return -y


def solution_to_exp_example_1(y):
    """
    The solution to y in case the derivative is given by the function `derivative_of_exp_example_1`
    :param y: np.array, The dependent variable range to calculate the solution for.
    """
    return np.exp(-y)
