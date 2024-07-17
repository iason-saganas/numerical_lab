import numpy as np
import matplotlib.pyplot as plt
from styles.matplotlib_style import *


# First, we define some test functions


def f1(x):
    # Note: The integral of this function from 0 to pi/2 is 1.
    return np.sin(x)


def f2(x):
    # Note: The integral of this function from 0 to pi is 2.
    return np.sin(x)


def f3(x):
    # Note: The integral of this function from 0 to 1 is 0.746824132812427.
    return np.exp(-x ** 2)


def f4(x):
    # Note: The integral of this function from 0 to 4 is 1.32581766366803.
    return 1 / (1 + x ** 2)


# Attach the correct values to the functions
f1.correct_val = 1
f2.correct_val = 2
f3.correct_val = 0.746824132812427
f4.correct_val = 1.32581766366803


def trapez_method(func: callable, a: float, b: float, n: int) -> float:
    """

    :param func:    Callable,       The function to integrate.
    :param a:       float,          The lower bound of the integration.
    :param b:       float,          The upper bound of the integration.
    :param n:       int,            The number of intervals to partition the integration interval into.
    :return:
    """
    # Calculate the step length
    h = (b - a) / n

    # Code up array containing weights
    w = np.ones(n+1)  # because the real w_i goes from w_0, w_n (length: n+1) so
    # the code w_i should go from w_0 to w_n+1 (length: n+1)
    w[0] = w[-1] = 1 / 2

    interval_areas = []
    for i in range(n+1):
        # iterate over all n+1 weights
        interval_area = w[i] * func(a + i * h)
        interval_areas.append(interval_area)

    area = h * np.sum(np.array(interval_areas))
    return area


def simpson_method(func: callable, a: float, b: float, n: int) -> float:
    """

    :param func:    Callable,       The function to integrate.
    :param a:       float,          The lower bound of the integration.
    :param b:       float,          The upper bound of the integration.
    :param n:       int,            The number of intervals to partition the integration interval into.
    :return:
    """

    if n % 2 != 0:
        raise ValueError("Number n of intervals must be even for simpson integration!")

    # Calculate the step length
    h = (b - a) / n

    # Code up array containing weights
    w = np.ones(n+1)
    w[0] = w[-1] = 1 / 3
    for i in range(1, n):
        # Start point: w[0], end point: w[n]
        # starts i=1 and ends at n-1 => Start and end point excluded
        is_even = (i % 2 == 0)
        if is_even:
            w[i] = 2 / 3
        else:
            w[i] = 4 / 3

    interval_areas = []
    for i in range(n+1):
        interval_area = w[i] * func(a + i * h)
        interval_areas.append(interval_area)

    area = h * np.sum(np.array(interval_areas))
    return area


def custom_analyze(func: callable, method: callable, a: float, b: float, print_table=False):
    """
    `method` must be one of `simpson_method` or `trapez_method`.
    `func` must be a python callable that has an attribute `correct_val` attached to it, that represents
     the correct (possibly analytic) value for the integral over some assumed interval.

     Returns an array of areas, errors, and error ratios.


    For parameter documentation see docstrings of functions `simpson_method` and `trapez_method`.
    This function takes the method and the integrand provided and iterates the numerical integration,
    varying the interval number as powers of two.

    Using powers of two is convenient for the interpretation of the error ratio, which is the quotient of
    error in step n-1 and step n.

    If this ratio is big, the error in step n is smaller than the error in step n-1, meaning we got closer to
    the solution.

    """

    areas = []
    errors = []
    ratios = []
    if print_table:
        print("\nPrinting table for method ", method)
    for index, j in enumerate(range(1, 12)):
        n = 2 ** j
        area = method(func, a, b, n)
        error = np.abs(func.correct_val - area)
        errors.append(error)
        areas.append(area)
        ratio = errors[index - 1] / errors[index]
        ratios.append(ratio)
        if print_table:
            print(f"n={n}:", f"area={area}", f"error={error}", f"ratio={ratio}")
    return np.array([areas, errors, ratios])


def make_error_plot(errors, ratios, names_of_methods, latex_title):
    """
    Use like:

    x = np.array([2**j for j in range(1,12)])
    trapez_areas, trapez_errors, trapez_ratio = table_I1_trapez
    simpson_areas, simpson_errors, simpson_ratio = table_I1_simpson

    make_error_plot(errors = [simpson_errors, trapez_errors],
                ratios = [simpson_ratio, trapez_ratio],
                names_of_methods=["Simpson", "Trapez"],
                latex_title=r"$\int_0^{\pi/2} \hspace{1mm} \mathrm{sin}(x)\hspace{1mm} \mathrm{d}x$")

    :param errors:  The absolute value of the absolute difference between an approximation at step n and the true value.
    :param ratios:  The error ratio between step n-1 and step n.
    :param names_of_methods: An array consisting of strings of the names of the methods that are compared
    :param latex_title: The latex code to be rendered as the plot
    representing the integral currently under investigation.
    :return:
    """

    x = np.array([2**j for j in range(1,12)])

    plt.subplot(2, 1, 1)
    for i, ratio in enumerate(ratios):
        plt.plot(x, ratio, ".-", label=names_of_methods[i])
    plt.ylabel(r"Error ratio $\tilde{\varepsilon}$")
    plt.legend()

    plt.subplot(2, 1, 2)
    for i, error in enumerate(errors):
        plt.plot(x, error, ".-", label=names_of_methods[i])
    plt.xlabel("Step $n$")
    plt.ylabel(r"Error $\varepsilon$")
    plt.yscale("log")
    plt.legend()

    plt.suptitle(latex_title, fontsize=20)
    plt.tight_layout()


def make_area_plot(areas, correct_value, names_of_methods, latex_title):
    """
    Use like:

    areas = [areas_trapez, areas_simpson]

    :param areas:   An array consisting of arrays that represent the n-iteration under different methods.
    :param correct_value: The correct integral value.
    :param names_of_methods: An array consisting of strings of the names of the methods that are compared
    :param latex_title: The latex code to be rendered as the plot
    representing the integral currently under investigation.
    :return:
    """
    x = np.array([2 ** j for j in range(1, 12)])

    for i, area in enumerate(areas):
        plt.plot(x, area, ".-", label=names_of_methods[i])

    plt.hlines(correct_value, 0, max(x), ls="--", color="black", label="Correct value")

    plt.xlabel("Step $n$")
    plt.ylabel("Value of integral")
    plt.yscale("log")
    plt.xscale("log")
    plt.title(latex_title)
    plt.legend()
    plt.tight_layout()


table_I1_trapez = custom_analyze(func=f1, method=trapez_method, a=0, b=np.pi/2, print_table=True)
table_I1_simpson = custom_analyze(func=f1, method=simpson_method, a=0, b=np.pi/2, print_table=True)

table_I3_trapez = custom_analyze(func=f3, method=trapez_method, a=0, b=1, print_table=True)
table_I3_simpson = custom_analyze(func=f3, method=simpson_method, a=0, b=1, print_table=True)

table_I4_trapez = custom_analyze(func=f4, method=trapez_method, a=0, b=4, print_table=True)
table_I4_simpson = custom_analyze(func=f4, method=simpson_method, a=0, b=4, print_table=True)

# Plot a comparison plot of the trapez and simpson method integrating the function f1

x = np.array([2**j for j in range(1,12)])
trapez_areas, trapez_errors, trapez_ratio = table_I1_trapez
simpson_areas, simpson_errors, simpson_ratio = table_I1_simpson

make_error_plot(errors = [simpson_errors, trapez_errors],
                ratios = [simpson_ratio, trapez_ratio],
                names_of_methods=["Simpson", "Trapez"],
                latex_title=r"$\int_0^{\pi/2} \hspace{1mm} \mathrm{sin}(x)\hspace{1mm} \mathrm{d}x$")
plt.show()


make_area_plot(areas=[simpson_areas, trapez_areas],
               correct_value=1,
               names_of_methods=["Simpson", "Trapez"],
               latex_title=r"$\int_0^{\pi/2} \hspace{1mm} \mathrm{sin}(x)\hspace{1mm} \mathrm{d}x$")
plt.show()
