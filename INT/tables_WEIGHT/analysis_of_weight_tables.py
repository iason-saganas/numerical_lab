import matplotlib.pyplot as plt
import numpy as np
from styles.matplotlib_style import *

def make_error_plot(errors, ratios, names_of_methods, latex_title, double=True):
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

    if double:


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

    else:
        for i, error in enumerate(errors):
            plt.plot(x, error, ".-", label=names_of_methods[i])
        plt.xlabel("Step $n$")
        plt.ylabel(r"Error $\varepsilon$")
        plt.yscale("log")
        plt.legend()
        plt.title(latex_title)
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
    plt.suptitle(latex_title, fontsize=20)
    plt.legend()
    plt.tight_layout()


weight_I1_areas, weight_I1_errors, weight_I1_ratios = np.loadtxt("I1_table_weight.txt")
weight_I2_areas, weight_I2_errors, weight_I2_ratios = np.loadtxt("I2_table_weight.txt")

spline_I1_areas, spline_I1_errors, spline_I1_ratios = np.loadtxt("/Users/iason/PycharmProjects2024/NumLab/INT/unequal steps/splines/tables_SPLINES/classic_table_splines_sin_function_integrated_upper_bound_pi_over_2.txt")
spline_I2_areas, spline_I2_errors, spline_I2_ratios = np.loadtxt("/Users/iason/PycharmProjects2024/NumLab/INT/unequal steps/splines/tables_SPLINES/classic_table_splines_sin_function_integrated_upper_bound_pi.txt.txt")

make_error_plot(errors = [weight_I1_errors, spline_I1_errors],
                ratios = [weight_I1_ratios, spline_I1_ratios],
                names_of_methods=["Constant weights", "Spline"],
                latex_title=r"$\int_0^{\pi/2} \hspace{1mm} \mathrm{sin}(x)\hspace{1mm} \mathrm{d}x$")
plt.show()


make_error_plot(errors = [weight_I2_errors, spline_I2_errors],
                ratios = [weight_I2_ratios, spline_I2_ratios],
                names_of_methods=["Constant weights", "Spline"],
                latex_title=r"$\int_0^{\pi} \hspace{1mm} \mathrm{sin}(x)\hspace{1mm} \mathrm{d}x$")
plt.savefig("Error comparison spline and constant weights sine with upper bound pi")
plt.show()


plt.subplot(2, 1, 1)

make_area_plot(areas = [weight_I1_areas, spline_I1_areas],
               correct_value=1,
                names_of_methods=["Constant weights", "Spline"],
                latex_title=r"$\int_0^{\pi/2} \hspace{1mm} \mathrm{sin}(x)\hspace{1mm} \mathrm{d}x$")

plt.subplot(2, 1, 2)

make_area_plot(areas = [weight_I2_areas, spline_I2_areas],
               correct_value=2,
                names_of_methods=["Constant weights", "Spline"],
                latex_title=r"$\int_0^{\pi} \hspace{1mm} \mathrm{sin}(x)\hspace{1mm} \mathrm{d}x$ (lower) and $\int_0^{\pi/2} \hspace{1mm} \mathrm{sin}(x)\hspace{1mm} \mathrm{d}x$ (upper)")
plt.savefig("Comparison of areas spline and constant weight sine.png")
plt.show()



