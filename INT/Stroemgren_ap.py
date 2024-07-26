import math
import numpy as np
import matplotlib.pyplot as plt

h = 6.6260755E-27
c = 2.99792458E+10
k = 1.380658E-16


def f(x):
   return math.sin(x)


def another_example(x):
    return np.exp(-x**2)

def yet_another_example(x):
    return 1/(1+x**2)


def yet_another_another_example(x):
    return x*np.exp(2*x)


def Trap(a, b, n, func):
    # N IS THE NUMBER OF INTERVALS
    h = (b - a) / n
    area = (func(a) + func(b))/2.0

    for i in range(1, n):
        x = a + i*h  # increment the x
        area = area + func(x)
    return area * h


def Romberging(a, b, N, func):

    first_value = Trap(a, b, 1, func)
    first_column = [first_value]

    for i in range(1, N+1):
        odd_numbers = [num for num in range(1, 2**i-1+2) if num%2!=0]  # flos contribution
        samples = [func(a+(b-a)/2**i*k) for k in odd_numbers]

        to_sum = np.array(samples)
        term = 1/2 * ((first_column[i-1]) + (b-a)/2**(i-1)*np.sum(to_sum))
        first_column.append(term)

    # alis approach
    romberg_matrix = np.zeros((N + 1, N + 1))
    romberg_matrix[:, 0] = first_column

    # alis approach
    for j in range(1, N + 1):
        for i in range(0, N + 1 - j):
            # is the row index
            # print("j, i ", j, i)
            # romberg_matrix[i, j] = area_contribution
            I = romberg_matrix
            area = (4 ** (j) * I[i + 1, j - 1] - I[i, j - 1]) / (4 ** (j) - 1)
            romberg_matrix[i, j] = area
        # j is the column index

    return romberg_matrix


def black_body_intensity(nu, T):
    # This is B_ν
    return (2*h*nu**3 / c**2) * 1/(np.exp(h*nu/(k*T))-1)


def black_body_luminosity(nu, T, R):
    # This is L_ν
    return 4*np.pi**2 * R**2 * black_body_intensity(nu, T)


def specific_luminosity_planck(nu):
    # Τhis is L_ν / (hν)
    stellar_temp = 40_000
    stellar_radius = 9*6.9599e10
    return black_body_luminosity(nu=nu, T=stellar_temp, R=stellar_radius) / (h*nu)

plot_integrand = False

if plot_integrand:
    nu_continuous = np.linspace(0, 2*1e16, 1000)
    plt.plot(nu_continuous, specific_luminosity_planck(nu=nu_continuous), "b.")
    plt.xlabel("Frequency")
    plt.ylabel("Integrand of $Q(H^0)$")
    plt.show()


lower_bound = 1.0973731569E+05 * c
upper_bound = 1.25e16*5
step_size = 15

integral_value_Q_H0_matrix = Romberging(a=lower_bound, b=upper_bound, N=step_size, func=specific_luminosity_planck)
best_romberg_value = integral_value_Q_H0_matrix[0, -1]

print("integral val: ", best_romberg_value)

def calculate_r_s(QH0, T, n_H):
    # stromgren radius
    return 3.15e-15 * (QH0/n_H**2)**(1/3) * T**(0.28)

strömgren_radius = calculate_r_s(QH0=best_romberg_value, T=0.75, n_H=10)

print("The strömgren radius is ", strömgren_radius)
