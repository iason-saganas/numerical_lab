from sys import stdin
import math as mt
import numpy as np


def f(x):
   return_val = mt.sin(x)
   return return_val

#def f(x):
#    return np.exp(-x**2)

#def f(x):
#    return 1/(1+x**2)


def Splines(dat):
    n1 = len(dat)
    n = n1 - 1
    X, Y = list(zip(*dat))
    X = [float(t) for t in X]
    Y = [float(t) for t in Y]
    a = Y[:]
    b = [0.0] * (n)
    d = [0.0] * (n)
    h = [X[i + 1] - X[i] for i in range(n)]
    al = [0.0] * n
    for i in range(1, n):
        al[i] = 3 / h[i] * (a[i + 1] - a[i]) - 3 / h[i - 1] * (a[i] - a[i - 1])
    c = [0.0] * n1
    v = [0.0] * n1
    u = [0.0] * n1
    z = [0.0] * n1
    v[0] = 1.0
    u[0] = z[0] = 0.0
    for i in range(1, n):
        v[i] = 2 * (X[i + 1] - X[i - 1]) - h[i - 1] * u[i - 1]
        u[i] = h[i] / v[i]
        z[i] = (al[i] - h[i - 1] * z[i - 1]) / v[i]
    v[n] = 1.0
    z[n] = c[n] = 0.0
    for j in range(n - 1, -1, -1):
        c[j] = z[j] - u[j] * c[j + 1]
        b[j] = (a[j + 1] - a[j]) / h[j] - (h[j] * (c[j + 1] + 2 * c[j])) / 3
        d[j] = (c[j + 1] - c[j]) / (3 * h[j])

    spl = []
    for i in range(n):
        spl.append((a[i], b[i], c[i], d[i]))
    return spl
    # End of Splines


def I_Spline(a, b, n):

    x = [0.0]*n
    y = [0.0]*n
    h = (b-a)/(n-1)
    for i in range(0, n):
        x[i] = a + i*h
        y[i] = f(x[i])
        # print("y", i, y[i])
    dat = list(zip(x,y))
    spl = Splines(dat)
    Y,A,B,C = list(zip(*spl))
    Y = [float(y) for y in Y]
    A = [float(y) for y in A]
    B = [float(y) for y in B]
    C = [float(y) for y in C]
    area = 0.0
    for i in range(0, n-1):
        h = x[i+1]-x[i]
        area = area + h*(Y[i]+h*(A[i]/2+h*(B[i]/3+h*C[i]/4)))
    return area
    # End of I_Spline


a = 0  # lower bound of integration interval
b = np.pi  # upper bound of integration interval
n = 11  # number of integration points
area = I_Spline(a, b, n)
print("With ", n-1, "splines, our estimate")
print("of the area from", a, "to", b, "= %.15f" % area)

errors = []
areas = []
ratios = []
for index, j in enumerate(range(1, 12)):
    n = 2**j
    print("n ", n)
    area = I_Spline(a, b, n)
    error = np.abs(2 - area)
    errors.append(error)
    areas.append(area)
    ratio = errors[index-1]/errors[index]
    ratios.append(ratio)

np.savetxt("tables_SPLINES/classic_table_splines_sin_function_integrated_upper_bound_pi.txt.txt", np.array([areas, errors, ratios]))

from simpsons_copy import SIMP

errors = []
areas = []
ratios = []
for index, j in enumerate(range(1, 12)):
    n = 2**j
    print("n ", n)
    area = SIMP(a, b, n, func=f)
    error = np.abs(2 - area)
    errors.append(error)
    areas.append(area)
    ratio = errors[index-1]/errors[index]
    ratios.append(ratio)

np.savetxt("tables_SPLINES/classic_sin_function_integrated_upper_bound_pi_with_SIMPS.txt", np.array([areas, errors, ratios]))


