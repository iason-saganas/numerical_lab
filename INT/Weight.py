# Compute the integral of f(x)g(x)dx from a to b using n intervals. # use arrays and trigonometric functions from numpy:
import numpy as np
import math as mt


def I1(x):
   # From 0 to np.pi/2 should return 1
   return np.sin(x)

def I2(x):
   # From 0 to np.pi should return 2
   return np.sin(x)


I1.correct_value = 1
I2.correct_value = 2


#def f(x):
#    return np.exp(-x**2)

# our function to be integrated, f(x)=1:
#def f(x):
#    return np.ones(len(x))

#def f(x):
#    return 1/(1+x**2)


# weights for weighting function g(x)=sin(x):
def w_sin(xl,xr,dx,xm):
    G = (np.cos(xl)-np.cos(xr))/dx
    H = ((np.sin(xr)-np.sin(xl))/dx*2.-(np.cos(xr)+np.cos(xl)))/dx
    return G, H

# weights for weighting function g(x)=1 (cf. trapezoid integration):
def w_1(xl,xr,dx,xm):
    G = np.ones(len(dx))
    H = np.zeros(len(dx))
    return G,H

def weights(xl,xr,g):
    dx =  xr-xl
    xm = (xr+xl)*.5
    G,H = g(xl,xr,dx,xm)
    Wm = dx*(G-H)*.5
    Wp = dx*(G+H)*.5
    return Wm,Wp

def integrate(a,b,n,f,g):
    x = np.array(list(range(n+1))) * (b-a)/n + a
    xl = x[0:n] # left sides of intervals
    xr = x[1:n+1] # right sides of intervals
    Wm,Wp = weights(xl,xr,g)
    return sum(f(xl)*Wm + f(xr)*Wp)


errors = []
areas = []
ratios = []
for index, j in enumerate(range(1, 12)):
    f = I2
    n = 2**j
    print("n ", n)
    area = integrate(0, np.pi, n, f, w_1)
    error = np.abs(f.correct_value - area)
    errors.append(error)
    areas.append(area)
    ratio = errors[index-1]/errors[index]
    ratios.append(ratio)

np.savetxt("tables_WEIGHT/I2_table_weight.txt", np.array([areas, errors, ratios]))



