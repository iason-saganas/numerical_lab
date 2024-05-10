import numpy as np
import matplotlib.pyplot as plt
from styles.matplotlib_style import *


def euler_method(dydt: callable, t0: float, y0: float, t_final: float, step_size: float) -> np.array:
    """
    Euler's approximation method for differential equation.
    This function implements

        y_n+1 = y_n + step_size * dydt(t_n, y_n, step_size).

    :param dydt: callable,        The derivative of y(t).
    :param t0:  float,            The initial time.
    :param y0:  float,            The initial y value at t0.
    :param t_final:               The final time to evolve to.
    :param step_size:             The step size.
    :return: results, tuple       The found approximate solutions at given times
    """
    results = [y0]
    times = [t0]
    # y_1 = y0 + ( t0 + step_size ) * dydt(y0)
    idx = -1
    t = t0
    while t <= t_final:
        idx += 1  # Increment counting index to grab old y value
        y_old = results[idx]  # Grab old y value
        y_update = y_old + (t0 + step_size) * dydt(y_old)  # Update y value according to Euler scheme
        results.append(y_update)  # Append the results
        t = t + step_size  # Update the time
        times.append(t)  # Append the new time to return at the end

    return np.array(times), np.array(results)


def y_der(y_vals):
    return y_vals


def real_sol(x_vals):
    return np.exp(x_vals)


x = np.linspace(0, 5, 100)

for h, col in zip([0.1, 0.001, 0.0001], ["r", "g", "b"]):
    times, y_sol = euler_method(dydt=y_der, t0=0, y0=1, t_final=5, step_size=h)
    plt.plot(times, y_sol, ".", color=col, label=f"Euler Method, step-size {h}")

# plt.plot(x, real_sol(x), "b-", label="Real solution")
plt.xlabel("$x$ values")
plt.ylabel("$y$ values")
plt.legend()
plt.show()


