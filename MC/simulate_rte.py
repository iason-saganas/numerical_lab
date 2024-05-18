import numpy as np


def draw_tau():
    x = np.random.rand()
    return -np.log(x)


def draw_mu(mode: str):
    x = np.random.rand()
    if mode == "emit":
        return np.sqrt(x)
    elif mode == "scatter":
        return 2 * x - 1
    else:
        raise ValueError("Unknown mode")


def draw_projection(mode: str):
    return draw_tau(), draw_mu(mode)


def simulate(num_of_photons: int, tau_max: float):
    mu_store = []
    for photon in range(num_of_photons):
        # Emit the first photon
        delta_tau, delta_mu = draw_projection("emit")
        vertical_tau = tau_max - delta_tau * delta_mu
        while vertical_tau >= 0:
            # Photon is inside the atmosphere, but might have already been scattered back;
            # If that is the case, re-emit photon.
            # Else, scatter photon.
            # Append the result to `mu_store`.
            if vertical_tau > tau_max:
                delta_tau, delta_mu = draw_projection("emit")
                vertical_tau = tau_max - delta_tau * delta_mu
            else:
                # The Photon is actually inside the atmosphere => Scatter
                delta_tau, delta_mu = draw_projection("scatter")
                vertical_tau = vertical_tau - delta_tau * delta_mu
        mu_store.append(delta_mu)

    # Sanity-check:
    if num_of_photons != len(mu_store):
        raise ValueError("num phot vs len mu: \n"
                         f"{num_of_photons}, {len(mu_store)}")

    return np.array(mu_store)