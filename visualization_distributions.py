import numpy as np
from styles.matplotlib_style import *

show_mu_vs_cosine_plot = True
show_distribution_at_lower_boundary = True
show_distribution_inside_atmosphere = True

if show_distribution_at_lower_boundary:
    plt.subplot(2, 1, 1)

    plt.title(r"Distribution of $\mu$ and $\vartheta$ at lower boundary")

    x = np.random.rand(int(10000))
    mu = x**(1/2)
    theta = np.arccos(mu)

    plt.hist(mu,  histtype="step", color=blue, lw=2, label=r"$\mu \hookleftarrow P(\mu)\propto \mu$", bins="auto")
    plt.legend()

    plt.subplot(2, 1, 2)

    plt.hist(theta*(360/(2*np.pi)),  histtype="step", color=red, lw=2,
             label=r"Resulting distribution in $\vartheta[^{\circ}]$", bins="auto")

    plt.xlabel("Value of random variable")
    plt.ylabel("Frequency")
    plt.legend()

    plt.tight_layout()
    plt.savefig("figures/2_mu_distribution_at_lower_boundary.png")

    plt.show()

if show_mu_vs_cosine_plot:

    theta_continuous = np.linspace(0, 2*np.pi, 300)
    mu_continuous = np.linspace(-1, 1, 300)
    plt.plot(mu_continuous, np.arccos(mu_continuous)*360/(2*np.pi), color="purple")
    plt.xlabel(r"Value of $\mu$")
    plt.ylabel(r"$\vartheta = \mathrm{arccos}(\mu)$ in $^{\circ}$")
    plt.savefig("figures/1_mu_vs_cosine.png")
    plt.show()

if show_distribution_inside_atmosphere:
    plt.subplot(2, 1, 1)

    plt.title(r"Distribution of $\mu$ and $\vartheta$ inside atmosphere")

    x = np.random.rand(int(10000))
    mu = 2*x-1
    theta = np.arccos(mu)

    plt.hist(mu,  histtype="step", color=blue, lw=2, label=r"$\mu \hookleftarrow P(\mu)= 1/2$", bins="auto")
    plt.legend()

    plt.subplot(2, 1, 2)

    plt.hist(theta*(360/(2*np.pi)),  histtype="step", color=red, lw=2,
             label=r"Resulting distribution in $\vartheta[^{\circ}]$", bins="auto")

    plt.xlabel("Value of random variable")
    plt.ylabel("Frequency")
    plt.legend()

    plt.tight_layout()
    plt.savefig("figures/3_distribution_of_mu_inside_atmosphere.png")

    plt.show()
