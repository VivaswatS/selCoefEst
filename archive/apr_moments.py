## apragsdale script (Jan 7, 2022)
import moments
import numpy as np
import matplotlib.pylab as plt


def run_sim(n, gamma, a, N=5000):
    # initialize the spectrum, with 1 in singleton bin
    fs = moments.Spectrum(np.zeros(n + 1))
    fs[1] = 1
    # simulate a generations
    T = a / 2 / N
    # set relative size to 1, theta to 0 to forbid new mutations
    fs.integrate([1], T, gamma=gamma, theta=0)
    return fs


if __name__ == "__main__":
    n = 40  # 40 haploid samples
    a = 50  # 50 generations

    fig = plt.figure(figsize=(8, 4))
    ax1 = plt.subplot(1, 2, 1)

    for gamma in [-10, -50, -100, -500, -1000]:
        fs = run_sim(n, gamma, a)
        print("minimum value, with gamma =", gamma)
        print(fs.min())
        print(fs[1])
        ax1.plot(fs, ".--", label=gamma)
    ax1.legend()
    ax1.set_yscale("log")

    ax2 = plt.subplot(1, 2, 2)
    # run with larger sample size, then project down to n
    print("Larger sample sizes:")
    n_large = 200
    for gamma in [-10, -50, -100, -500, -1000]:
        fs = run_sim(n_large, gamma, a)
        fs = fs.project([n]) * n_large / n # adjust for sample size differences
        print("minimum value, with gamma =", gamma)
        print(fs.min())
        print(fs[1])
        ax2.plot(fs, ".--", label=gamma)
    ax2.legend()
    ax2.set_yscale("log")

    ax1.set_ylabel("Density")
    ax1.set_xlabel("Allele count")
    ax2.set_xlabel("Allele count")
    ax1.set_ylim(1e-60, 1)
    ax2.set_ylim(1e-60, 1)
    plt.tight_layout()
    plt.savefig("sample-size-test.pdf")
