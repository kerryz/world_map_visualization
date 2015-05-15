import matplotlib.pyplot as plt


def main():
    plot_markers()


def plot_markers():
    # define corner points
    x = 1
    y = 1

    fig, ax = plt.subplots()
    ax.plot(x, y, 'o', markerfacecolor="blue", markersize=100, markeredgecolor=None)
    ax.plot(100, 100, 'o', markerfacecolor="blue", markersize=100, markeredgecolor="black")

    plt.show()


if __name__ == "__main__":
    main()
