import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation



global lines
lines = []
global data
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

def gen_line_plots(data, ax, fig, names=None, two_color=None, perspective=True):
    # Generate line plots
    X = np.linspace(-1, 1, data.shape[-1])
    for i in range(len(data)):
        # Small reduction of the X extents to get a cheap perspective effect
        if perspective is not None:
            xscale = 1 - i / 200.
            # Same for linewidth (thicker strokes on bottom)
            lw = 1.5 - i / 100.0
        else:
            xscale = 1
            lw = 1.5
        if two_color is not None:
            cmap = ListedColormap(['k', 'r'])
            norm = BoundaryNorm([-1, 0, 1], cmap.N)
            bound = two_color[i]
            line, = ax.plot(
                xscale * X[:bound + 1], i + data[i][:bound + 1], color="k", lw=lw)
            lines.append(line)
            line, = ax.plot(
                xscale * X[bound:], i + data[i][bound:], color="r", lw=lw)
            lines.append(line)
            # segments = np.concatenate(points_left, points_right, axis=1)
            # lc = LineCollection(segments, cmap=cmap, norm=norm)
            # lc.
            if (names is not None) and (len(names) > i):
                ax.text(1, i, names[i])
        else:
            line, = ax.plot(xscale * X, i + data[i], color="k", lw=lw)
            if (names is not None) and (len(names) > i):
                ax.text(1, i, names[i])
            lines.append(line)
    return lines


def update(*args):
    # Shift all data to the right
    data[:, 1:] = data[:, :-1]

    # Fill-in new values
    data[:, 0] = np.random.uniform(0, 1, len(data))

    # Update data
    for i in range(len(data)):
        lines[i].set_ydata(i + data[i])

    # Return modified artists
    return lines


def make_fig(data, output_filename='figs/pulsar.pdf', names=None,
             two_color=None, perspective=True):
    if len(data) == 0:
        return
    if type(data[0]) is type([]):
        data = np.array([np.array(x) for x in data])
    if type(data) is type([]):
        data = np.array([np.array(x) for x in data])
    print type(data)
    print type(data[0])
    print data.shape[:]
    globals()['data'] = data
    # Create new Figure with black background
    fig = plt.figure(figsize=(8, 8), facecolor='black')
    # Add a subplot with no frame
    ax = plt.subplot(111, frameon=False)
    gen_line_plots(
        data, ax, fig, names=names,
        two_color=two_color, perspective=perspective
    )
    # Set y limit (or first line is cropped because of thickness)
    ax.set_ylim(-1, len(data) + 5)
    # No ticks
    ax.set_xticks([])
    ax.set_yticks([])
    # Construct the animation, using the update function as the animation
    # director.
    if two_color is None:
        anim = animation.FuncAnimation(fig, update, interval=10)
    plt.savefig(output_filename, format='pdf')


if __name__ == '__main__':
    data = np.random.uniform(0, 1, (15, 75))
    make_fig(data, names=["***" + str(x) for x in range(len(data))],
             two_color=[4,2,29, 0, 75] + [5] * 10)

