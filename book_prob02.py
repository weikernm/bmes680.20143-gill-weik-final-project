import sys

import numpy as np

import matplotlib as mlib
mlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from cobweb import *

def func_finite_a(x, r=0.5):
    return r * x * (1 - x)
def func_finite_b(x):
    return -(x*x) * (1 - x)
def func_finite_c(x):
    return 1./(2 + x)
def func_finite_d(x):
    return x * np.log(x**2)

fig = Figure(figsize=(6,6))
canvas = FigureCanvas(fig)
ax_a = fig.add_subplot(2,2,1)
cobweb_plot(ax_a, func_finite_a, 0.05, np.linspace(-0.1, 0.1, 500),
        n_iter=50)
ax_a.set_xlim(-0.1, .10)
ax_a.set_title('a: $x_{n+1} = 0.5 x_n (1 - x_n)$')
ax_b = fig.add_subplot(2,2,2)
cobweb_plot(ax_b, func_finite_b, 1.62, np.linspace(1, 2, 500),
        n_iter=5)
ax_b.set_xlim(1.4, 2)
ax_b.set_title('b: $x_{n+1} = -x_n^2 (1-x_n)$')
ax_c = fig.add_subplot(2,2,3)
cobweb_plot(ax_c, func_finite_c, 0.7, np.linspace(0.2, 0.75, 500),
        n_iter=50)
ax_c.set_title('c: $x_{n+1} = 1/(2+x_n)$')
ax_d = fig.add_subplot(2,2,4)
cobweb_plot(ax_d, func_finite_d, 2, np.linspace(0.1, 4, 500),
        n_iter=4)
ax_d.set_xlim(1, 4)
ax_d.set_title('d: $x_{n+1} = x_n \mathrm{ln}x_n^2$')
for ax in [ax_a, ax_b, ax_c, ax_d]:
    ax.grid()
    ax.set_xlabel('$x_n$')
    ax.set_ylabel('$x_{n+1}$')

fig.tight_layout()

canvas.print_figure(sys.argv[1])
