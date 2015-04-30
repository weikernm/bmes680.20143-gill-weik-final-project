import numpy as np

def iter_finite(func_finite, x_0, n=100):
    '''
    Iterate a finite difference equation given a starting point `x_0'
    (can be a vector), and a function which accepts input of the form of
    `x_0'.
    The equation will be iterated `n' times (default: 100).
    '''
    x_i = np.empty(np.shape(x_0) + (n,))
    x_i[..., 0] = x_0
    for t in np.arange(1, n):
        x_i[..., t] = func_finite(x_i[..., t-1])
        pass
    return x_i

def cobweb_plot(ax, func_finite, x_0, x_n, n_iter=100,
        func_kwarg=dict()):
    '''
    Create a cobweb plot on the supplied axes `ax', of the finite
    difference function `func_finite', given an initial value `x_0',
    an array of values over which to compute which value would be next
    `x_n', and a number of times to iterate the finite difference
    function `n_iter'.
    '''
    x_i = iter_finite(func_finite, x_0, n=n_iter, **func_kwarg)
    x_n1 = func_finite(x_n)
    ax.plot(x_n, x_n)
    ax.plot(x_n, x_n1)
    x_cobweb = np.repeat(x_i[:-1], 2)
    y_cobweb = np.hstack((x_i[0], np.repeat(x_i[1:-1], 2), x_i[-1]))
    ax.plot(x_cobweb, y_cobweb)

