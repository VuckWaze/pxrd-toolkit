def pksgrab(file):
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    DF = pd.read_csv(file)
    x,y = DF['2 Theta'], DF['Intensities']
    plt.plot(x,y)
    plt.yticks(())
    plt.xlabel(r'$\mathrm{Angle (2\theta)}$')
    plt.title('Click peaks to find their positions.\nMiddle-click when done.')
    peaks = np.asarray(plt.ginput(n=0, timeout=0))
    xs, ys = peaks[:,0], peaks[:,1]
    dx = 0.01 * (max(x) - min(x))
    Y = np.zeros_like(ys)
    X = np.zeros_like(xs)

    for n in np.arange(0,len(xs)):
        low = xs[n] - dx
        high = xs[n] + dx
        Y[n] = max(y[(low < x) & (x < high)])
        X[n] = x[(y==Y[n])]
    output = np.asarray([X,Y])
    return output

