import numpy as np
import matplotlib.pyplot as plt

def click_peaks(ax=None, marker='v',color='black', write_x=True):
    if ax == None:
        ax = plt.gca()
    cds = ax.lines[len(plt.gca().lines) - 1] # cds = current data set
    x, y = cds.get_xdata(), cds.get_ydata()
    peaks = np.asarray(plt.ginput(n=0, timeout=0))
    xs, ys = peaks[:,0], peaks[:,1]
    # Make sure that the marker typ is a valid matplotlib marker type
    # https://matplotlib.org/3.1.1/api/markers_api.html

    # The following finds the index for the x-value lying closest
    #     to what was clicked:
    nearIdx = []
    for i in xs:
        nearIdx.append( ((x-i)**2).argmin() )

    plt.plot(x[nearIdx],
             y[nearIdx] + np.diff(plt.ylim())*0.02,
             marker=marker,
             linestyle='',
             color=color
             )
    if write_x:
        text = [str(np.around(X, 1))for X in x[nearIdx]]
        for n in range(len(text)):
    	    plt.text(x[nearIdx][n], y[nearIdx][n]+0.05*plt.ylim()[1],
    			text[n],
    			fontsize=9,
    			fontweight='bold',
    			horizontalalignment='right',
    			rotation=-45)
    plt.draw()
    return np.asarray((x[nearIdx], y[nearIdx]))

def click_free_peaks(marker='v',color='black', write_x=True):
    """
    Use this to draw markers anywhere in the current figure.

    Unlike click_peaks(), which finds the closest values for x and y
    in the current figure to the coordinates which were cicked, this
    function draws a marker where the user has clicked.

    """
    #plt
    peaks = np.asarray(plt.ginput(n=0, timeout=0))
    try:
        xs, ys = peaks[:,0], peaks[:,1]
    except:
        xs, ys = [],[]
    finally:
        pass
    # Make sure that the marker typ is a valid matplotlib marker type
    # https://matplotlib.org/3.1.1/api/markers_api.html

    # The following finds the index for the x-value lying closest
    #     to what was clicked:
    plt.plot(xs,
             ys + np.diff(plt.ylim())*0.02, #
             marker=marker,
             linestyle='',
             color=color
             )
    if write_x:
        text = [str(np.around(X, 1))for X in xs]
        for n in range(len(text)):
    	    plt.text(xs[n], ys[n]+0.05*plt.ylim()[1],
    			text[n],
    			fontsize=9,
    			fontweight='bold',
    			horizontalalignment='right',
    			rotation=-45)
    plt.draw()
    return np.asarray((xs, ys))

def single_free_peaks(marker='v',color='black', write_x=True):
    """
    Use this to draw markers anywhere in the current figure.
    Updates the figure with a marker everytime it is clicked.

    Unlike click_peaks(), which finds the closest values for x and y
    in the current figure to the coordinates which were cicked, this
    function draws a marker where the user has clicked.

    """
    #plt
    peaks = append(np.asarray(plt.ginput(n=1, timeout=0)))
    xs, ys = peaks[:,0], peaks[:,1]
    # Make sure that the marker typ is a valid matplotlib marker type
    # https://matplotlib.org/3.1.1/api/markers_api.html

    # The following finds the index for the x-value lying closest
    #     to what was clicked:
    plt.plot(xs,
             ys + np.diff(plt.ylim())*0.02, #
             marker=marker,
             linestyle='',
             color=color
             )
    if write_x:
        text = [str(np.around(X, 1))for X in xs]
        for n in range(len(text)):
    	    plt.text(xs[n], ys[n]+0.05*plt.ylim()[1],
    			text[n],
    			fontsize=9,
    			fontweight='bold',
    			horizontalalignment='right',
    			rotation=-45)
    plt.draw()
    return np.asarray((xs, ys))