import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from normalize import normalize
import numpy as np
from numpy.polynomial import polynomial as P
from tt2d import tt2d
from smooth import smooth
from scipy import interpolate

plt.style.use('seaborn-white')
mpl.use('Qt5Agg')
fPath = r'C:/Users/Master/OneDrive - Stockholm University/Data/SVOH/'


def loadpxrd_data(fPath):
    files = os.listdir(fPath)
    files = [f for f in files if f[-3:] == 'csv']
    names = [f[:-3] for f in files]
    idx = np.arange(len(names))
    DF = {}
    for n in idx:
        DF[names[n]] = pd.read_csv(fPath + files[n])
    return DF

def getxy(df):
    x,y = df['2 Theta'].values, df['Intensities'].values
    return x,y

def plotTextPowder():
    plt.xlabel(r'$\mathrm{Angle (2\theta)}$')
    plt.minorticks_on()
    plt.tick_params(axis='x', length=4)
    plt.tick_params(axis='x', which='minor', length=2)
    plt.ylabel(r'$\mathrm{Counts (a.u.)}$')
    plt.yticks(())
    return

DF = loadpxrd_data(fPath)

# Individual plots
# for key in DF.keys():
#     plt.cla()
#     x,y = getxy(DF[key])
#     y = smooth(y, 5)
#     plt.plot(x, y,
#              c = [0.1, 0.25, 0]
#              )
#     plotTextPowder()
#     plt.savefig(fPath + key)

x, y = {}, {}
keys = ('1600',
        '1660',
        '3200',
        '3260',
        '3220'
        )
titles = dict(zip(keys, ('1.6H, 1 bar Ar',
                         '1.6H, 60 bar $\mathrm{H_2}$',
                         '3.2H, 1 bar Ar',
                         '3.2H, 60 bar $\mathrm{H_2}$',
                         '3.2D, 20 bar $\mathrm{D_2}$')))
labs = ('48h','96h','Washed')

colors = (
    (0.1, 0., 0.1),
    (0.6, 0.2, 0.6),
    (0., 0.7, 1)
          )
sp = 1

for k in keys:
    plt.figure(figsize=(3,5))
    i = 0
    e = 1
    lab = []
    for key in DF.keys():
        if key[4:8] == k:
            if key[-1] == 'W':
                lab = 'Post wash'
            else:
                lab = labs[i]
                i += 1
            x,y = getxy(DF[key])
            y = smooth(y,size=np.round(0.0025*len(y)))
            plt.plot(x, y + 0.5*e,
                     c = colors[i-1],
                     label = lab
                      )
            plt.xlim((31.4,34.4))
            plt.title(titles[k])
            e += 1
    plotTextPowder()
    plt.legend()
    plt.savefig(fPath + k + '110-101')
    #plt.title(titles[sp-1])
    sp += 1

plt.show()
