# -*- coding: utf-8 -*-
"""
A module created for working with powder x-ray diffraction data.

The most useful content thus far is the pxrdPattern() class which can take
 xrdml-files (or .csv) and has built-in tools and functions for data analysis
 and processing.

See doc-strings for info about the individual components

@author: Master
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import xrayutilities as xu  # Tools for reading xrdml-files
# from crystals import Crystal, Element # Class for crystal structures, https://crystals.readthedocs.io
import scipy.constants as __physical_constants__


# Default directories
__data_dir__ = r'%USERPROFILE%\Documents\PXRD'
__cif_dir__= __data_dir__ + '\CIF-files'
__pks_dir__ = __data_dir__ + '\PXRD Peak lists'
__pxrd_dir__ = __data_dir__ + '\PXRD'

# Other

CuKa = (1.5406e-10, 'm') #wavelength

"""
_elements_ = {} # Table of the elements
for a in Element.valid_symbols:
    atom = Element(a)
    _elements_[a] = (atom.element_full,
                     atom.atomic_number,
                     atom.mass,
                     atom.magnetic_moment_ground)

"""

# Atomic form factors
#
#_aff_ = pd.read_csv(__data_dir__ + '\atomic-form-factors-xrays.txt')


##############################################################################
# Utilities
#

def powder2csv(file, scale_intensities=True):
    """
    Convert an xrdml-scan to a csv-file using pandas built-in converter.
      Scales intensities by default according to:
          Y = (Y - min(Y))/(max(Y) - min(Y))
    """
    powder = xu.io.XRDMLFile(file)                                      # Read xrdml-file
    tt, I = powder.scan.scanmot, powder.scan.int                        # Get 2Θ and counts from file
    DF = pd.DataFrame(columns=('2 Theta', 'Intensities'))               # Create an empty pandas dataframe and set column namnes
    if scale_intensities: I = (I - min(I))/(max(I) - min(I))             # If scaling was requested, scale intensities
    DF['2 Theta'], DF['Intensities'] = tt, I                            # Append data to dataframe
    DF.to_csv(file[:-6] + '.csv', index=False)                          # Save as csv
    return

def tt2d(tt, wavelength='CuKa1'):
    """
    Convert a 2theta value to d-spacings from Braggs law:

        wl = 2d*sin(Θ)
        d = wl/(2*sin(Θ))

    """
    CuKa = xu.wavelength('CuKa1')# 1.5405980 Å
    #BRAGGS LAW:
    # wl = 2d sin(theta)
    if type(tt) == 'float':
        t = tt/2

    else:
        try:
            t = float(tt)/2
        except ValueError:
            print('Could not convert input value to float type.')
            print(f'Input value: {tt}\nType: {type(tt)}')

        t = np.radians(t)
    d = CuKa/(2*np.sin(t))
    return(d)

def d2tt(d, wavelength='CuKa1'):
    """
    Convert  d-spacings to 2Theta:

        wl = 2d*sin(Θ)
        Θ = arcsin(wl/2d)

    wl: Wavelength of radiation source.
    """
    CuKa = xu.wavelength(wavelength)# 1.5405980 Å
    #BRAGGS LAW:
    # wl = 2d sin(theta)
    t = np.arcsin(CuKa/(2*d)) # t = theta
    return(2*t)



##############################################################################
##############################################################################
# Signal treatment functions


##############################################################################
##############################################################################
def calculate_pattern():
    """
    To calculate the intensity of a reflection:

    Structure factor equation:
        Sum of atomic scattering factors times wave functions for planes
        with indicies [h k l] and for positions x,y,z:

            ____
            \
    F(s) =   >    f_n(ELEMENT) × exp(2πi[hx + ky + lz])
            /___
              n

    To determine 2Θ positions, calculate reflections
        """

    return

# Plot x & y with pxrd formatting. For use in the pxrdPattern class, mainly.
def __plot__(pattern):
    x,y = pattern.tt, pattern.I
    plt.plot(x,y)
    plt.xlabel('Angle (2θ)')
    plt.ylabel('')
    plt.yticks(())
    plt.minorticks_on()
    return


# Used to format an existing plot with PXRD relevant axes.
def format_plot():
    ax = plt.gca()
    ax.set_xlabel('Angle (2θ)')
    ax.set_ylabel('')
    ax.set_yticks(())

##############################################################################
##############################################################################
##############################################################################

class pxrdPattern():
    """
    A class for working with experimental PXRD patterns in an object-oriented
      manner so that routine analysis and processing can be performed much
      more rapidly.

    My intention is to create built-in tools for processing, visualizing and
      accessing data by simple means.
    """

    def __init__(self, file=None):
        # 2 Theta and intensity data is enough for now.
        self.data= []
        self.__scan__ = []
        self.properties = {}
        self.tt = []
        self.I = []
        self.name = ''
        if type(file) != type(None):
            self.load(file)
        self.methods = __methods__(self)

    def load(self, file):
        """
        Use to load data into a pxrdPattern object.
        """
        self.name, fileExt = os.path.splitext(file)
        if fileExt == '.xrdml':
            xrdml = xu.io.XRDMLFile(file)
            self.__scan__ = xrdml.scan
            self.tt, self.I = xrdml.scan.scanmot, xrdml.scan.int
            self.I = (self.I- self.I.min()) / (self.I.max() - self.I.min()) # Scale intensities
        elif fileExt == '.csv':
            self.__scan__ = np.array([])
            df = pd.read_csv(file)
            self.tt = df['2 Theta'].values
            self.I = df['Intensities'].values

    def add_indices(self, hkl, tt):
        self.hkl =  pd.DataFrame(np.asarray([tt, hkl]),
                                 columns=('2 Theta', 'hkl'))

    def plot(self):
        __methods__.plot()
        self.plot = __methods__.plot()
        self.plot()
        format_plot()

    def draw_pattern(self):
        self.methods.draw_pattern()


#class CrystalStructure():
#    """
#    Supposed to be simpler alternative of the Crystal() class.
#    """
#    __args__ = ['name',
#                'parameters', # [a, b, c, α, ß, g]
#                'space_group', # {}
#                'atoms' # np.array[ATOM_TYPE, X, Y, Z]
#                ]
#    def __init__(self, *args):
#        self.assign
#        assert(len(args) < =)
#        assert(type(arg[0]) == str)


class __methods__():
    """
    Class for initializing powder methods.
    Use this for inheritance
    """

    def __init__(self, data):
        self.tt = data.tt
        self.I = data.I
        self.__name__ = data.name

        self.__hasAtoms__ = False
        self.__has__ = False

    def draw_pattern(self):
        plt.figure()
        x,y = self.tt, self.I
        self.line = __plot__(self)
        format_plot()

    def find_pks(self):
        from scipy.signal import find_peaks as fp
        x,y = self.tt, self.I
        peaks, _ = fp(y, # Intensities
                      distance=round(0.025*len(x)) # N points between each peak.
                      )

    def draw_pks(self):
        try :
            p = self.line.properties()
        except :
            self.line = plt.Line2D(self.tt, self.I)
            p = self.line.properties()
        from scipy.signal import find_peaks as fp
        x,y = self.tt, self.I
        peaks, _ = fp(y, distance=0.025)
        self.pks_line = plt.plot(x,peaks,lw=0, color=p['color'], marker='|')
    # Add Lines2D-object into current axes. Create waterfall plot by default.
    def plot(self):
        stack = True
        if stack ==True & len(plt.gca().lines) != 0:
            __yData__ = np.array([0])
            __diff__ = 0
            for l in plt.gca().lines:
                l = l.get_ydata()
                if l.max() > __yData__.max():
                    __yData__ = l
            __diff__ = __yData__.max() - __yData__.min()
            __offset__ = __yData__.max()
        else:
            __offset__ = 0
        x,y = self.tt, self.I
        self.line = plt.plot(x,y+__offset__,label=self.__name__)[0]
        try : del __offset__, __yData__, __diff__, l
        except : pass

class reflectionTable():
    """
    Loads a reference table of reflections as an object.
    """
    def __init__(self, name=None):
        self.table_data = pd.DataFrame()
        if type(name) == str:
            self.table_data = read_reflections(name) # defined below
            self.name = name
        elif type(name) == type(None):
            self.table_data = pd.DataFrame(columns=read_reflections('SrVO3').columns)
        else:
            try: self.table_data = name.table_data
            except: TypeError(f'Could not read {name} as a reflection table.')
        self.I = self.table_data['I']
        self.tt = self.table_data['2Theta']
        self.values = self.table_data.values

    # Return table without K_a2
    def strip_k_alpha2(self):
        new = reflectionTable()
        new.table_data = self.table_data.loc[self.table_data['ID(alpha)'] == 1].copy()
        #print(new.table_data)
        return new
    # Only show table values where 2 theta is between "low" and "high"
    def tt_filter(self, low, high):
        new = reflectionTable()
        new.table_data = self.table_data.loc[(low < self.table_data['2Theta']) & (self.table_data['2Theta'] < high)].copy()
        #print(new.table_data)
        return new

    # Only return 2 theta within range of axes x-lim in current plot
    def in_xlim(self):
        new = reflectionTable()
        low, high = plt.xlim()
        new.table_data = self.table_data.loc[(low < self.table_data['2Theta']) & (self.table_data['2Theta'] < high)].copy()
        #print(new.table_data)
        return new

    # Return table of hkl
    def hkl(self):
        new = reflectionTable()
        table = self.table_data[['h', 'k', 'l']].copy()
        table = table['h']
        #print(new.table_data)
        return new

    # Plot the index of a powder pattern in an axes
    # Uses the last line in the axes object

##############################################################################
##############################################################################
##############################################################################
# Reflection data


def read_reflections(structure):
    """Fetch exported VESTA powder pattern for a structure"""

    path = r'C:\Users\Master\OneDrive - Stockholm University\Data\PXRD Peak lists' # Default path
    try:
        pks =  pd.read_csv(path + '/' + structure + '.rflx')
    except:
        raise FileNotFoundError(f'Could not find "{structure}" file.')
    finally:
        pass

    try:
        pks = pks.apply(pd.to_numeric)
    except ValueError: "Could not convert all values to numeric"
    return pks
# Function for drawing lines in the pattern

def draw_reflection_lines(structure, intensity_threshold=5, c='blue', ls=':'):
    # Set path
    path = __pks_dir__
    pks = read_reflections(structure)
    mapPks = pks.loc[(pks['I'] > float(intensity_threshold)) & (pks['ID(alpha)']==1)] # mapped peaks
    tt, I = mapPks['2Theta'].copy(), mapPks['I'].copy()
    for p in tt.values:
        a = 0.5# I.loc[tt==p] # Transparency of line
        plt.axvline(p,
                   alpha=a,
                   c=c,
                   ls=ls,
                   lw=0.75,
                   label=f'{structure}_{p}'
                   )
    y = plt.ylim()[1]
    y = y + 0.1*y
    x = tt.min()
    prettyName = ''
    for s in structure:
        if s.isdigit():
            s = '_' + s
        prettyName = prettyName + s
    plt.text(x,y,
             f'${prettyName}$',
             horizontalalignment='left',
             rotation=45,
             color=c)
    return mapPks,

def phase_markers(structure, intensity_threshold=5, c='blue', ls=':'):
    pks = read_reflections(structure)
    mapPks = pks.loc[(pks['I'] > float(intensity_threshold)) & (pks['ID(alpha)']==1)] # mapped peaks
    tt, I = mapPks['2Theta'].copy(), mapPks['I'].copy()
    for p in tt.values:
        a = 0.5# I.loc[tt==p] # Transparency of line
        plt.axvline(p, alpha=a, c=c, ls=ls, lw=0.75)
    y = plt.ylim()[1]
    y = y + 0.1*y
    x = tt.min()
    prettyName = ''
    for s in structure:
        if s.isdigit():
            s = '_' + s
        prettyName = prettyName + s
    plt.text(x,y,
             f'${prettyName}$',
             horizontalalignment='left',
             rotation=45,
             color=c)
    return mapPks
