# -*- coding: utf-8 -*-
"""
Created on Fri May  3 14:21:26 2019

@author: mita3616
"""

def readXRDML(file, normalize=True):
    import pandas as pd
    import xrayutilities as xu    
    
    powder = xu.io.XRDMLFile(file)
    tt, I = powder.scan.scanmot, powder.scan.int
    DF = pd.DataFrame(columns=('2 Theta', 'Intensities'))
    if normalize:
        DF['2 Theta'], DF['Intensities'] = tt, (I - min(I))/(max(I) - min(I))
    else:
        DF['2 Theta'], DF['Intensities'] = tt, I
    
    return(DF)