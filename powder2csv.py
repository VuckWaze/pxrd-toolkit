def powder2csv(file):
    import pandas as pd
    import xrayutilities as xu    
    
    powder = xu.io.XRDMLFile(file)
    tt, I = powder.scan.scanmot, powder.scan.int
    DF = pd.DataFrame(columns=('2 Theta', 'Intensities'))
    DF['2 Theta'], DF['Intensities'] = tt, (I - min(I))/(max(I) - min(I))
    DF.to_csv(file[:-6] + '.csv', index=False)
    
    return