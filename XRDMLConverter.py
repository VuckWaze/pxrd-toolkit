import pandas as pd
import os
import xrdtools as xt
import glob
print('Searching for XRDML files...')
files = glob.glob('/**/**/**/**/*.xrdml', recursive=True)
print('Found ' +  str(len(files)) + ' xrdml files.')
print('\nStarting conversion...')
c = 0
for file in files:
    fName = os.path.split(file)[1][:-6]
    fPath = os.path.split(file)[0]
    try:
        pxrd = xt.read_xrdml(fPath + '\\' + fName + '.xrdml')
        x = pxrd['2Theta']
        y = pxrd['data']
        y = (y-min(y))/(max(y)-min(y))
        DF = pd.DataFrame(data=[x, y])
        DF = DF.T
        DF.columns = ['2 Theta', 'Intensities']
        fList = fPath + '\\' + fName + '.csv'
        DF.to_csv(fList, index=False)
        c += 1
    except:
        print('Could not convert ' + fPath + '\\' + fName + '.')

print('\n===============================================\n\nConverted ' + str(c) + ' files.')
input()