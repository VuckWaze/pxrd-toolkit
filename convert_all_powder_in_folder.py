import os
from powder2csv import powder2csv

files = [f for f in os.listdir() if f[-5:] == 'xrdml']
for f in files:
	powder2csv(f)
	