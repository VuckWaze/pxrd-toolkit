# -*- coding: utf-8 -*-
"""
Created on Fri May  3 14:58:04 2019

@author: mita3616
"""

def tt2d(tt):
    from numpy import sin, radians
    CuKa = 1.5405980e-10
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

        t = radians(t)
    d = CuKa/(2*sin(t))
    return(d)
