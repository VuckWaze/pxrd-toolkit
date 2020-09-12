# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 10:55:35 2019

@author: Master
"""

def ads2fds(twoTheta, intensities, alpha=1):
    from numpy import sin, deg2rad
    a = deg2rad(alpha) # = 1
    R = 240 # (mm) Goniometer radius
    t = deg2rad(twoTheta/2)

    SL =  10

    A = R * sin(a/2)/SL
    B = 1/sin(t+a/2) + 1/sin(t-a/2)
    return intensities * A * B
