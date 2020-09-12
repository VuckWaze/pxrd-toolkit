def d2tt(d):
    from numpy import arcsin, radians, degrees
    CuKa = 1.5405980e-10
    #BRAGGS LAW:
    # wl = 2d sin(theta)
    t = degrees(arcsin( CuKa / (2*d) ) )
    tt = 2 * t
    return(tt)
