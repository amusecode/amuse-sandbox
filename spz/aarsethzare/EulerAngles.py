from math import *
from numpy import *
import sys

# Rotates vector given wrt (i,j,k) through Euler angles (i,w,Omega).
def EulerAngles(angles, x) :

    cosi=cos(angles[0])
    sini=sin(angles[0])
    cosw=cos(angles[1])
    sinw=sin(angles[1])
    cosOm=cos(angles[2])
    sinOm=sin(angles[2])

    B = matrix(
        [
         [cosw*cosOm-sinw*cosi*sinOm, cosw*sinOm+sinw*cosi*cosOm, sinw*sini], 
         [-sinw*cosOm-cosw*cosi*sinOm, -sinw*sinOm+cosw*cosi*cosOm, cosw*sini],
         [sini*sinOm, -sini*cosOm, cosi]
        ])

#        B(1,1) = cosw*cosOm - sinw*cosi*sinOm
#        B(1,2) = cosw*sinOm + sinw*cosi*cosOm
#        B(1,3) = sinw*sini
#        B(2,1) = -sinw*cosOm - cosw*cosi*sinOm
#        B(2,2) = -sinw*sinOm + cosw*cosi*cosOm
#        B(2,3) = cosw*sini
#        B(3,1) = sini*sinOm
#        B(3,2) = -sini*cosOm
#        B(3,3) = cosi

    temp = []
    for i in range(len(x)) :
        sum=0
        for j in range(len(x)) :
            sum += B[j,i]*x[j]
        temp.append(sum)
    y = array([temp[0], temp[1], temp[2]])
#    for ti in temp :
#        y.append(ti)
    return y

if __name__=="__main__":

    angles = []
    x = []
    if len(sys.argv) > 6 :
        angles.append(radians(float(sys.argv[1])))
        angles.append(radians(float(sys.argv[2])))
        angles.append(radians(float(sys.argv[3])))
        x.append(float(sys.argv[4]))
        x.append(float(sys.argv[5]))
        x.append(float(sys.argv[6]))
        print "Euler Angles: ", EulerAngles(angles, x) 
    else :
        print "give angles and vector on the command line."
        print "   The angles are: inclination"
        print "                   argument of periastron"
        print "           and the longitude of line of nodes (of inner binary)."
        print "   The vector is a position vector or velocity vector"
        print "                                 (in arbitrary units)."


