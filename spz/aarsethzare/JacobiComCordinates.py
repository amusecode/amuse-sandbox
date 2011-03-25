from math import *
from numpy import *
from EulerAngles import EulerAngles
import sys

#Transform from Jacobi to centre of mass coordinates
def JacobiToCoMCoordinates(m, rin, vin, rout, vout) :

    m1 = m[0]
    m2 = m[1]
    m3 = m[2]
    m12 = m[0]+m[1]
    m123 = m12 + m[2]
    x = zeros((3,3))
    xdot = zeros((3,3))
    for i in range(len(rin)) :
        x[i,0]=-(m3/m123)*rout[i]-(m2/m12)*rin[i]
        xdot[i,0]=-(m3/m123)*vout[i]-(m2/m12)*vin[i]
        x[i,1]=-(m3/m123)*rout[i]+(m1/m12)*rin[i]
        xdot[i,1]=-(m3/m123)*vout[i]+(m1/m12)*vin[i]
        x[i,2]=(m12/m123)*rout[i]
        xdot[i,2]=(m12/m123)*vout[i]

    for i in range(len(rin)) :
        rin[i]=x[i,1]-x[i,0]
        vin[i]=xdot[i,1]-xdot[i,0]

        rout[i]=(m123/m12)*x[i,2]
        vout[i]=(m123/m12)*xdot[i,2]

    print "out:", rout[2], vout[2]
    # Transform from Jacobi to centre of mass coordinates
    pos = [[0.,0,0],[0,0,0],[0,0,0]]
    vel = [[0.,0,0],[0,0,0],[0,0,0]]
    for i in range(len(rin)) :
        #pos[i][0]=-(m3/m123)*rout[i]-(m2/m12)*rin[i]
        #vel[i][0]=-(m3/m123)*vout[i]-(m2/m12)*vin[i]
        #pos[i][1]=-(m3/m123)*rout[i]+(m1/m12)*rin[i]
        #vel[i][1]=-(m3/m123)*vout[i]+(m1/m12)*vin[i]
        #pos[i][2]= (m12/m123)*rout[i]
        #vel[i][2]= (m12/m123)*vout[i]
        pos[0][i]=-(m3/m123)*rout[i]-(m2/m12)*rin[i]
        vel[0][i]=-(m3/m123)*vout[i]-(m2/m12)*vin[i]
        pos[1][i]=-(m3/m123)*rout[i]+(m1/m12)*rin[i]
        vel[1][i]=-(m3/m123)*vout[i]+(m1/m12)*vin[i]
        pos[2][i]= (m12/m123)*rout[i]
        vel[2][i]= (m12/m123)*vout[i]

    return pos, vel

def extract_double_vector(line) :
    x = zeros(3)
    y = zeros(3)
    sline = line.split()
    for i in range(len(x)) :
        x[i] = float(sline[i])
    sline = line.split()
    for i in range(len(y)) :
        y[i] = float(sline[i+3])
    return x, y

def extract_single_vector(line) :
    x = zeros(3)
    sline = line.split()
    for i in range(len(x)) :
        x[i] = float(sline[i])
    return x

if __name__=="__main__":

    filename = 'RelPosVector.init'
    if len(sys.argv) > 1 :
        filename = sys.argv[1]
    fptr = open(filename)
    
    line = fptr.readline()
    m = extract_single_vector(line)
    line = fptr.readline()
    rin, vin = extract_double_vector(line)
    line = fptr.readline()
    rout, vout = extract_double_vector(line)
    pos, vel = JacobiToCoMCoordinates(m, rin, vin, rout, vout) 
    print "Jacobi to CoM coordinates: ", 
    print "Pos=", pos
    print "Vel=", vel

