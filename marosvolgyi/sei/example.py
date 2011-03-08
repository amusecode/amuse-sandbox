import os
import sys
import numpy as np
import math

from interface import SeiInterface

OMEGA = 1.0


if __name__ == '__main__':
    instance = SeiInterface()
    instance.initialization()
    impactparameter = 8.0 * 0.69336127435063
    instance.set_state(impactparameter, 3.1415926535 * 300.0/2.0*OMEGA*impactparameter,0,0,-3.0/2.0*OMEGA*impactparameter,0)
    for t in np.arange(0, 600,1):
        instance.evolve(t)
        print instance.get_state()['x'],instance.get_state()['y'], instance.get_state()['z'] 

    instance.stop()
    
