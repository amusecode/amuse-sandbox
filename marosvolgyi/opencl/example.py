from amuse.community import *

from interface import OpenclInterface
import numpy
import time

if __name__ == '__main__':
    instance = OpenclInterface(redirection='none')
    print instance.initialization()
    s =raw_input()
    instance.stop()
