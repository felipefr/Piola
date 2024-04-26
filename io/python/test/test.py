from Iconverter import *
import numpy as np

PK = homogenisationfibres("Param.txt",192,101,30,2,4,6,20).reshape((101,2,2),order = 'F')

print PK[0,:,:]
