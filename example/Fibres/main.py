import os, sys
#~ import numpy as np
#~ import matplotlib.pyplot as plt

SOLVERGP = "/home/frocha/sources/Piola/build/Piola.x"
CONVERTER = "/home/frocha/sources/Piola/io/build/converter.x"

    
os.system("cp IniFile000.txt IniFile.txt")
os.system("cp Param000.txt Param.txt")

print('Simulation has begun')
os.system(SOLVERGP)
print('Simulation finished')

inputConv = " 3 '2 0 2 0 F F F' " # kind of converter, NtotalField , NParamCell , NvolGroups , NFace, isPeriodicNoLagrange , isXFEM ,  AddGlobalNode
inputConv += "'u 0 T T' 'l 2 T F' "# label , iShift , isPoint? , isVec?
 
os.system(CONVERTER + inputConv) 
 

