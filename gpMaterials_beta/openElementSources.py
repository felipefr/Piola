import os, sys
import numpy as np

SOLVERGPPsrc = "/home/felipe/PToolV201311/SolverGPP2016/src/"

f = open("Basparam.txt")

found = False
elements = []
libs = []
for line in f:
	
	if found:
		r = line.split(' ')
		for e in r:
			v = int(e)
			if(v>0):
				elements = elements + [v]
			else:
				libs = libs + [v]
		
		break
		 
	elif "*Used Elements" in line:
		found = True
		
d = {}
d.update({100 : 'others/NullElement'})
d.update({101 : 'others/PolyElement'})
d.update({502 : 'Multiscale/enforcePeriodic2D'})
d.update({506 : 'Multiscale/computeTotalDisp2D'})
d.update({508 : 'others/PostProcessingElem2D'})
d.update({511 : 'Multiscale/computeTangHom'})
d.update({512 : 'Multiscale/minRestrictionBC2DExact'})
d.update({519 : 'others/DirichletNode'})
d.update({520 : 'damage/DamageGeneric'})
d.update({522 : 'FiniteStrain/finiteStrainNew'}) 
d.update({523 : 'Multiscale/flucTangHom'})
d.update({524 : 'others/PostProcessingRVEboundary'})
d.update({525 : 'FiniteStrain/finiteStrainDelta'})
d.update({529 : 'FiniteStrain/finiteStrainGen'})
d.update({530 : 'FiniteStrain/FbarReduced'})
d.update({531 : 'FiniteStrain/FiniteStrainSpatial'})
d.update({532 : 'FiniteStrain/FbarMaterial'})
d.update({533 : 'FiniteStrain/NeumannFiniteStrain'})
d.update({534 : 'FiniteStrain/NeumannFiniteStrainSpatial'})
d.update({538 : 'FiniteStrain/FbarSpatial'})
d.update({539 : 'others/NodalForce'})
d.update({540 : 'others/computePressure'})
d.update({541 : 'Multiscale/flucTangHomGen'})
d.update({542 : 'Multiscale/computeTangHomGen'})
d.update({600 : 'Multiscale/DecisionBifurcation'})
d.update({601 : 'Multiscale/LocDomainBC'})
d.update({700 : 'Multiscale/computeMinDetQ'})
d.update({701 : 'Multiscale/MarkLocPoint'})
d.update({702 : 'Multiscale/LocPointsBC'})
d.update({800 : 'FiniteStrain/FSgen'})
d.update({801 : 'FiniteStrain/FSFbar'}) 
d.update({802 : 'Multiscale/TotalDisp'})
d.update({803 : 'damage/DamageGenericGen'})
d.update({804 : 'Multiscale/canonicalProblem'})
d.update({805 : 'Multiscale/TangentHom'}) 
d.update({806 : 'others/posProcElem'})  
d.update({900 : 'others/insertDeformation'})
d.update({901 : 'others/setMaterialParam'})
d.update({902 : 'others/setDamageParam'})

lib = {}
lib.update({-1 : 'others/executer'})
lib.update({-2 : 'others/loadingLib'})
lib.update({-3 : 'others/globalVariables'})
lib.update({-4 : 'others/ptsGaussLib'})
lib.update({-5 : 'others/ptsGaussLib2'})
lib.update({-6 : 'FiniteStrain/materialLib'})
lib.update({-7 : 'FiniteStrain/finiteStrainLib'})
lib.update({-8 : 'damage/damageLib'})
lib.update({-9 : 'damage/damageNewLib'})
lib.update({-10 : 'Multiscale/multiscaleLib'})

for e in elements:
	os.system('geany ' + SOLVERGPPsrc + d[e] + '.f90 &')

for l in libs:
	os.system('geany ' + SOLVERGPPsrc + lib[l] + '.f90 &')


