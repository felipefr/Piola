rm Iconverter.pyf
f2py -m Iconverter -h Iconverter.pyf Iconverter.f90
#~ f2py -c --fcompiler=gfortran --debug-capi Iconverter.pyf Iconverter.f90 -L./converterLib/ -I./converterLib/  -lconverter 
f2py -c --fcompiler=gfortran Iconverter.pyf Iconverter.f90 -L./converterLib/ -I./converterLib/  -lconverter 
 
