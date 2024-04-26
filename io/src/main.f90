!!! Instructions (18/08)
! arg1 = 1 , 2 , 3 , 4 = TETGEN2SGP, GMSH2SGP, SGP2VTK, Write2Param
! * TETGEN2SGP and GMSH2SGP
! arg2 = Ndirichlet,Nneumann,NPeriodicParam
! next the NDirichlet conditions in the format : 'Face dof value'
! next the NNeumman conditions in the format : 'Face1 ... FaceNNeumman' , ps: write '' for NNeumman = 0
! next the NPeriodic conditions  in the format : 'Face1 Face2 axis' 
! argLast = Ndof , NParamCell , Nsubsteps , IisPxPy, FindInterior, isPeriodicNoLagrange (logical for the last ones)
! ps : IisPxPy : 0 P1 or P2, 1 P1P1(stab) , 2 P2P1 , 3P1bP1
! * GMSH2SGP
! arg2 = NtotalField , NParamCell , NvolGroups , IisPxPy , isPeriodicNoLagrange (logical for the last ones)
! next NtotalField in the format : 'labelField , IshiftField , isPoint , isVec' (Logical the last ones) 
! * Write2Param
! arg2 = NGP,NdamageParam2,kindDamage, NinputDamage 
! arg3 = inputDamage(1) ... inputDamage(NinputDamage) 
! ps: kindDamage : 1 constantDamage , 2 slopeDamage, 3 constantDamagePlate , 4 slopeDamagePlate

! MAIN PROGRAM


program converter    
use coreConverterFuncs
use writeParam

implicit none

character(len=32) :: arg
integer :: op 

call init()

if (iargc() == 0) then
	write(*,*) "Invalid option : To TETGEN2SGP type 1 , to GMSH2SGP type 2, to SGP2VTK type 3, to writeParam type 4"
else
	call getarg(1, arg)
	read(arg,*) op
	
	select case(op)
		case(1,2)
			call MESH2SGP(op)
		case(3)
			call SGP2VTK	
		case(4)
			call write2Param
		case(5)
			call SGP2VTK_inifile
		case(6)
			call homogenisationFibres
		case(7)
			call MESH2SGP_withRotation
		case default
			write(*,*) "Invalid option : To TETGEN2SGP type 1 , to GMSH2SGP type 2, to SGP2VTK type 3,&
						to writeParam type 4"
	end select
		
end if

end program    



