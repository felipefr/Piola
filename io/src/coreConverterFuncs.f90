module coreConverterFuncs
use utils
use auxConverterFuncs
use ConverterData

implicit none

private readFromGMSH, write2SGP, write2VTK
public MESH2SGP , SGP2VTK , readFromSGP , readSGPMesh

contains

! changes MESH2SGP to rotate mesh (Very messy impllementation )
subroutine MESH2SGP_withRotation()
	character(len=32) :: arg
	integer :: inputFace , inputDof, i , k , ip ,  NPeriodicParam
	integer, allocatable :: PeriodicParam(:,:)
	real*8 :: inputUd, rotation
		
	!======  Setting Neumann and Dirichlet conditions ===============
	if(iargc() /= 1) then
		call getarg(2, arg)
		read (arg,*) rotation, Ndirichlet,Nneumann,NPeriodicParam
		if(iargc() > 3 + Ndirichlet + NPeriodicParam) then
			call getarg(4 + Ndirichlet + NPeriodicParam , arg)
			read (arg,*)  Ndof , NParamCell , Nsubsteps , FindInterior, isPeriodicNoLagrange, addGlobalNode
			
			if(NParamCell > 0) withParamCell =.True. 
			
		end if
	end if
	write(0,*) "All the parameters set"

	call readFromGMSH()
	
	write(0,*) "Mesh Read"
	
	write(0,*) "Dirichlet conditions running"
	do i = 1, Ndirichlet
		call getarg(i+2, arg)
		read (arg,*) inputFace,inputDof,inputUd
		call numberingDirichlet(inputFace,inputDof,inputUd)
	end do
	write(0,*) "Dirichlet conditions set"

	write(0,*) "Neumann conditions running"
	if(Nneumann /=0) then
		allocate(BMneumannFace(Nneumann),NNeumannFace(Nneumann),NeumannFace(NFace))
		call getarg(3+Ndirichlet, arg)
		read (arg,*) BMneumannFace
		call numberingNeumann()
	end if
	write(0,*) "Neumann conditions set"
	
	write(0,*) "Periodic conditions running"
	do i = 1,NPeriodicParam
		if(i==1) allocate(PeriodicParam(NPeriodicParam,NodPer+1))
		call getarg(3+Ndirichlet+i, arg)
		read(arg,*) PeriodicParam(i,:)
		if(i==NperiodicParam) call setPeriodicConditions(PeriodicParam,NPeriodicParam)
	end do
	write(0,*) "Periodic conditions set"
	
!~ 	! ================================================================

	write(0,*) "Interior Found (if demanded) running"
	if(FindInterior) then 
		if(Cell_type == 10 .or. Cell_type == 24) then 
			call FindInteriorNode()
		else if(Cell_type == 5 .or. Cell_type == 22) then
			call FindInteriorNode2D()
		else if(Cell_type == 9 ) then
			call FindInteriorNode2DQuad()
		end if
	end if
	write(0,*) "Interior Found (if demanded) set"

!!! Bubbles 
!~ 	write(0,*) "Single to double field (if demanded) running"
!~ 	if(isP1bP1) call Single2DoubleFieldP1bP1()
!~ 	write(0,*) "Single to double field (if demanded) set"
	
	write(0,*) "Periodic Mesh No Lagrangean (if demanded) running"
	if(isPeriodicNoLagrange) call setPeriodicMesh()
	write(0,*) "Periodic Mesh No Lagrangean (if demanded) set"

	write(0,*) "Checking Jacobian"
!~ 	call checkJacobian()
	call checkJacobian_tri()
	write(0,*) "Jacobian checked"
	
	write(0,*) "Checking Jacobian"
!~ 	call checkJacobian()
	call checkJacobian_tri()
	write(0,*) "Jacobian checked"
	
	write(0,*) "Adding Global Node On Boundary (if demanded) running"
	if(addGlobalNode) call setGlobalNode()
	write(0,*) "Adding Global Node On Boundary (if demanded) set"	
	
	call rotateCoordinates(rotation)
	
	write(0,*) "SGP files written running"
	call write2SGP()
	write(0,*) "SGP files written finished"

end subroutine



! Converts TETGEN or GMSH format to SolverGP format
subroutine MESH2SGP(op)
	character(len=32) :: arg
	integer, intent(in) :: op
	integer :: inputFace , inputDof, i , k , ip ,  NPeriodicParam
	integer, allocatable :: PeriodicParam(:,:)
	real*8 :: inputUd
		
	!======  Setting Neumann and Dirichlet conditions ===============
	if(iargc() /= 1) then
		call getarg(2, arg)
		read (arg,*) Ndirichlet,Nneumann,NPeriodicParam
		if(iargc() > 3 + Ndirichlet + NPeriodicParam) then
			call getarg(4 + Ndirichlet + NPeriodicParam , arg)
			read (arg,*)  Ndof , NParamCell , Nsubsteps , FindInterior, isPeriodicNoLagrange, addGlobalNode
			
			if(NParamCell > 0) withParamCell =.True. 
			
		end if
	end if
	write(0,*) "All the parameters set"
		
	select case(op)
		case(1,2)
			call readFromGMSH()
	end select
	
	write(0,*) "Mesh Read"
	
	write(0,*) "Dirichlet conditions running"
	do i = 1, Ndirichlet
		call getarg(i+2, arg)
		read (arg,*) inputFace,inputDof,inputUd
		call numberingDirichlet(inputFace,inputDof,inputUd)
	end do
	write(0,*) "Dirichlet conditions set"

	write(0,*) "Neumann conditions running"
	if(Nneumann /=0) then
		allocate(BMneumannFace(Nneumann),NNeumannFace(Nneumann),NeumannFace(NFace))
		call getarg(3+Ndirichlet, arg)
		read (arg,*) BMneumannFace
		call numberingNeumann()
	end if
	write(0,*) "Neumann conditions set"
	
	write(0,*) "Periodic conditions running"
	do i = 1,NPeriodicParam
		if(i==1) allocate(PeriodicParam(NPeriodicParam,NodPer+1))
		call getarg(3+Ndirichlet+i, arg)
		read(arg,*) PeriodicParam(i,:)
		if(i==NperiodicParam) call setPeriodicConditions(PeriodicParam,NPeriodicParam)
	end do
	write(0,*) "Periodic conditions set"
	
!~ 	! ================================================================

	write(0,*) "Interior Found (if demanded) running"
	if(FindInterior) then 
		if(Cell_type == 10 .or. Cell_type == 24) then 
			call FindInteriorNode()
		else if(Cell_type == 5 .or. Cell_type == 22) then
			call FindInteriorNode2D()
		else if(Cell_type == 9 ) then
			call FindInteriorNode2DQuad()
		end if
	end if
	write(0,*) "Interior Found (if demanded) set"

!!! Bubbles 
!~ 	write(0,*) "Single to double field (if demanded) running"
!~ 	if(isP1bP1) call Single2DoubleFieldP1bP1()
!~ 	write(0,*) "Single to double field (if demanded) set"
	
	write(0,*) "Periodic Mesh No Lagrangean (if demanded) running"
	if(isPeriodicNoLagrange) call setPeriodicMesh()
	write(0,*) "Periodic Mesh No Lagrangean (if demanded) set"

	write(0,*) "Checking Jacobian"
!~ 	call checkJacobian()
	call checkJacobian_tri()
	write(0,*) "Jacobian checked"
	
	write(0,*) "Checking Jacobian"
!~ 	call checkJacobian()
	call checkJacobian_tri()
	write(0,*) "Jacobian checked"
	
	write(0,*) "Adding Global Node On Boundary (if demanded) running"
	if(addGlobalNode) call setGlobalNode()
	write(0,*) "Adding Global Node On Boundary (if demanded) set"	
	
	write(0,*) "SGP files written running"
	call write2SGP()
	write(0,*) "SGP files written finished"

end subroutine

! Converts from SolverGP format to VTK format
subroutine SGP2VTK()
	integer i
	character(len=32) :: arg
	
	if(iargc() > 1) then 
		
		call getarg(2, arg)
		read(arg,*) NtotalField , NParamCell , NvolGroups , NFaceGroups, isPeriodicNoLagrange , isXFEM, addGlobalNode
		
		if(NParamCell > 0) withParamCell = .True.
		
		allocate( iShiftFields(NtotalField), KindFields(NKindFields,NtotalField) , labelFields(NtotalField) )
		
		do i=1, NtotalField
			call getarg(2+i, arg)
			read(arg,*) labelFields(i) , IshiftFields(i), KindFields(1,i) , KindFields(2,i) !! KindFields(1,*) Point or Cell, KindFields(2,*) Vec or Scalar  
		end do
		
	end if

	call readFromSGP
	
!~ 	call Debug
	
	
	if(isXFEM) call modifyXFEM
	if(addGlobalNode) call subtractAdditionalGlobalNodes
	
	do i=1, NtotalField
		if(IshiftFields(i)<0) then
			IshiftFields(i) = -IshiftFields(i)  
			call TetraLin2TetraQuad(IshiftFields(i))
		end if 
	end do
	
	
	do i = 1, TimeSteps
		call write2VTK(i)
	end do
	
end subroutine

!writes SOLVERGP format (Mesh.txt and IniFile000.txt)
subroutine write2SGP()
	integer , parameter :: OUnitMesh = 22 , OUnitInit = 23, OUnitParam = 24
	integer :: i , ip , k , iPeriodic
	character(len=30) :: formatDOF, IformatDOF , IformatElem, IformatFace, formatNdim, IformatEdge, formatSetVar, IformatPer

	iPeriodic = 0
	if(NPeriodic>0) iPeriodic = 1
	
	open (OUnitMesh, file=MeshFile)  
	open (OUnitInit, file=IniFile0)
	open (OUnitParam, file=ParamFile0)
	
	write(formatDOF,"(a,I2,a)") "(", NDof, "(E24.16,2X))"
	write(IformatDOF,"(a,I2,a)") "(", NDof, "(I5))" 
	write(IformatElem,"(a,I2,a)") "(", NodEl, "(I6,2X))"
	write(IformatFace,"(a,I2,a)") "(", NodFace, "(I6,2X))"
	write(IformatEdge,"(a,I1,a)") "(", NodEdge, "(I5))"
	write(formatNdim,"(a,I1,a)") "(", Ndim, "(E24.16,2X))"
	write(IformatPer,"(a,I2,a)") "(", NodPer, "(I6,2X))"
	formatSetVar = "(a,/,I8)"
	
	write(OUnitMesh,formatSetVar) "*NODAL DOFs",NDof
	write(OUnitMesh,formatSetVar) "*DIMEN",NDim

	write(OUnitMesh,formatSetVar) "*COORDINATES",NNodes
	write(OUnitMesh,*) 
	write(OUnitMesh, formatNdim ) X ! Ndim assumed 3

	write(OUnitMesh,*) 
	if(AddGlobalNode) then 
		write(OUnitMesh,formatSetVar) "*ELEMENT GROUPS", size(NElemGroup) + Nneumann + iPeriodic + 1
	else
		write(OUnitMesh,formatSetVar) "*ELEMENT GROUPS", size(NElemGroup) + Nneumann + iPeriodic
	end if
	
	write(OUnitMesh,"(I2,1x,I8,1x,a)") ( i , NElemGroup(i) , "Generic" , i = 1, size(NElemGroup) ) 
	if(NNeumann /=0) write(OUnitMesh,"(I2,1x,I8,2x,a)") ( (size(NElemGroup) + i) , NneumannFace(i) , "Generic" , i=1,NNeumann) 
	if(NPeriodic /=0) write(OUnitMesh,"(I2,I5,2x,a)") size(NElemGroup) + NNeumann + 1 , NPeriodic , "Generic"  
	if(AddGlobalNode) write(OUnitMesh,"(I2,I5,2x,a)") size(NElemGroup) + NNeumann + iPeriodic + 1 , 1 , "Generic"  !!! Global element
	
	write(OUnitMesh,"(I2)") ( NodEl , i = 1, NElem)
	if(Nneumann /=0) write(OUnitMesh,"(I2)") ( NodFace , i = 1, sum(NneumannFace))
	if(NPeriodic /=0) write(OUnitMesh,"(I2)") ( NodPer , i = 1, NPeriodic)
	if(AddGlobalNode) write(OUnitMesh,"(I2)") 1

	write(OUnitMesh,*) 
	write(OUnitMesh,"(a)") "*INCIDENCE"
	write(OUnitMesh, IformatElem ) Elem ! assumed 4

	if(Nneumann /=0) then
		do k = 1, NFace 
			if(NeumannFace(k) /=0) then
				ip = (NeumannFace(k) - 1)*NodFace
				write(OUnitMesh,IformatFace) Face(ip + 1 : ip+NodFace)
			end if
		end do
		
	end if

	if(NPeriodic /=0) write(OUnitMesh,IformatPer) Periodic
	if(AddGlobalNode) write(OUnitMesh,"(I4)") Nnodes
	
	write(OUnitMesh,*) 
	write(OUnitMesh,"(a)") "*ELEMENT TYPE"
	write(OUnitMesh,"(I2)") ( ElemGroup(i) , i = 1, NElem)

	do i = 1 , Nneumann
		write(OUnitMesh,"(I2)") ( (i + size(NElemGroup)) , k = 1, NneumannFace(i))
	end do
	
	if(NPeriodic /=0) write(OUnitMesh,"(I2)") ( Nneumann + size(NElemGroup) + 1 , k=1 , NPeriodic )
	
	if(AddGlobalNode) write(OUnitMesh,"(I4)") size(NElemGroup) + Nneumann + iPeriodic + 1
	
	write(OUnitMesh,*) 
	write(OUnitMesh,"(a)") "*ELEMENT MAT"
	
	if(withParamCell) then
		write(OUnitMesh,"(I7)") ( 1+i , i = 1, NElem)
	else
		write(OUnitMesh,"(I1)") ( 1 , i = 1, NElem)
	end if
	
	do i = 1 , Nneumann
		write(OUnitMesh,"(I1)") ( 1 , k = 1, NneumannFace(i))
	end do

	if(NPeriodic /=0) write(OUnitMesh,"(I1)") ( 1 , k=1 , NPeriodic )
	if(AddGlobalNode) write(OUnitMesh,"(I1)") 1
	
	write(OUnitMesh,*) 
	write(OUnitMesh,"(a)") "*DIRICHLET CONDITIONS"
	do i = 1 , Nsubsteps
		write(OUnitMesh,IformatDOF) IDD
		write(OUnitMesh,*) 
	end do
	
	do i = 1 , Nsubsteps
		write(OUnitMesh,formatDOF) DValues
		write(OUnitMesh,*)
	end do
	
	write(OUnitInit,"(a)") "*Initial Conditions"
	do i = 1 , Nsubsteps
		write(OUnitInit,formatDOF) DValues
		write(OUnitInit,*)
	end do
	
	if(withParamCell) then
		write(OUnitParam,formatSetVar) "*Parameter Groups" , NElem + 1  ! it has always a neutral in order to link to others elements
		write(OUnitParam,"(a)") "*Real Parameters"
		write(OUnitParam,"(I3)") 0
		write(0,*) "Gauss point number", NGP
		write(OUnitParam,"(I3)") ( NParamCell*NGP , i = 1 , Nelem ) 
		write(OUnitParam,"(E24.16)")  ( 0.0, i = 1, NParamCell*NGP*NElem)
		write(OUnitParam,"(a)") "*Integer Parameters"
		write(OUnitParam,"(I6,'*',I1)") NElem + 1, 0
	else 
		write(OUnitParam,formatSetVar) "*Parameter Groups" , 1  ! it has always a neutral in order to link to others elements
		write(OUnitParam,"(a)") "*Real Parameters"
		write(OUnitParam,"(I3)") 0
		write(OUnitParam,"(a)") "*Integer Parameters"
		write(OUnitParam,"(I6,'*',I1)") 1, 0
	end if
	
	close(OUnitMesh)
	close(OUnitInit)
	close(OUnitParam)

end subroutine

subroutine readSGPMesh()
	integer , parameter :: IUnitMesh = 22
	integer :: i , j , elementsGroup, Idummy
	real*8 :: dummy
	
	open (IUnitMesh, file=MeshFile, status= 'old' )  
	
	! Reading Mesh
	call FindKeyword(IUnitMesh,"*NODAL DOFs")
	read(IUnitMesh,*) NDof
	call FindKeyword(IUnitMesh,"*DIMEN")
	read(IUnitMesh,*) Ndim
	call FindKeyword(IUnitMesh,"*COORDINATES")
	read(IUnitMesh,*) NNodes
	allocate(X(NNodes*Ndim))
	read(IUnitMesh,*) X
	call FindKeyword(IUnitMesh,"*ELEMENT GROUPS")
	read(IUnitMesh,*) elementsGroup
	NElem = 0
	NFace = 0
	write(*,*) NvolGroups
	do i = 1 , NvolGroups
		read(IUnitMesh,*) dummy , Idummy
		NElem = NElem + Idummy
	end do 
	do i = 1 , NfaceGroups
		read(IUnitMesh,*) dummy , Idummy 
		NFace = NFace + Idummy 
	end do
	
	do i = 1 , elementsGroup - NvolGroups - NfaceGroups
		read(IUnitMesh,*) 
	end do
	
	read(IUnitMesh,*) NodEL
	do i = 1, Nelem -1
		read(IUnitMesh,*)
	end do
	
	if(NfaceGroups /=0) read(IUnitMesh,*) NodFace
	
	allocate(Elem(NElem*NodEL),ElemGroup(NElem),Face(NFace*NodFace))
	call FindKeyword(IUnitMesh,"*INCIDENCE")
	read(IUnitMesh,*) Elem
	read(IUnitMesh,*) Face
	
	if(isPeriodicNoLagrange) call doubledMesh2Single()
	
	call FindKeyword(IUnitMesh,"*ELEMENT TYPE")
	read(IUnitMesh,*) ElemGroup
	
	close(IUnitMesh)

end subroutine


! Reads solverGP format (Mesh.txt and DataOut.txt (if any) )
subroutine readFromSGP()
	integer , parameter :: IUnitDataOut = 23, IUnitParam = 24
	integer :: i,j
	real*8 :: dummy
	
	call readSGPMesh()
	
	if(addGlobalNode) call subtractGlobalNodeInVolume()
	
	call selectCellType()
	
	open (IUnitDataOut, file=DataOutFile, status= 'old' )
	open (IUnitParam, file=ParamFile, status= 'old' )
	
	call countKeyword(TimeSteps,IUnitDataOut,"*Time")
		
	allocate(Sol(NNodes*NDof,TimeSteps),Dt(TimeSteps),Tini(TimeSteps))
	allocate(Param(NElem*NParamCell*NGP,TimeSteps))
	
	write(0,*) "NElem = " , NElem , "NFace =" , NFace, "NPeriodic =" , NPeriodic
	
	do i = 1 , TimeSteps
		call FindKeyword(IUnitDataOut,"*Time")
		read(IUnitDataOut,*) dummy , Tini(i) , Dt(i) 
		read(iUnitDataOut,*) Sol(:,i)
		
		if(withParamCell) then
			call FindKeyword(IUnitParam,"*Real Parameters")
			read(IUnitParam,*) ( dummy , j=1,NElem + 1)
			read(iUnitParam,*) Param(:,i)
		end if
	end do
	
	close(IUnitDataOut)
	close(IUnitParam)
	
end subroutine

! Writes VTK format
subroutine write2VTK(step) 
	integer , parameter :: OUnitVTK = 45
	character(len=80) :: formatDOF, IformatElem, formatNdim,  formatSetVar, formatSetVar2 , nameFile
	integer :: i , step, ip , j , k
	real*8 :: rAux
	
	write(IformatElem,"(a,I2,a)") "(", NodEl_vtk + 1, "(I15))"
	
	if(Ndim == 3) then
		write(formatNdim,"(a,I1,a)") "(", Ndim, "(E24.16,2X))"
	else if(Ndim == 2) then 
		write(formatNdim,"(a,I1,a)") "(", Ndim, "(E24.16,2X),'0.0000')"
	else
		write(*,*) "Dimension must be 2 or 3"
	end if
	
	formatSetVar = "(a,1X,I8)"
	formatSetVar2 = "(a,1X,I8,1X,a)"
		
	namefile = trim(VTKfilePrefix) // trim(str(step)) // ".vtk"
	
	open (OUnitVTK, file=nameFile )
	
	call selectCellType() !!! there is already in readFromSGP, see if it is really necessary
			
	write(OUnitVTK,"(a)") "# vtk DataFile Version 2.0"
	write(OUnitVTK,"(a)") "Unstructured Grid"
	write(OUnitVTK,"(a)") "ASCII"
	write(OUnitVTK,"(a)") "DATASET UNSTRUCTURED_GRID"
	
!~ 	if(isP1bP1) NNodes = NNodes - NElem !!! bubbles
	
	write(OUnitVTK,formatSetVar2) "POINTS", NNodes , "double"
	write(OUnitVTK,formatNdim) X(1:Ndim*NNodes)
	
	write(OUnitVTK,*)
	write(OUnitVTK,"(a,2X,I15,2X,I15)" ) "CELLS", NElem , NElem*(1+NodEl_vtk)
	write(OUnitVTK,IformatElem) ( NodEL_vtk , Elem( (i-1)*NodEL + 1 : (i-1)*NodEL + NodEl_vtk ) - 1, i = 1 , Nelem ) 
	
!~ 	if(isP1bP1) then !!! Bubbles ( if - else starting two lines before )
!~ 		write(OUnitVTK,"(a,2X,I15,2X,I15)" ) "CELLS", NElem , NElem*(1+NodEl-1)  ! ghost node in the end
!~ 		write(OUnitVTK,IformatElem) ( NodEL-1 , Elem( (i-1)*NodEL + 1 : (i-1)*NodEL + NodEl-1 ) - 1, i = 1 , Nelem )		

	write(OUnitVTK,*) 
	write(OUnitVTK,formatSetVar) "CELL_TYPES", Nelem
	write(OUnitVTK,"(I5)") ( Cell_type , i = 1 , Nelem )
	write(OUnitVTK,*) 

	!! Writing Fields	
	write(OUnitVTK,formatSetVar) "CELL_DATA", Nelem
	
	!! First field just the group of elements
	write(OUnitVTK,"(a)") "SCALARS groupId float"
	write(OUnitVTK,"(a)") "LOOKUP_TABLE default"
	write(OUnitVTK,"(I5)") ( ElemGroup(i) , i= 1 , Nelem ) 
	
	do k = 1 , NtotalField
		if(KindFields(1,k)) cycle !! if point get out, must be false to continue as CELL
		if(KindFields(2,k)) cycle !! if Vec get out, must be false to continue as SCALAR ( not yet implemented for Vec in CELLS)
		
		write(OUnitVTK,"(a,a,a)") "SCALARS ", labelFields(k), " float"
		write(OUnitVTK,"(a)") "LOOKUP_TABLE default"
			
		write(OUnitVTK,"(E24.16)") ( Param((i-1)*NParamCell+1+IshiftFields(k), step) , i= 1 , Nelem ) 
		
	end do
			
	write(OUnitVTK,formatSetVar) "POINT_DATA", NNodes
	
	do k = 1 , NtotalField
		if(.not. KindFields(1,k)) cycle !! have to be point
		
		if(KindFields(2,k)) then !! is vector field
			write(OUnitVTK,"(a,a,a)") "VECTORS ", labelFields(k) , " float"
!~ 			write(OUnitVTK,formatNdim) ( Sol((i-1)*Ndof + IshiftFields(k) + 1 : &
!~ 											  (i-1)*Ndof + IshiftFields(k) + Ndim, step) , i= 1 , NNodes ) !! Original
			write(OUnitVTK,formatNdim) ( Sol((i-1)*Ndof + IshiftFields(k) + 1 : &
										(i-1)*Ndof + IshiftFields(k) + 2, step) , 0.0d0  , i= 1 , NNodes ) !! aconchabramento para 2D
											  
		else
			write(OUnitVTK,"(a,a,a)") "SCALARS ", labelFields(k), " float"
			write(OUnitVTK,"(a)") "LOOKUP_TABLE default"
			write(OUnitVTK,"(E24.16)") ( Sol((i-1)*Ndof + IshiftFields(k) + 1 , step) , i= 1 , NNodes )
		end if
	
	end do
	
	close(OUnitVTK)
	
end subroutine


subroutine readFromGMSH() 
	real*8 :: dummy
	integer , parameter :: IUnitMesh = 12 , MaxParamElem = 5
	integer :: i , Idummy , NelemAux , kFace,kEdge,kElem, kNElemGroup, j , k
	character(len=32) :: arg
	integer, allocatable:: iMaux(:,:)
	integer :: countElem(5) !! Line (1) , Triangle(2) , Tetrahedra(4) , Points(15)  
	
	!! IMaux(1,:) = element index (dummy)
	!! IMaux(2,:) = element type 
	!! IMaux(3,:) = element number of nodes (dummy)  
	!! IMaux(4,:) = element physical group
	!! IMaux(5,:) = element entities (original label of the lines) (dummy) 
	
	
	open (IUnitMesh, file=GMSHFile, status= 'old' )  
	Ndim = 3 !! always 3D
	
	! Reading nodes
	call FindKeyword(IUnitMesh,"$Nodes")
	read(IUnitMesh,*) NNodes
	write(*,*) NNodes,Ndim
	allocate(X(NNodes*Ndim),IDD(NNodes*NDof),DValues(NNodes*NDof))
	Read (iUnitMesh,*) ( dummy , X((i-1)*Ndim + 1: i*Ndim ) , i=1 , NNodes )
	
	IDD = 0
	DValues = 0.0d0 
	
	call FindKeyword(IUnitMesh,"$Elements")
	read(IUnitMesh,*) NelemAux
	allocate(iMaux(MaxParamElem,NelemAux))
	
	do i=1,NelemAux
		read(iUnitMesh,*) iMaux(1:MaxParamElem,i)		
	end do
		
	if(any(iMaux(2,:)==4)) then ! linear tetrahedra
		NElem = count(iMaux(2,:)==4)
		NFace = count(iMaux(2,:)==2) ! linear triangles
		NEdge = count(iMaux(2,:)==1) !! no use by the moment
		NodEl = 4
	elseif(any(iMaux(2,:)==11)) then ! quadratic tetrahedras
		NElem = count(iMaux(2,:)==11) 
		NFace = count(iMaux(2,:)==9)! quadratic triangles
		NEdge = count(iMaux(2,:)==1) !! no use by the moment
		NodEl = 10
	elseif(any(iMaux(2,:)==2)) then ! linear triangles
		NElem = count(iMaux(2,:)==2) 
		NFace = count(iMaux(2,:)==1) ! linear lines
		NodEl = 3
	elseif(any(iMaux(2,:)==9)) then ! quadratic triangles
		NElem = count(iMaux(2,:)==9) 
		NFace = count(iMaux(2,:)==8) ! quadratic lines
		NodEl = 6
	elseif(any(iMaux(2,:)==3)) then ! linear quadrangles
		NElem = count(iMaux(2,:)==3) 
		NFace = count(iMaux(2,:)==1) ! linear lines
		NodEl = 4
	elseif(any(iMaux(2,:)==10)) then ! quadratic quadrangles
		NElem = count(iMaux(2,:)==10) 
		NFace = count(iMaux(2,:)==8) ! quadratic lines
		NodEl = 9
	elseif(any(iMaux(2,:)==1)) then ! linear lines
		write(0,*) "linear lines"
		call myPrint(iMaux(2,:))
		NElem = count(iMaux(2,:)==1) 
		NFace = count(iMaux(2,:)==15) ! points
		NodEl = 2
	end if

	call selectCellType() !! sets NodFace, NodEdge, NodEl_vtk automatically
		
	write(0,*) "NElem = ", NElem, "NodEl = ", NodEl, "cellType = ", cell_Type , "NFace =", NFace, "NodFace =", NodFace
		
	allocate(Elem(NElem*NodEl),ElemGroup(NElem))
	allocate(Face(NFace*NodFace),BMFace(NFace))
	allocate(Edge(NEdge*NodEdge))

	rewind(IUnitMesh)
	call FindKeyword(IUnitMesh,"$Elements")
	read(IUnitMesh,*)
	
	kElem = 0
	kFace = 0
	kEdge = 0
	kNElemGroup = 0
	ElemGroup = -100
	
	do i = 1 , NelemAux
		
		if(cell_type == 10 .or. cell_type == 24) then !Linear or quadratic tetrahedra as elements
			select case(iMaux(2,i))
				case(4,11,15)
	!~ 				write(*,*) "Element" 
					if(kElem==0 .or. count(ElemGroup(1:kElem) == iMaux(4,i))==0 ) kNElemGroup = kNElemGroup + 1	
					kElem = kElem + 1
					ElemGroup(kElem) = iMaux(4,i)
	!~ 				write(*,*) ElemGroup(kElem)
					read(IUnitMesh,*) (dummy, j=1,MaxParamElem) , Elem( (kElem-1)*NodEl + 1 :  kElem*NodEl )  
				case(2,9)
	!~ 				write(*,*) "Face"
					kFace = kFace + 1
					BMFace(kFace) = iMaux(4,i) 
					read(IUnitMesh,*) (dummy, j=1,MaxParamElem) , Face( (kFace-1)*NodFace + 1 : kFace*NodFace )
				case(1)
	!~ 				write(*,*) "Edge"
					kEdge = kEdge + 1
					read(IUnitMesh,*) (dummy, j=1,MaxParamElem) , Edge( (kEdge-1)*NodEdge + 1 : kEdge*NodEdge )
					write(*,*) Edge( (kEdge-1)*NodEdge + 1 : kEdge*NodEdge )
				case default
					write(*,*) "GMSH Element type unknown", iMaux(2,i) 
					read(IUnitMesh,*) 
			end select
		
		elseif(cell_type == 5) then ! Linear triangle as elements
			select case(iMaux(2,i))
				case(2)
					if(kElem==0 .or. count(ElemGroup(1:kElem) == iMaux(4,i))==0 ) kNElemGroup = kNElemGroup + 1	
					kElem = kElem + 1
					ElemGroup(kElem) = iMaux(4,i)
					read(IUnitMesh,*) (dummy, j=1,MaxParamElem) , Elem( (kElem-1)*NodEl + 1 :  kElem*NodEl )  
				case(1)
					kFace = kFace + 1
					BMFace(kFace) = iMaux(4,i) 
					read(IUnitMesh,*) (dummy, j=1,MaxParamElem) , Face( (kFace-1)*NodFace + 1 : kFace*NodFace )
				case default
					write(*,*) "GMSH Element type unknown", iMaux(2,i) 
					read(IUnitMesh,*) 
			end select		
		elseif(cell_type == 22) then ! quadratic triangles as elements
			select case(iMaux(2,i))
				case(9)
					if(kElem==0 .or. count(ElemGroup(1:kElem) == iMaux(4,i))==0 ) kNElemGroup = kNElemGroup + 1	
					kElem = kElem + 1
					ElemGroup(kElem) = iMaux(4,i)
					read(IUnitMesh,*) (dummy, j=1,MaxParamElem) , Elem( (kElem-1)*NodEl + 1 :  kElem*NodEl )  
				case(8)
					kFace = kFace + 1
					BMFace(kFace) = iMaux(4,i) 
					read(IUnitMesh,*) (dummy, j=1,MaxParamElem) , Face( (kFace-1)*NodFace + 1 : kFace*NodFace )
				case default
					write(*,*) "GMSH Element type unknown", iMaux(2,i) 
					read(IUnitMesh,*) 
			end select		
		elseif(cell_type == 3) then ! Linear Lines as elements
			select case(iMaux(2,i))
				case(1)
					if(kElem==0 .or. count(ElemGroup(1:kElem) == iMaux(4,i))==0 ) kNElemGroup = kNElemGroup + 1	
					kElem = kElem + 1
					ElemGroup(kElem) = iMaux(4,i)
					read(IUnitMesh,*) (dummy, j=1,MaxParamElem) , Elem( (kElem-1)*NodEl + 1 :  kElem*NodEl )  
				case(15)
					kFace = kFace + 1
					BMFace(kFace) = iMaux(4,i) 
					read(IUnitMesh,*) (dummy, j=1,MaxParamElem) , Face( (kFace-1)*NodFace + 1 : kFace*NodFace )
				case default
					write(*,*) "GMSH Element type unknown", iMaux(2,i) 
					read(IUnitMesh,*) 
			end select
		elseif(cell_type == 9) then ! Linear quadrangles as elements
			select case(iMaux(2,i))
				case(3)
					if(kElem==0 .or. count(ElemGroup(1:kElem) == iMaux(4,i))==0 ) kNElemGroup = kNElemGroup + 1	
					kElem = kElem + 1
					ElemGroup(kElem) = iMaux(4,i)
					read(IUnitMesh,*) (dummy, j=1,MaxParamElem) , Elem( (kElem-1)*NodEl + 1 :  kElem*NodEl )  
				case(1)
					kFace = kFace + 1
					BMFace(kFace) = iMaux(4,i) 
					read(IUnitMesh,*) (dummy, j=1,MaxParamElem) , Face( (kFace-1)*NodFace + 1 : kFace*NodFace )
				case default
					write(*,*) "GMSH Element type unknown", iMaux(2,i) 
					read(IUnitMesh,*) 
			end select		
		elseif(cell_type == 23) then ! quadratic quadrangles as elements
			select case(iMaux(2,i))
				case(10)
					if(kElem==0 .or. count(ElemGroup(1:kElem) == iMaux(4,i))==0 ) kNElemGroup = kNElemGroup + 1	
					kElem = kElem + 1
					ElemGroup(kElem) = iMaux(4,i)
					read(IUnitMesh,*) (dummy, j=1,MaxParamElem) , Elem( (kElem-1)*NodEl + 1 :  kElem*NodEl )  
				case(8)
					kFace = kFace + 1
					BMFace(kFace) = iMaux(4,i) 
					read(IUnitMesh,*) (dummy, j=1,MaxParamElem) , Face( (kFace-1)*NodFace + 1 : kFace*NodFace )
				case default
					write(*,*) "GMSH Element type unknown", iMaux(2,i) 
					read(IUnitMesh,*) 
			end select		
		end if
	end do
	
	allocate(NElemGroup(kNElemGroup))
	
	write(*,*) kNElemGroup	

	k = 1
	do i = 1 , kNElemGroup
		if(k==1) then 
			Idummy = ElemGroup(k)
			k = k + 1
		else
			do j = k , NElem
				Idummy = ElemGroup(j)
				if( .not. ( any(ElemGroup(1:j-1) == Idummy ) ) ) exit
			end do
			k = j + 1
		end if
		NElemGroup(i) = count(ElemGroup == Idummy)
	end do
			
!~ 	call myPrint(NElemGroup)

	if (Cell_type == 24) call convert2StandardNumeration() !! Maybe it does not need
	
	k = 0
	call CountKeyword(k,IUnitMesh,"$Periodic")
	close (IUnitMesh)
	
	if(k>0 .and. .not. isPeriodicNoLagrange) call setPeriodicConditionsAuto()
	
end subroutine

subroutine SGP2VTK_inifile()
	integer i
	character(len=32) :: arg
	
	if(iargc() > 1) then 
		
		call getarg(2, arg)
		read(arg,*) NtotalField , NParamCell , NvolGroups , NFaceGroups, isPeriodicNoLagrange , isXFEM, addGlobalNode
		
		if(NParamCell > 0) withParamCell = .True.
		
		allocate( iShiftFields(NtotalField), KindFields(NKindFields,NtotalField) , labelFields(NtotalField) )
		
		do i=1, NtotalField
			call getarg(2+i, arg)
			read(arg,*) labelFields(i) , IshiftFields(i), KindFields(1,i) , KindFields(2,i) !! KindFields(1,*) Point or Cell, KindFields(2,*) Vec or Scalar  
		end do
		
	end if

	call readFromSGP_inifile
	
	if(isXFEM) call modifyXFEM
	if(addGlobalNode) call subtractAdditionalGlobalNodes
	
	do i=1, NtotalField
		if(IshiftFields(i)<0) then
			IshiftFields(i) = -IshiftFields(i)  
			call TetraLin2TetraQuad(IshiftFields(i))
		end if 
	end do
	
	TimeSteps = 1

	call write2VTK(1)	
	
end subroutine

subroutine readFromSGP_inifile()
	integer , parameter :: IUnitInifile = 23, IUnitParam = 24
	integer :: i,j
	real*8 :: dummy
	
	call readSGPMesh()
	
	if(addGlobalNode) call subtractGlobalNodeInVolume()
	
	call selectCellType()
	
	open (IUnitInifile, file=IniFile, status= 'old' )
	open (IUnitParam, file=ParamFile, status= 'old' )
		
	TimeSteps = 1
	
	allocate(Sol(NNodes*NDof,TimeSteps),Dt(TimeSteps),Tini(TimeSteps))
	allocate(Param(NElem*NParamCell*NGP,TimeSteps))
	
	write(0,*) "NElem = " , NElem , "NFace =" , NFace, "NPeriodic =" , NPeriodic
	
	do i = 1 , TimeSteps
		call FindKeyword(IUnitInifile,"*Initial Conditions")
		read(iUnitInifile,*) Sol(:,i)
		
		if(withParamCell) then
			call FindKeyword(IUnitParam,"*Real Parameters")
			read(IUnitParam,*) ( dummy , j=1,NElem + 1)
			read(iUnitParam,*) Param(:,i)
		end if
	end do
	
	close(IUnitInifile)
	close(IUnitParam)
	
end subroutine

end module 
