module ConverterData

implicit none
public Init, KillData

real*8 , allocatable :: X(:) , DValues(:) , Sol(:,:) , Dt(:) , Tini(:) , Param(:,:)
integer , allocatable :: Edge(:) ,  Elem(:) , Face(:) , BMface(:) , IDD(:) , Periodic(:,:)
integer , allocatable ::  BMneumannFace(:) , NeumannFace(:) , NneumannFace(:), ElemGroup(:) , NElemGroup(:) , iShiftFields(:) 
logical , allocatable :: KindFields(:,:) !! isPoint , isVec
character(len=10) , allocatable :: labelFields(:) 

integer :: NNodes, NEdge, NElem, NodEl, NodEl_vtk, NodEdge, NodFace, NFace, Ndim , NdimE, NDof, TimeSteps , Cell_type
integer :: Ndirichlet, Nneumann , Nsubsteps, NGP,NParamCell, NvolGroups, NfaceGroups , Nperiodic ,NtotalField 
logical :: FindInterior, withParamCell, isQuadratic, isPeriodicNoLagrange, isXFEM, addGlobalNode
integer , parameter :: NGPlinear = 1 , NGPquad = 4 ,  NKindFields = 2 , NodPer = 2
character(len = 50) :: ParamFile, IniFile, MeshFile, DataOutFile, VTKfilePrefix, ParamFile0, IniFile0, GMSHfile

contains
subroutine Init()

	FindInterior = .False.
	withParamCell = .False.
	isPeriodicNoLagrange = .false.
	isXFEM = .false.
	addGlobalNode = .false.
	
	Nsubsteps = 1
	NParamCell = 0
	NvolGroups = 1
	NfaceGroups = 0
	NtotalField = 0
	NNodes = 0
	NEdge = 0
	NElem = 0
	NFace = 0
	NodEl = 0
	NodEl_vtk = 0
	NodEdge = 2
	NodFace = 3
	
	Ndim = 3
	NdimE = 2
	NDof = 1
	TimeSteps = 0 
	Cell_type = 0

	Ndirichlet = 0
	Nneumann = 0
	Nperiodic = 0
	
	NGP = 1		

	MeshFile = "Mesh.txt"
	IniFile = "IniFile.txt"
	ParamFile = "Param.txt"
	DataOutFile = "DataOut.txt"
	VTKfilePrefix = "mesh"
	IniFile0 = "IniFile000.txt"
	ParamFile0 = "Param000.txt"
	GMSHfile = "mesh.msh"

end subroutine

subroutine KillData()

	if (allocated(X)) deallocate(X)
	if (allocated(Sol)) deallocate(Sol)
	if (allocated(Dt)) deallocate(Dt)
	if (allocated(Tini)) deallocate(Tini)
	if (allocated(Param)) deallocate(Param)
	if (allocated(Elem)) deallocate(Elem)
	if (allocated(Face)) deallocate(Face)
	if (allocated(BMFace)) deallocate(BMFace)
	if (allocated(dvalues)) deallocate(dvalues)
	if (allocated(Edge)) deallocate(Edge)
	if (allocated(Periodic)) deallocate(Periodic)
	if (allocated(BMneumannFace)) deallocate(BMneumannFace)
	if (allocated(NeumannFace)) deallocate(NeumannFace)
	if (allocated(ElemGroup)) deallocate(ElemGroup)
	if (allocated(NElemGroup)) deallocate(NElemGroup)
	if (allocated(IshiftFields)) deallocate(IshiftFields)
	if (allocated(KindFields)) deallocate(KindFields)
	
end subroutine

end module
