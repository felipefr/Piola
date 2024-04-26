module auxConverterFuncs
use utils
use converterData

implicit none

private findSym, SubtractTriFromTet, SubtractLineFromTri, SubtractLineFromQuad, setGlobalNodeInFace, setGlobalNodeInVolume, fixNode

public selectCellType, convert2StandardNumeration, numberingDirichlet, numberingNeumann, setPeriodicConditions, &
	   FindInteriorNode, setPeriodicMesh, getXele, findBaricentric, TetraLin2TetraQuad, checkJacobian, modifyXFEM, &
	   FindInteriorNode2D, setPeriodicConditionsAuto, subtractGlobalNodeInVolume, & !! Single2DoubleFieldP1bP1  ( for bubbles )
	   getArraySolEl, getArrayXele, Debug, checkJacobian_tri, FindInteriorNode2DQuad, rotateCoordinates

contains

subroutine Debug
!~ 
!~ real*8 , allocatable :: X(:) , DValues(:) , Sol(:,:) , Dt(:) , Tini(:) , Param(:,:)
!~ integer , allocatable :: Edge(:) ,  Elem(:) , Face(:) , BMface(:) , IDD(:) , Periodic(:,:)
!~ integer , allocatable ::  BMneumannFace(:) , NeumannFace(:) , NneumannFace(:), ElemGroup(:) , NElemGroup(:) , iShiftFields(:) 
!~ logical , allocatable :: KindFields(:,:) !! isPoint , isVec
!~ character(len=10) , allocatable :: labelFields(:) 
!~ 
!~ integer :: NNodes, NEdge, NElem, NodEl, NodEdge, NodFace, NFace, Ndim , NDof, TimeSteps , Cell_type
!~ integer :: Ndirichlet, Nneumann , Nsubsteps, NGP,NParamCell, NvolGroups, NfaceGroups , Nperiodic ,NtotalField 
!~ logical :: FindInterior, withParamCell, isQuadratic, isPeriodicNoLagrange, isXFEM, addGlobalNode
!~ integer , parameter :: NGPlinear = 1 , NGPquad = 4 ,  NKindFields = 2 , NodPer = 2


write(*,*) "NNodes = " , NNodes
write(*,*) "NElem = " , NElem
write(*,*) "NodEl = " , NodEl
write(*,*) "NodFace = " , NodFace
write(*,*) "NFace = " , NFace

write(*,*) "Elem =  "
call myPrint(Elem)

write(*,*) "X =  "
call myPrint(X)

write(*,*) "Face =  "
call myPrint(Face)

end subroutine

subroutine rotateCoordinates(thetaRot)
	real*8, intent(in) :: thetaRot
	real*8 :: Qrot(Ndim,Ndim)
	integer :: i , ipDim
	
	Qrot = 0.0d0
	Qrot(1,1) = dcos(thetaRot)
	Qrot(2,2) = dcos(thetaRot)
	Qrot(1,2) = -dsin(thetaRot)
	Qrot(2,1) = dsin(thetaRot)
	Qrot(3,3) = 1.0d0
	
	do i = 1, Nnodes
		ipDim = (i-1)*Ndim 
		X(ipDim + 1 : ipDim + Ndim) = matmul(Qrot, X(ipDim + 1 : ipDim + Ndim))
	end do

end subroutine



subroutine subtractAdditionalGlobalNodes
	integer :: i , j, ip, ipp, ipOld, ippOld, NNodesOld
	real*8 :: XOld(NNodes*Ndim), SolOld(NNodes*Ndof,TimeSteps)
	
	Xold = X
	SolOld = Sol

	NnodesOld = Nnodes
 	
	Nnodes = NnodesOld - 2

	deallocate(X,Sol)
	allocate(X(NNodes*Ndim),Sol(NNodes*Ndof,TimeSteps))
	
	X = XOld(1:Nnodes*Ndim)
	Sol = SolOld(1:Nnodes*Ndof,:)
	
end subroutine

subroutine subtractGlobalNodeInVolume()
	integer :: i , j, ip, ipp, ipOld, ippOld, NodElOld
	integer :: ElemOld(NElem*NodEl)
	
	ElemOld = Elem

	NodElOld = NodEl
 	
	NodEl = NodElOld - 1
	
	deallocate(Elem)
	allocate(Elem(NElem*NodEl))

	do i = 1 , NElem
		ip = (i-1)*NodEl
		ipOld = (i-1)*NodElOld
		
		do j = 1,NodEl
			ipp = ip + j
			ippOld = ipOld + j
			Elem(ipp) = ElemOld(ippOld) 
		end do
	end do	
	
end subroutine

subroutine setGlobalNodeInVolume(idGlobal)
	integer ,  intent(in):: idGlobal
	integer :: i , j, ip, ipp, ipOld, ippOld, NodElOld
	integer :: ElemOld(NElem*NodEl)
	
	ElemOld = Elem

	NodElOld = NodEl
 	
	NodEl = NodElOld + 1
	
	deallocate(Elem)
	allocate(Elem(NElem*NodEl))

	do i = 1 , NElem
		ip = (i-1)*NodEl
		ipOld = (i-1)*NodElOld
		
		do j = 1,NodElOld
			ipp = ip + j
			ippOld = ipOld + j
			Elem(ipp) = ElemOld(ippOld) 
		end do
		Elem(ip + NodEl) = idGlobal
	end do	
	
end subroutine



subroutine setGlobalNodeInFace(idGlobal)
	integer , intent(in) :: idGlobal
	integer :: i , j, ip, ipp, ipOld, ippOld, NodFaceOld
	integer :: FaceOld(NFace*NodFace)
	
	FaceOld = Face
	
	NodFaceOld = NodFace
	
	NodFace = NodFaceOld + 1
	
	deallocate(Face)
	allocate(Face(NFace*NodFace))

	do i = 1 , NFace
		ip = (i-1)*NodFace
		ipOld = (i-1)*NodFaceOld
		
		do j = 1,NodFaceOld
			ipp = ip + j
			ippOld = ipOld + j
			Face(ipp) = FaceOld(ippOld) 
		end do
		Face(ip + NodFace) = idGlobal
	end do	
	
end subroutine

subroutine fixNode(idNode,dofs)
	integer , intent(in) :: idNode, dofs(:)
	integer :: iSize, i
	
	iSize = ubound(dofs,1)
	
	write(*,*) iSize
	
	do i = 1 , iSize
		if(dofs(i) <= Ndof) IDD((idNode-1)*Ndof + dofs(i) ) = -1
	end do
	
end subroutine

subroutine setGlobalNode()
	integer :: i , j, ip, ipp, ipOld, ippOld, NNodesOld
	integer :: IDDOld(NNodes*Ndof)
	real*8 :: XOld(NNodes*Ndim), DValuesOld(NNodes*Ndof)
	integer :: FixedDofs(12) 	
	
!~ 	!! with tangent homogenisation
	FixedDofs = (/ 1,2,5,6,7,8,9,10,11,12,15,16 /)
	call fixNode(1, FixedDofs)
!~ 	call fixNode(2, FixedDofs)
!~ 	call fixNode(3, FixedDofs)
!~ 	call fixNode(4, FixedDofs)

!! without tangent homogenisation
!~ 	FixedDofs(1) = 1 ; FixedDofs(2) = 2
!~ 	call fixNode(1, FixedDofs(1:2))

	
	
	Xold = X
	IDDOld = IDD
	DValuesOld = DValues

	NnodesOld = Nnodes
 	
	Nnodes = NnodesOld + 2

	deallocate(X,IDD,DValues)
	allocate(X(NNodes*Ndim),IDD(NNodes*Ndof),DValues(NNodes*Ndof))
	
	X(1:NnodesOld*Ndim) = XOld
	X(NnodesOld*Ndim + 1: Nnodes*Ndim) = 9999.0d0

	IDD(1:NnodesOld*Ndof) = IDDOld 
	IDD(NnodesOld*Ndof + 1: Nnodes*Ndof) = 0

	DValues(1:NnodesOld*Ndof) = DValuesOld 
	DValues(NnodesOld*Ndof + 1: Nnodes*Ndof) = 0

	call setGlobalNodeInFace(NnodesOld + 1)
	call setGlobalNodeInVolume(NnodesOld + 2)
	
end subroutine 

subroutine checkJacobian() !! just for triangles (I suppose, for quadrangles proceed without checking)
	integer :: i , j , ip , Element(NodEl)
	real*8 :: Xele(NodEl,Ndim), xG(Ndim) , VnFace(Ndim), normVnFace, JacFac
	
	write(0,*) "Verifying Elements"
	do i = 1 , Nelem
		ip = (i-1)*NodEl
		Element = Elem(ip+1:ip+NodEl)
		call getXele(Xele,X,Element,Ndim,NodEl)
		call findBaricentric(XG,X,Element,Ndim,NodEl)
		
		call crossProduct(VnFace,normVnFace,Xele(2,:)-Xele(1,:),Xele(3,:)-Xele(1,:))
		JacFac = normVnFace*dot_product(VnFace,Xele(4,:)-XG)
		if(JacFac<0.0d0) write(0,*) i, "JacFac" , JacFac

	end do
	
	write(0,*) "Verifying Faces"
	do i = 1 , NFace
		ip = (i-1)*NodFace
		Element = Face(ip+1:ip+NodEl)
		call getXele(Xele,X,Element,Ndim,NodFace)
		call findBaricentric(XG,X,Element,Ndim,NodFace)
		
		call crossProduct(VnFace,normVnFace,Xele(2,:)-Xele(1,:),Xele(3,:)-Xele(1,:))
		JacFac = normVnFace*dot_product(VnFace,Xele(4,:)-XG)
		if(JacFac<0.0d0) write(0,*) i, "JacFac" , JacFac

	end do

	
end subroutine

subroutine checkJacobian_tri()
	integer :: i , j , ip , Element(NodEl)
	real*8 :: Xele(NodEl,Ndim), xG(Ndim) , VnFace(Ndim), normVnFace, JacFac
	real*8 :: ElemNew(Nelem*NodEl)
	
	if(NodEl /= 3) then
		write(*,*) NodEl , "it's not a triangle, we cannot check jacobian"
		return 
	end if
	
	write(0,*) "Verifying Elements"
	do i = 1 , Nelem
		ip = (i-1)*NodEl
		Element = Elem(ip+1:ip+NodEl)
		call getXele(Xele,X,Element,Ndim,NodEl)
		call findBaricentric(XG,X,Element,Ndim,NodEl)
		
		call crossProduct(VnFace,normVnFace,Xele(2,:)-Xele(1,:),Xele(3,:)-Xele(1,:))
	
		if(VnFace(3)>0.0d0) then 
			ElemNew(ip+1:ip+NodEl) = Element
		else
			ElemNew(ip+1) = Element(1)
			ElemNew(ip+2) = Element(3)
			ElemNew(ip+3) = Element(2)
			write(0,*) i, "JacFac" , JacFac, VnFace(3)
		end if
	end do
	
	deallocate(Elem)
	allocate(Elem(NodEl*Nelem))
	
	Elem = ElemNew
	
end subroutine

subroutine TetraLin2TetraQuad(iShift)
	integer , intent(in) :: iShift
	integer :: i , j , ip , ElemSol(NodEl)
	
	do i = 1 , TimeSteps
		do j = 1 , NElem
			ip = (j-1)*NodEl
			ElemSol = (Elem(ip+1:ip+NodEl)-1)*NDof
			
			Sol(ElemSol(5) + iShift + 1 , i) = 0.5d0 * ( Sol(ElemSol(1) + iShift + 1, i) + Sol(ElemSol(2) + iShift + 1, i) ) 
			Sol(ElemSol(6) + iShift + 1 , i) = 0.5d0 * ( Sol(ElemSol(2) + iShift + 1, i) + Sol(ElemSol(3) + iShift + 1, i) ) 
			Sol(ElemSol(7) + iShift + 1 , i) = 0.5d0 * ( Sol(ElemSol(1) + iShift + 1, i) + Sol(ElemSol(3) + iShift + 1, i) ) 
			Sol(ElemSol(8) + iShift + 1 , i) = 0.5d0 * ( Sol(ElemSol(1) + iShift + 1, i) + Sol(ElemSol(4) + iShift + 1, i) ) 
			Sol(ElemSol(9) + iShift + 1 , i) = 0.5d0 * ( Sol(ElemSol(2) + iShift + 1, i) + Sol(ElemSol(4) + iShift + 1, i) ) 
			Sol(ElemSol(10) + iShift + 1, i) = 0.5d0 * ( Sol(ElemSol(3) + iShift + 1, i) + Sol(ElemSol(4) + iShift + 1, i) )
			 
		end do
	end do
	
end subroutine

subroutine doubledMesh2Single()
	integer :: i , j , ip , ipp, ipOld, ippOld, NodElOld , ElemOld(NElem*NodEl)

	NodElOld = NodEl
	ElemOld = Elem
	
	NodEl = NodElOld/2
	
	deallocate(Elem)
	allocate(Elem(NElem*NodEl))
	
	do i = 1 , NElem
		ip = (i-1)*NodEl
		ipOld = (i-1)*NodElOld
		
		Elem(ip+1:ip+NodEl) = ElemOld(ipOld+1:ipOld+NodEl)
	end do
	
end subroutine

subroutine setPeriodicMesh() !!! duplicates the mesh, being the first NodEl nodes geometrics and the others algebrics
	integer :: mapping(NNodes), i , j , k , ip, ipp, ipOld, ippOld, NodElOld, NodFaceOld , NodEdgeOld
	integer :: ElemOld(NElem*NodEl) , FaceOld(NFace*NodFace), EdgeOld(NEdge*NodEdge) 
	integer , parameter :: maxAdjust = 5
	integer :: FixedDofs(12)
	
	mapping = (/(i, i=1,NNodes)/)
	
	if(NPeriodic>0) then
		do i = 1, NPeriodic
			mapping(Periodic(NodPer,i)) = Periodic(1,i)
		end do
		
		do k = 1 , maxAdjust
		do i = 1 , NNodes
		do j = 1 , NPeriodic
			if(mapping(i) == Periodic(NodPer,j)) mapping(i) = Periodic(1,j)
		end do
		end do
		end do
	end if
	
	ElemOld = Elem
	FaceOld = Face
	EdgeOld = Edge
	
	NodElOld = NodEl
	NodFaceOld = NodFace
	NodEdgeOld = NodEdge
	
	NodEl = 2*NodEl
	NodFace = 2*NodFace
	NodEdge = 2*NodEdge
	NPeriodic = 0
	
	deallocate(Elem,Face,Edge)
	if(allocated(Periodic)) deallocate(Periodic)
	allocate(Elem(NElem*NodEl),Face(NFace*NodFace),Edge(NEdge*NodEdge))
	
	do i = 1 , Nelem
		ip = (i-1)*NodEl
		ipOld = (i-1)*NodElOld
		
		do j = 1,NodEl
			ipp = ip + j
			if(j > NodElOld) then
				ippOld = ipOld + j - NodElOld
				Elem(ipp) = mapping(ElemOld(ippOld)) 
			else
				ippOld = ipOld + j
				Elem(ipp) = ElemOld(ippOld)
			end if
		end do
	end do

	do i = 1 , NFace
		ip = (i-1)*NodFace
		ipOld = (i-1)*NodFaceOld
		
		do j = 1,NodFace
			ipp = ip + j
			if(j > NodFaceOld) then
				ippOld = ipOld + j - NodFaceOld
				Face(ipp) = mapping(FaceOld(ippOld)) 
			else
				ippOld = ipOld + j
				Face(ipp) = FaceOld(ippOld)
			end if
		end do
	end do
	
	do i = 1 , NEdge
		ip = (i-1)*NodEdge
		ipOld = (i-1)*NodEdgeOld
		
		do j = 1,NodEdge
			ipp = ip + j
			if(j > NodEdgeOld) then
				ippOld = ipOld + j - NodEdgeOld
				Edge(ipp) = mapping(EdgeOld(ippOld)) 
			else
				ippOld = ipOld + j
				Edge(ipp) = EdgeOld(ippOld)
			end if
		end do
	end do
	
!! with tangent homogenisation
!~ 	FixedDofs = (/ 1,2,5,6,7,8,9,10,11,12,15,16 /)
!~ 	call fixNode(1, FixedDofs)
!~ 	call fixNode(2, FixedDofs)
!~ 	call fixNode(3, FixedDofs)
!~ 	call fixNode(4, FixedDofs)

!! without tangent homogenisation
	FixedDofs(1) = 1 ; FixedDofs(2) = 2
	call fixNode(1, FixedDofs(1:2))

		 
end subroutine 

subroutine findSym(SymNode, refNode,bmFaceResearch,normalAxis)
	integer, intent(in) :: refNode, bmFaceResearch, normalAxis
	integer, intent(out) :: SymNode
	integer :: i , j , ip , ipp, ippDim , kppDim, axis1, axis2 
	real*8 ,parameter :: tol = 1.0e-10
	 
	SymNode = 0
	axis1=0
	axis2=0
	
	select case(normalAxis)
		case(1)
			axis1 = 2 ; axis2 = 3
		case(2)
			axis1 = 1 ; axis2 = 3
		case(3)
			axis1 = 1 ; axis2 = 2
	end select
	
	kppDim = (refNode - 1)*Ndim
	
	do i = 1 , NFace 
		if(BMFace(i) == bmFaceResearch) then
			ip = (i-1)*NodFace
			do j = 1 , NodFace
				ipp = ip + j
				ippDim = (Face(ipp) - 1)*Ndim
				if( abs(X(kppDim+axis1) - X(ippDim+axis1))<tol .and. abs(X(kppDim+axis2) - X(ippDim+axis2))<tol) then  
					SymNode = Face(ipp)
					exit 
				end if
			end do
		end if
		if(SymNode/=0) exit
	end do

end subroutine 

subroutine numberingDirichlet(idFace,dof,uD)
	integer, intent(in) :: idFace,dof
	real*8, intent(in) :: uD
	integer i,ip,n,j
	
	n=size(BMface)
	do i=1,n
		if(BMface(i)==idFace) then
			do j=1,NodFace 
			ip=(Face((i-1)*NodFace +j) - 1)*NDof
			IDD(ip+dof)=-1
			DValues(ip+dof)=uD
			end do
		end if
	end do
	
end subroutine

subroutine numberingNeumann()
	integer i,j, k
	
	NeumannFace = 0
	NNeumannFace = 0
	k=0
	do i = 1 , NNeumann
		do j=1,NFace
			if(BMFace(j)==BMNeumannFace(i)) then
				NNeumannFace(i) = NNeumannFace(i) + 1
				k = k + 1
				NeumannFace(k) = j
			end if
		end do
	end do

end subroutine

subroutine FindInteriorNode()
	integer , parameter :: Usize = 4 
	integer  ::  i, j, k, ip, ipp, ippp
	integer :: Face3D(NFace*NodEl) , u(Usize)
 	
 	
 	write(*,*) Cell_type, NodEl
 	
	do i = 1 , Nelem
		ip = (i-1)*NodEl
		do j = 1, NFace
			ipp = (j-1)*NodFace
			if (SubtractTriFromTet(Face(ipp+1:ipp+NodFace),Elem(ip+1:ip+NodEl),NodFace,NodEl,u) ) then
				ippp = (j-1)*NodEl
				
				
				!! Convention 1 : Maintain normal v12 X v13 (standard, but negative Jacobian)
				
!~ 				Face3D( ippp + 1: ippp + 3) = Face(ipp + 1:ipp + 3 ) 
!~ 				Face3D( ippp + 4) = u(1)
!~ 				
!~ 				if(Cell_type == 24) then !! just quadratic
!~ 					Face3D( ippp + 5: ippp + 7) =  Face(ipp + 4:ipp + 6)
!~ 					Face3D( ippp + 8: ippp + 10) = u(2:Usize)
!~ 				end if 
				
				!! Convertion 2 : normal as v13 X v12 (different convention, but have positive jacobian)
				Face3D( ippp + 1) = Face(ipp + 1)
				Face3D( ippp + 2) = Face(ipp + 3)
				Face3D( ippp + 3) = Face(ipp + 2)
				Face3D( ippp + 4) = u(1)
				
				if(Cell_type == 24) then !! just quadratic
					Face3D( ippp + 5) = Face(ipp + 6)
					Face3D( ippp + 6) = Face(ipp + 5)
					Face3D( ippp + 7) = Face(ipp + 4)
					Face3D( ippp + 8) = u(1 + 1 )
					Face3D( ippp + 9) = u(1 + 3)
					Face3D( ippp + 10) = u(1 + 2)
				end if 
			
			end if
		end do
	end do
 	
	deallocate(Face)
	allocate(Face(NFace*NodEl))
	Face = Face3D
	NodFace = NodEl
 	
end subroutine

Logical function SubtractTriFromTet(v,w,n,m,u) ! U = W-V, m > n assumed
	integer, parameter :: Usize = 4 , TriSize = 3 
	integer, intent(in) :: v(n),w(m),n,m
	integer, intent(out) :: u(Usize)
	integer :: i,j,k ,TriaVert(TriSize) , InteriorVert, p1, p2
	
	u(1) = -1
	InteriorVert = -1 !! prevent warning
 	
	SubtractTriFromTet = .False.
	
	outer: do i = 1,4
				do j = 1,3
					if(w(i) == v(j)) then
						triaVert (j) = i
						cycle outer
					end if
				end do
				
				if ( u(1)==-1 ) then
					u(1) = w(i)
					InteriorVert = i
				else
					return
				end if
			end do outer
			

	SubtractTriFromTet = .True.
	if(Cell_type == 10) return !! Tetrahedra Linear
	
	do i = 1,3
		if (interiorVert < triaVert(i)) then
			p1 = interiorVert
			p2 = triaVert(i)
		else
			p2 = interiorVert
			p1 = triaVert(i)
		end if
		
		if(p1==1 .and. p2==2) then
			u(i+1) = w(5)
		elseif(p1==2 .and. p2==3) then
			u(i+1) = w(6)
		elseif(p1==1 .and. p2==3) then
			u(i+1) = w(7)
		elseif(p1==1 .and. p2==4) then
			u(i+1) = w(8)
		elseif(p1==2 .and. p2==4) then
			u(i+1) = w(9)
		elseif(p1==3 .and. p2==4) then
			u(i+1) = w(10)
		end if
	end do			

end function

!~ subroutine Single2DoubleFieldP1bP1()
!~ 	integer ::  i, j , k , ip, ipp, ipDof, ipNod, ipOld
!~ 	real*8 :: XOld(NNodes*Ndim) , DValuesOld(NNodes*Ndof)
!~ 	integer :: ElemOld(NElem*NodEl), IDDold(NNodes*Ndof)
!~ 	integer :: NNodesOld, NodElOld
!~ 	real*8 :: pAux(Ndim)
!~ 	
!~ 	write(0,*) "Modifying Mesh to P1bP1"
!~ 	
!~ 	NNodesOld = NNodes ; NodElOld = NodEl
!~ 	XOld = X ; DValuesOld = DValues; ElemOld = Elem; IDDold = IDD
!~ 	
!~ 	deallocate(X,Elem,IDD,DValues)  
!~ 	
!~ 	NodEl = NodElOld + 1
!~ 	NNodes = NNodesOld + NElem
!~ 	
!~ 	allocate(X(NNodes*Ndim),DValues(NNodes*Ndof),Elem(NElem*NodEl), IDD(NNodes*Ndof))
!~ 	
!~ 	X(1:NNodesOld*Ndim) = Xold(:)
!~ 	
!~ 	DValues(1:NNodesOld*Ndof) = DValuesOld(:)
!~ 	IDD(1:NNodesOld*Ndof) = IDDold(:)
!~ 	
!~ 	DValues(NNodesOld*Ndof + 1 :) = 0.0d0
!~ 	IDD(NNodesOld*Ndof + 1 :) = 0
!~ 	
!~ 	do i = 1,NElem
!~ 		ip = (i-1)*NodEl
!~ 		ipOld = (i-1)*NodElOld
!~ 		
!~ 		pAux = 0.0d0
!~ 		do j = 1 , NodElOld
!~ 			Elem(ip+j) = ElemOld(ipOld + j)
!~ 			
!~ 			do k = 1, Ndim
!~ 				pAux(k) = pAux(k) + Xold((ElemOld(ipOld+j)-1)*Ndim + k ) 
!~ 			end do
!~ 			
!~ 		end do
!~ 		
!~ 		pAux = pAux/NodElOld
!~ 		
!~ 		do k = 1, Ndim
!~ 			X( (NNodesOld + i - 1)*Ndim + k) = pAux(k) 
!~ 		end do
!~ 		
!~ 		Elem(ip + NodElOld + 1) = i + NNodesOld 
!~ 	
!~ 	end do
!~ 	
!~ end subroutine

subroutine convert2StandardNumeration()
	integer , parameter :: NodTetra = 4 , NodTriang = 3  
	integer :: mapTetra(NodEl - NodTetra) , mapTriang(NodFace - NodTriang) ,&
				TetraAux(NodEl - NodTetra) , triAux(NodFace - NodTriang)
	integer :: swap1, swap2, i , ip , idiff , j

! TETGEN Standard
!~ 	mapTetra(1) = 7  ! 5
!~ 	mapTetra(2) = 8 ! 6 
!~ 	mapTetra(3) = 10! 7 
!~ 	mapTetra(4) = 6 ! 8
!~ 	mapTetra(5) = 9 ! 9
!~ 	mapTetra(6) = 5 ! 10
!~ 	
!~ 	mapTriang(1) = 6   
!~ 	mapTriang(2) = 4  
!~ 	mapTriang(3) = 5   
!~ 	
	mapTetra(1) = 5  ! 5
	mapTetra(2) = 6 ! 6 
	mapTetra(3) = 7 ! 7 
	mapTetra(4) = 8 ! 8
	mapTetra(5) = 10 ! 9
	mapTetra(6) = 9 ! 10	

	mapTriang(1) = 4   
	mapTriang(2) = 5  
	mapTriang(3) = 6   
	
	idiff = NodEl - NodTetra

	write(0,*) size(Elem) , Nelem, NodEl

	do i = 1 , Nelem
		ip = (i-1)*NodEl + NodTetra
		TetraAux = Elem(ip+1:ip+idiff)
		do j = 1, idiff
			Elem(ip+j) = TetraAux(mapTetra(j)-NodTetra)
		end do
	end do
	
	idiff = NodFace - NodTriang
	
	write(0,*) size(Face) , NFace, NodFace, size(TriAux)
	
	do i = 1 , NFace
		ip = (i-1)*NodFace + NodTriang
		TriAux = Face(ip+1:ip+idiff)
		do j = 1, idiff
			Face(ip+j) = TriAux(mapTriang(j)-NodTriang)
		end do
	end do
	
	
end subroutine

subroutine selectCellType()
	
	Cell_Type = 0
	
	select case(NodEL)
		case(2) 
			Cell_type = 3 ! Line
		case(3)
			Cell_type = 5 ! triangle
		case(4)
			if(NdimE == 3) Cell_type = 10 ! Tetahedron
			if(NdimE == 2) Cell_type = 9 ! quad
		case(5)
			Cell_type = 10 ! Tetahedron ! with bubble
		case(6) 
			Cell_type = 22 ! Quadratic Triangle
		case(8) 
			Cell_type = 12 ! Hexahedron
		case(9)	
			Cell_type = 23 ! quadratic quad
		case(10)
			Cell_type = 24 ! Quadratic Tetrahedra
		case default
			Cell_type = 0
			write(*,*) "Element unknown"
	end select
	
	select case(Cell_type)
		case(3) ! linear line
			NodEdge = 1
			NodFace = 1 ! The face is the point
			NGP = 1
			NodEl_vtk = 2
		case(5) ! triangle linear
			NodEdge = 2
			NodFace = 2 ! The face is the Edge
			NGP = 1
			NodEl_vtk = 3
		case(9) ! Quad
			NodEdge = 2
			NodFace = 2 ! The face is the Edge
			NGP = 1 
			NodEl_vtk = 4
		case(10) ! Tetrahedra
			NodEdge = 2
			NodFace = 3 ! No face
!~ 			if(FindInterior) NodFace = 4
			NGP = 1
			NodEl_vtk = 4
		case(12) ! Linear Hexahedra
			NodEdge = 2
			NodFace = 4 ! No face
!~ 			if(FindInterior) NodFace = 8
			NGP = 1 ! I don't know, but not useful at moment
			NodEl_vtk = 8
		case(22) ! Quadratic Triangle
			NodEdge = 3
			NodFace = 3 ! The face is the Edge
			NGP = 1 
			NodEl_vtk = 6
		case(23) ! Quadratic Quad
			NodEdge = 3
			NodFace = 3 ! The face is the Edge
			NGP = 1
			NodEl_vtk = 8
		case(24) ! Quadratic Tetrahedra
			NodEdge = 3
			NodFace = 6 ! No face
!~ 			if(FindInterior) NodFace = 10
			NGP = 1 !!! Considering damage in the middle
!~ 			NGP = 4
			NodEl_vtk = 10
		case default
			write(0,*) "You Must specify a valid element"
	end select
	
end subroutine

subroutine findBaricentric(XG,Xnodes,Element,Ndim,NodEl)
	integer, intent(in) :: Ndim , NodEl, Element(:)
	real(8), intent(in) :: Xnodes(:)
	real(8), intent(out) :: XG(Ndim)
	real(8) :: Xele(NodEl,Ndim)	
	integer :: i
	
	call getXele(Xele,Xnodes,Element,Ndim,NodEl)
	do i = 1,Ndim
		xG(i) = sum(Xele(:,i))/real(NodEl)
	end do
	
end subroutine


subroutine getXele(Xele,Xnodes,Element,Ndim,NodEl)
	integer, intent(in) :: Ndim , NodEl, Element(:)
	real(8), intent(in) :: Xnodes(:)
	real(8), intent(out) :: Xele(NodEl,Ndim)
	integer :: i, j , ipp 
	
	Xele = 0.0d0
	
!~ 	write(*,*) NodEl
!~ 	call myPrint(Element)
	
	
	do i=1,Ndim
		do j=1,NodEl
			ipp = (Element(j) - 1)*Ndim 
			Xele(j,i) = Xnodes(ipp+i)
		end do
	end do
	
!~ 	pause
	
end subroutine

subroutine getArrayXele(Xele,e)
	integer, intent(in) :: e
	real(8), intent(out) :: Xele(NodEl*Ndim)
	integer :: i, j , ipp , ip, iShiftEl
	
	Xele = 0.0d0
	iShiftEl = (e-1)*NodEl	

!~ 	write(*,*) e, NodEl
!~ 	call myPrint(Elem(1:NodEl))
	
	do i=1,NodEl
		ip = (i-1)*Ndim
		ipp = (Elem(iShiftEl + i) - 1)*Ndim
		do j=1,Ndim
			Xele(ip+j) = X(ipp+j)
		end do
	end do
	
end subroutine

subroutine getArraySolEl(SolEl,e,timeStep)
	integer, intent(in) :: e, timeStep
	real(8), intent(out) :: SolEl(NodEl*Ndof)
	integer :: i, j , ip, ipp , iShiftEl
	
	SolEl = 0.0d0
	iShiftEl = (e-1)*NodEl
	
	do i=1,NodEl
		ip = (i-1)*Ndof
		ipp = (Elem(iShiftEl + i) - 1)*Ndof
		do j=1,Ndof
			SolEl(ip+j) = Sol(ipp+j, timeStep)
		end do
	end do
	
end subroutine


subroutine setPeriodicConditions(PeriodicParam,NPeriodicParam) !! Create list of periodic correspondences (not informed in gmsh)
	integer , intent(in) :: NPeriodicParam , PeriodicParam(NPeriodicParam,NodPer + 1)
	integer :: i , j , k , kp, jp, ip, ipp , cont, flagAux(NNodes,NPeriodicParam) , flagAux2(NNodes), SymNode , axis1 , axis2, flag
	real*8 ,parameter :: tol = 1.0e-10
	
	Nperiodic = 0
	flagAux = 0
	
	do i = 1 , NPeriodicParam
	
		do j = 1 , NFace 
			flag = 0 
			if(BMFace(j) == PeriodicParam(i,1)) then
				flag = 1
			else if(BMFace(j) == PeriodicParam(i,2)) then
				flag = 2
			end if
			
			if(flag /=0 ) then
				do k = 1 , NodFace
					flagAux(Face((j-1)*NodFace + k),i) = flag
				end do
			end if		
		end do
		
		Nperiodic = Nperiodic + sum(flagAux(:,i))/3 !!! Because there is a flag 2
		
	end do
		
	allocate(Periodic(NodPer,Nperiodic))
	
	cont = 0
	do i = 1 , NPeriodicParam 
		select case(PeriodicParam(i,3))
			case(1)
				axis1 = 2 ; axis2 = 3
			case(2)
				axis1 = 1 ; axis2 = 3
			case(3)
				axis1 = 1 ; axis2 = 2
			case default
				axis1 = 0 ; axis2 = 0
		end select

		do j = 1, NNodes
			if(flagAux(j,i) == 1) then
				jp = (j-1)*Ndim 
				cont = cont+1
				do k = 1, NNodes
					if(flagAux(k,i) == 2) then 
						kp = (k-1)*Ndim
						 
						if( abs(X(jp+axis1) - X(kp+axis1))<tol .and. abs(X(jp+axis2) - X(kp+axis2))<tol) then  
							Periodic(1,cont) = j
							Periodic(2,cont) = k
						end if
					end if
				end do
			end if
		end do
		
	end do

end subroutine 

subroutine setPeriodicConditionsAuto() !! Create list of periodic correspondences direct from readFromGMSH 
	integer :: iAux, iAux2, iAux3, i , j
	integer , parameter :: IUnitMesh = 12 
	integer :: FixedDofs(12) 	
	
	write(*,*) "setPeriodicConditionsAuto"
	
	open (IUnitMesh, file='mesh.msh', status= 'old' ) 
	call FindKeyword(IUnitMesh,"$Periodic")
	
	Nperiodic = 0
	read(IUnitMesh,*) iAux ! number of elements groups of periodic conditions
	!! set the total number of periodic conditions
	do i = 1 , iAux
		read(IUnitMesh,*) 
		read(IUnitMesh,*) iAux2
		write(*,*) "first loop , ", iAux2
		Nperiodic = Nperiodic + iAux2
		do j = 1,iAux2
			read(IUnitMesh,*) 
		end do 
	end do
	
!~ 	write(*,*) NodPer, Nperiodic
	
	allocate(Periodic(NodPer,Nperiodic))
	
	rewind(IUnitMesh) 
	call FindKeyword(IUnitMesh,"$Periodic")
	read(IUnitMesh,*)
	iAux3 = 1 !! pointer to the next periodic assigned
	do i = 1 , iAux
		read(IUnitMesh,*) 
		read(IUnitMesh,*) iAux2
		write(*,*) "second loop , ", iAux2
		do j = 1,iAux2
			read(IUnitMesh,*) Periodic(:,iAux3)
			iAux3 = iAux3 + 1 
		end do 
	end do
		
!! with tangent homogenisation
	FixedDofs = (/ 1,2,5,6,7,8,9,10,11,12,15,16 /)
!~ 	call fixNode(1, FixedDofs)
!~ 	call fixNode(2, FixedDofs)
!~ 	call fixNode(3, FixedDofs)
!~ 	call fixNode(4, FixedDofs)

!! without tangent homogenisation
!~ 	FixedDofs(1) = 1 ; FixedDofs(2) = 2
!~ 	call fixNode(1, FixedDofs(1:2))

	
	close(IUnitMesh)

end subroutine 


subroutine modifyXFEM() !
	integer :: i, j, k, ip, jp, kp, kcut, kpp, kp1, kp2, kpp4p, kpp4m, kpp5p, kpp5m, kppp4p, kppp4m, kppp5p, kppp5m 
	integer :: IshiftEnr, IshiftUe, step, countCut, NelemNew, NnodesNew, NdofNew
	integer :: ElemLoc(NodEl)
	real*8 :: g(NodEl), t1, t2 , Xele(NodEl,Ndim), x4(Ndim) , x5(Ndim), h
	real*8, allocatable :: XNew(:) , SolNew(:)
	integer, allocatable :: ElemNew(:) 
	logical :: isCuttedEl(Nelem), isScalar
	
	isScalar =.true.
!~ 	
	if(isScalar) then
		IshiftEnr = 2
		NdofNew = 4
		IshiftUe = 1
	else
		IshiftEnr = 4
		NdofNew = 7
		IshiftUe = 2
	end if

	step = 2
	countCut = 0
	isCuttedEl = .False.
	kcut = 0
	
	!! choose elements to be cutted
	do i = 1 , Nelem
		ip = (i-1)*NodEl
		
		ElemLoc = Elem(ip + 1:ip + NodEl)
		
		do j = 1 , NodEl
			jp = (ElemLoc(j)-1)*NDof + iShiftEnr + 1
			g(j) = Sol(jp,step) 
		end do 
		
		if(maxval(g)*minval(g)<0.0d0) then  
			isCuttedEl(i) = .True.
			
			! reduces to the case g(1)*g(2)>0
			if(g(2)*g(3) > 0.0d0) then 
				Elem(ip + 1) = ElemLoc(2)
				Elem(ip + 2) = ElemLoc(3)
				Elem(ip + 3) = ElemLoc(1)	
			elseif(g(3)*g(1) > 0.0d0) then
				Elem(ip + 1) = ElemLoc(3)
				Elem(ip + 2) = ElemLoc(1)
				Elem(ip + 3) = ElemLoc(2)
			end if
			
		end if
		
	end do

	countCut = count(isCuttedEl)
	NelemNew = Nelem + 2*countCut
	NnodesNew = Nnodes + 4*countCut
	
	allocate(ElemNew(NelemNew*NodEl), XNew(NnodesNew*Ndim) , SolNew(NnodesNew*NdofNew))
	SolNew = 0.0d0

	SolNew = 0.0d0
	XNew(1:Nnodes*Ndim) = X
	ElemNew(1:Nelem*NodEl) = Elem
	
	do i = 1, Ndof
		SolNew(i:Nnodes*NdofNew:NdofNew) = Sol(i:Nnodes*Ndof:Ndof,step)
	end do
	
	do i = 1 , Nelem
		ip = (i-1)*NodEl
		ElemLoc = Elem(ip + 1:ip + NodEl)
		
		if(isCuttedEl(i)) then
			kcut = kcut + 1
			kp1 = (Nelem + 2*(kcut - 1))*NodEl !! point to the first added element
			kp2 = kp1 + NodEl !! point to the second added element
			kpp = Nnodes + 4*(kcut-1)  !number of the to the first added node minus 1
			kpp4p = (Nnodes + 4*(kcut-1))*Ndim ! point to the first added node (dim)
			kpp4m = kpp4p + Ndim! point to the second added node (dim)
			kpp5p = kpp4p + 2*Ndim !point to the third added node (dim)
			kpp5m = kpp4p + 3*Ndim !point to the fouth added node (dim)
			kppp4p = (Nnodes + 4*(kcut-1))*NdofNew !point to the first added node (dof)
			kppp4m = kpp4p + NdofNew !point to the second added node (dof)
			kppp5p = kpp4p + 2*NdofNew !point to the third added node (dof)
			kppp5m = kpp4p + 3*NdofNew !point to the fouth added node (dof)
			
			do j = 1 , NodEl
				jp = (ElemLoc(j)-1)*NDof + iShiftEnr + 1
				g(j) = Sol(jp,step) 
			end do
			
			call getXele(Xele,X,ElemLoc,Ndim,NodEl)
			
			t1 = dabs(g(2))/( dabs(g(2)) + dabs(g(3)) ) !! edge 2-3
			t2 = dabs(g(3))/( dabs(g(3)) + dabs(g(1)) ) !! edge 3-1
	 
			x4 = (1.0 - t1)*Xele(2,:) + t1*Xele(3,:)
			x5 = (1.0 - t2)*Xele(3,:) + t2*Xele(1,:)
			
			!~ triangleFlag = (/ 3, 5, 4, 1, 4, 5, 1, 2, 4 /)
			ElemNew(ip + 1) = ElemLoc(3)
			ElemNew(ip + 2) = kpp + 3 ! 5+
			ElemNew(ip + 3) = kpp + 1 ! 4+
			ElemNew(kp1 + 1) = ElemLoc(1) 
			ElemNew(kp1 + 2) = kpp + 2 ! 4-
			ElemNew(kp1 + 3) = kpp + 4 ! 5-
			ElemNew(kp2 + 1) = ElemLoc(1)
			ElemNew(kp2 + 2) = ElemLoc(2)
			ElemNew(kp2 + 3) = kpp + 2 ! 4-

			
			do j = 1 , Ndof !! just for the new nodes
				SolNew(kppp4p + j) = (1.0 - t1)*Sol((ElemLoc(2)-1)*Ndof + j,step) + t1*Sol((ElemLoc(3)-1)*Ndof + j,step)      !4+
				SolNew(kppp4m ) = SolNew(kppp4p + j)     !4- 
				SolNew(kppp5p + j) = (1.0 - t2)*Sol((ElemLoc(3)-1)*Ndof + j, step) + t2*Sol((ElemLoc(1)-1)*Ndof + j, step)    ! 5+
				SolNew(kppp5m + j) = SolNew(kppp5p + j)	! 5-
			end do
		
			h = dsign(1.0d0,g(3))
			do j = Ndof + 1 , NdofNew 
				SolNew(kppp4p + j) = SolNew(kppp4p + j - NdofNew) + h*SolNew(kppp4p + j - NdofNew + iShiftUe)               !4+
				SolNew(kppp4m + j) = SolNew(kppp4m + j - NdofNew) - h*SolNew(kppp4m + j - NdofNew + iShiftUe)  				!4- 
				SolNew(kppp5p + j) = SolNew(kppp5p + j - NdofNew) + h*SolNew(kppp5p  + j - NdofNew + iShiftUe)  ! 5+
				SolNew(kppp5m + j) = SolNew(kppp5m + j - NdofNew) - h*SolNew(kppp5m + j - NdofNew + iShiftUe) 	! 5-
				SolNew((ElemLoc(1)-1)*NdofNew + j) = Sol((ElemLoc(1)-1)*Ndof + j - Ndof, step) &
													- h*Sol((ElemLoc(1)-1)*Ndof + j - Ndof + iShiftUe, step)
				SolNew((ElemLoc(2)-1)*NdofNew + j) = Sol((ElemLoc(2)-1)*Ndof + j - Ndof, step) &
													- h*Sol((ElemLoc(2)-1)*Ndof + j - Ndof + iShiftUe, step)
				SolNew((ElemLoc(3)-1)*NdofNew + j) = Sol((ElemLoc(3)-1)*Ndof + j - Ndof, step) &
													+ h*Sol((ElemLoc(3)-1)*Ndof + j - Ndof + iShiftUe, step)
			end do
			
						
			XNew(kpp4p + 1 : kpp4p + Ndim) = x4
			XNew(kpp4m + 1 : kpp4m + Ndim) = x4
			XNew(kpp5p + 1 : kpp5p + Ndim) = x5
			XNew(kpp5m + 1 : kpp5m + Ndim) = x5
				
		else
			
			do j = 1 , NodEl
				if( .not. dabs(Sol((ElemLoc(j)-1)*Ndof + IshiftEnr, step)) > 0.0d0 ) then
					do k = Ndof + 1, NdofNew
						SolNew((ElemLoc(j)-1)*NdofNew + k) = Sol((ElemLoc(j)-1)*Ndof + k-Ndof, step)
					end do
				end if
			end do  
		
		end if
					
	end do
	
	deallocate(Elem,X,Sol)
	allocate(Elem(NelemNew*NodEl), X(NnodesNew*Ndim) , Sol(NnodesNew*NdofNew,2))
	Elem = ElemNew
	X = XNew
	Sol = 0.0d0
	Sol(:,2) = SolNew
	
	
	Nelem = NelemNew
	Nnodes = NnodesNew
	Ndof = NdofNew
	 
	deallocate(ElemGroup)
	allocate(ElemGroup(NelemNew))
	ElemGroup = 1
	
end subroutine

subroutine FindInteriorNode2D()
	integer , parameter :: Usize = 3 
	integer  ::  i, j, k, ip, ipp, ippp, kp1,kp2, kp3
	integer :: Face2D(NFace*NodEl) , u(Usize)
	real*8 :: vecU(Ndim), vecV(Ndim), normal(Ndim), areaDummy 
 	
 	
 	write(*,*) Cell_type, NodEl,NFace
 	
	do i = 1 , Nelem
		ip = (i-1)*NodEl
		do j = 1, NFace
			ipp = (j-1)*NodFace
			if (SubtractLineFromTri(Face(ipp+1:ipp+NodFace),Elem(ip+1:ip+NodEl),NodFace,NodEl,u) ) then
				ippp = (j-1)*NodEl
				
				kp1 = (Face( ipp + 1)-1)*Ndim
				kp2 = (Face( ipp + 2)-1)*Ndim
				kp3 = (u(1)-1)*Ndim
				
				Face2D( ippp + 3) = u(1)
				
				do k = 1 , Ndim
					vecU(k) = X(kp1+k) - X(kp3+k)
					vecV(k) = X(kp2+k) - X(kp3+k)
				end do
				
				call crossProduct(normal,areaDummy,vecU,vecV)
				
				if(normal(3)>0.0) then
					Face2D( ippp + 1) = Face(ipp + 1)
					Face2D( ippp + 2) = Face(ipp + 2)
				else
					Face2D( ippp + 1) = Face(ipp + 2)
					Face2D( ippp + 2) = Face(ipp + 1)
				end if
				
			end if
		end do
	end do
 	
	deallocate(Face)
	allocate(Face(NFace*NodEl))
	Face = Face2D
	NodFace = NodEl
 	
end subroutine

subroutine FindInteriorNode2DQuad()
	integer , parameter :: Usize = 4 
	integer  ::  i, j, k, ip, ipp, ippp, kp1,kp2, kp3
	integer :: Face2D(NFace*NodEl) , u(Usize)
	real*8 :: vecU(Ndim), vecV(Ndim), normal(Ndim), areaDummy 
 	
 	write(*,*) Cell_type, NodEl,NFace, NodFace
 	
	do i = 1 , Nelem
		ip = (i-1)*NodEl
		do j = 1, NFace
			ipp = (j-1)*NodFace
			if (SubtractLineFromQuad(Face(ipp+1:ipp+NodFace),Elem(ip+1:ip+NodEl),NodFace,NodEl,u) ) then
				ippp = (j-1)*NodEl
				
				kp1 = (Face( ipp + 1)-1)*Ndim
				kp2 = (Face( ipp + 2)-1)*Ndim
				kp3 = (u(1)-1)*Ndim
				
				Face2D( ippp + 1) = Face(ipp + 1)
				Face2D( ippp + 2) = Face(ipp + 2)
				Face2D( ippp + 3) = u(1)
				Face2D( ippp + 4) = u(2)
				
			end if
		end do
	end do
 	
	deallocate(Face)
	allocate(Face(NFace*NodEl))
	Face = Face2D
	NodFace = NodEl
 	
end subroutine

Logical function SubtractLineFromTri(v,w,n,m,u) ! U = W-V , m > n assumed
	integer, parameter ::  LineSize = 2 , 	Usize = 3
	integer, intent(in) :: v(n),w(m),n,m
	integer, intent(out) :: u(Usize)
	integer :: i,j,k ,LineVert(LineSize) , InteriorVert, p1, p2

	u(1) = -1
	InteriorVert = -1 !! prevent warning
 
	SubtractLineFromTri = .False.
	
	outer: do i = 1,NodEl
				do j = 1, LineSize
					if(w(i) == v(j)) then
						LineVert (j) = i
						cycle outer
					end if
				end do
				
				if ( u(1)==-1 ) then
					u(1) = w(i)
					InteriorVert = i
				else
					return
				end if
			end do outer
			

	SubtractLineFromTri = .True.
	
end function

Logical function SubtractLineFromQuad(v,w,n,m,u) ! U = W-V , m > n assumed
	integer, parameter ::  Usize = 2
	integer, intent(in) :: v(n),w(m),n,m
	integer, intent(out) :: u(Usize)
	integer :: i
	
	if(isVinW(v,w,n,m)) then 
		SubtractLineFromQuad = .True.
		do i = 1, m 
			if(w(i) == v(1)) then
				u(1) = w(mod(i+1,NodEl)+1)
				u(2) = w(mod(i+2,NodEl)+1) 
			end if
		end do	
	else
		SubtractLineFromQuad = .False.
	end if	
	
end function

!~ subroutine generateEdgesFromTriangles(Edges,Elem,Nelem)


end module 

