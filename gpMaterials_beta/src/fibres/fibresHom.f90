Subroutine fibresHom(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
	use globalVariables, only : getF, NdimE, getMaterial, maxMatPar, addToPKhomGen
	use fibresLib
	use DETERMINANT , only : iElem
	
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT, NodElt! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

!~ 	!   =====   END ARGUMENTS  =======
         
	Real*8 :: PK(NdimE,NdimE), L, normal(NdimE), dVol, t1(NdimE), t2(NdimE),y1(NdimE), y2(NdimE)
	real*8 , allocatable ::  SolTraction(:) , Xel(:) !! all have dimension NodG*NdimE
	integer :: NodG, Nelem, iElemBegin
	integer :: iShiftTraction

	iShiftTraction = nint(commonPar(1)) 
	iElemBegin = nint(commonPar(2))
	Nelem = nint(commonPar(3))
	
!~ 	NodG = NodElt - 1 	
	NodG = 2
		
	call getSliceAllocate(SolTraction,Sol1,1,NodG,iShiftTraction + 1 ,iShiftTraction + NdimE, iDofT)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	
	normal(1) = Xel(4) - Xel(2) 
	normal(2) = -(Xel(3) - Xel(1)) 
	
	L = norm2(normal)
	normal = normal/L
 	
	t1 = SolTraction(1:2)
	t2 = SolTraction(3:4)
	y1 = Xel(1:2)
	y2 = Xel(3:4)

	PK(1,1) = 0.5d0*(t1(1)*y1(1) + t2(1)*y2(1))
	PK(1,2) = 0.5d0*(t1(1)*y1(2) + t2(1)*y2(2))
	PK(2,1) = 0.5d0*(t1(2)*y1(1) + t2(2)*y2(1))
	PK(2,2) = 0.5d0*(t1(2)*y1(2) + t2(2)*y2(2))	
	
	dVol  = L
		
!~ 	write(0,*) iElem, Nelem 
	call addToPKhomGen(PK,dVol,iElem,iElemBegin,Nelem)
	
	AE(1,1) = 1.0d0
	BE(1) = Sol1(1)
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine fibresHomS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
    
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)

	Coupling (1, 1) = 1

end Subroutine
