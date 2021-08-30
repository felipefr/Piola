Subroutine fibresHom2(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
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
         
	Real*8 :: PK(NdimE,NdimE), Lfib, afib(NdimE), dVol
	real*8 , allocatable ::  SolTraction(:) , Xel(:) !! all have dimension NodG*NdimE
	integer :: NodG, Nelem, iElemBegin
	integer :: iShiftTraction

	iShiftTraction = nint(commonPar(1)) 
	iElemBegin = nint(commonPar(2))
	Nelem = nint(commonPar(3))
	
	NodG = NodElt - 1 	
		
	call getSliceAllocate(SolTraction,Sol1,1,NodG,iShiftTraction + 1 ,iShiftTraction + NdimE, iDofT)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	
	afib(1) = Xel(3) - Xel(1)
	afib(2) = Xel(4) - Xel(2)  
	
	Lfib = norm2(afib)
	afib = afib/Lfib
	 	
	PK(1,1) = Lfib*SolTraction(3)*afib(1)
	PK(1,2) = Lfib*SolTraction(3)*afib(2)
	PK(2,1) = Lfib*SolTraction(4)*afib(1)
	PK(2,2) = Lfib*SolTraction(4)*afib(2)	
	
	dVol  = 1.0d0/real(Nelem - iElemBegin + 1) 
		
!~ 	write(0,*) iElem, Nelem 
	call addToPKhomGen(PK,dVol,iElem,iElemBegin,Nelem)
	
	AE(1,1) = 1.0d0
	BE(1) = Sol1(1)
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine fibresHom2S(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
    
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)

	Coupling (1, 1) = 1

end Subroutine
