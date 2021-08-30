Subroutine linearBarElementMS(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
	use globalVariables, only : getF, NdimE, getMaterial, maxMatPar
	use fibresLib
	use ptsGaussLib2, only : setNodG
	
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT, NodElt! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

!~ 	!   =====   END ARGUMENTS  =======
         
	Integer :: A, B, p, q, ApRow, BpCol, App, Bqp
	Real*8 :: MatPar(maxMatPar)
	Real*8 ::  Lfib, E, Area , afib(NdimE) , r, c2, s2, sc, F(NdimE,NdimE)
	real*8 , allocatable :: Xel(:) !! all have dimension NodG*NdimE
	integer :: NodG,iShiftUf, constLaw , iFemType,iFType, iMaterial

	iShiftUf = nint(CommonPar(1))
	iFemType = nint(commonPar(2))
	iFType = nint(commonPar(3))
	iMaterial = nint(commonPar(4))
!~ 	iDamageParType = nint(commonPar(5))
	    
    call setNodG(iFemtype, NodG)
    		
	call getMaterial(constLaw, matPar, iMaterial)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	
	call getF(F,iFtype) !!! with identity summed

	afib(1) = Xel(3) - Xel(1)
	afib(2) = Xel(4) - Xel(2)  
	
	Lfib = norm2(afib)
	
	afib = afib/Lfib
	
	c2 = afib(1)**2.0
	s2 = afib(2)**2.0
	sc = afib(1)*afib(2)
	
	E = matPar(1)
	Area = matPar(2)
	r = E*Area/Lfib
		
	AE(iShiftUf + 1,iShiftUf + 1) = c2
	AE(iShiftUf + 1,iShiftUf + 2) = sc
	AE(iShiftUf + 1,iShiftUf + 3) = -c2
	AE(iShiftUf + 1,iShiftUf + 4) = -sc

	AE(iShiftUf + 2,iShiftUf + 1) = sc
	AE(iShiftUf + 2,iShiftUf + 2) = s2
	AE(iShiftUf + 2,iShiftUf + 3) = -sc
	AE(iShiftUf + 2,iShiftUf + 4) = -s2

	AE(iShiftUf + 3,iShiftUf + 1) = -c2
	AE(iShiftUf + 3,iShiftUf + 2) = -sc
	AE(iShiftUf + 3,iShiftUf + 3) = c2
	AE(iShiftUf + 3,iShiftUf + 4) = sc
	
	AE(iShiftUf + 4,iShiftUf + 1) = -sc
	AE(iShiftUf + 4,iShiftUf + 2) = -s2
	AE(iShiftUf + 4,iShiftUf + 3) = sc
	AE(iShiftUf + 4,iShiftUf + 4) = s2
	
	AE = r*AE
	
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine linearBarElementMSS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
    
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
    
    Integer :: A,B,q,p, ApRow, BpCol,iShiftUf
    integer , parameter :: NdimE = 2, NodG = 2

    iShiftUf = nint(commonPar(1))

	Coupling = 1 !!! Wrong


end Subroutine
