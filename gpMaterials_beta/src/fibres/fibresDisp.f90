Subroutine fibresDisp(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
	use globalVariables, only : getF, NdimE, getMaterial, getDamagePar, maxMatPar, maxDamagePar
	use fibresLib
	use damageFibLib
	
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT, NodElt! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

!~ 	!   =====   END ARGUMENTS  =======
         
	Integer :: A, B, p, q, ApRow, BpCol, App, Bqp, ApDim
	Real*8 :: F(NdimE,NdimE), Uaux(NdimE)
	Real*8 ::  tfib(NdimE), signal
	real*8 , allocatable ::  Xel(:) !! all have dimension NodG*NdimE
	integer :: NodG, iFtype
	integer :: iShiftU  

	iShiftU = nint(commonPar(1))
	iFtype = nint(commonPar(2)) 
	
	NodG = NodElt - 1	
	
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	
	call getF(F,iFtype)

	Do A=1, NodG !LoopRow
		ApRow  = (A-1)*iDofT + iShiftU
		ApDim = (A-1)*NdimE
		
		Uaux = matmul(F,Xel(ApDim + 1: ApDim + NdimE)) - Xel(ApDim + 1: ApDim + NdimE)
		
		Do p = 1 , NdimE
			App = ApRow+p
			
			AE(App,App) = 1.0d0 
			BE(App) = Uaux(p) 
			
		end do ! loop Ap
	Enddo !LoopRow
	
	
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine fibresDispS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
    
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
    
    Integer :: A,B,q,p, ApRow, BpCol, iShiftU
    integer , parameter :: NdimE = 2, NodG = 2

    iShiftU = nint(commonPar(1))
	
	do A = 1,NodG
        ApRow = (A-1) * iDofT + iShiftU
		
		do p=1,nDimE
			Coupling (ApRow + p, ApRow + p) = 1
		enddo
	enddo


end Subroutine
