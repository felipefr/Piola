
!     ------------------------------------------------------------------
Subroutine zeroAverageFibre(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, &
						Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
	use globalVariables, only : NdimE
	use fibresLib
		
	IMPLICIT NONE
	
	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	! Reals Vectors and matrices
	Real*8  :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), &
		Param(*), JParam(*), CommonPar(*)
	
	!   =====   END ARGUMENTS  =======
     		
	Integer :: ApRow, BpCol , App,Bpp, p, A
	integer, parameter ::  NodG = 2
	real*8 , allocatable ::  Xel(:) 
	Real*8 :: eps, pen 
	integer :: iShiftFluc, iShiftLag 
		
	iShiftFluc = nint(CommonPar(1))
	iShiftLag = nint(CommonPar(2))
	pen = commonPar(3)
	eps = commonPar(4)
	
	BpCol = iShiftLag !!! located in first node
	
	Do A=1,NodG 
		ApRow  = (A-1)*iDofT + iShiftFluc
		
		Do p = 1 , NdimE 
			App = ApRow + p
			Bpp = BpCol + p !! it's a kind of identity 
												 	
			AE(Bpp,App) = pen
			AE(App,Bpp) = pen
				
			AE(Bpp,Bpp) = eps
				
			BE(Bpp) = eps*Sol1(Bpp)

		end do ! loop Ap
	Enddo !LoopRow

end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine zeroAverageFibreS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	use globalVariables, only : NdimE
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)
	integer ::  A, p , App, ApRow, BpCol, Bpp, iShiftFluc, iShiftLag
	integer , parameter :: NodG = 2
	
	iShiftFluc = nint(CommonPar(1))
	iShiftLag = nint(CommonPar(2))
	
	BpCol = iShiftLag !!! located in first node
	
	Do A=1,NodG 
		ApRow  = (A-1)*iDofT + iShiftFluc
		
		Do p = 1 , NdimE 
			App = ApRow + p
			Bpp = BpCol + p !! it's a kind of identity 
												 	
			Coupling(Bpp,App) = 1
			Coupling(App,Bpp) = 1
			Coupling(Bpp,Bpp) = 1

		end do ! loop Ap
	Enddo !LoopRow

end subroutine





