    !> Incorporates minRestriction for fibres network 
    !! 
	!! @param iShiftFluc = nint(CommonPar(1))
	!! @param iShiftLag = nint(CommonPar(2))
	!! @param iMaterial  = nint(CommonPar(3)) 
	!! @param pen  = commonPar(4)
	!! @param eps  = commonPar(5)

    !! @author Rocha, Felipe Figueredo
    
!     ------------------------------------------------------------------
Subroutine minRestrictionFibres(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, &
						Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
	use globalVariables, only : NdimE, getMaterial, maxMatPar
	use fibresLib
		
	IMPLICIT NONE
	
	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	! Reals Vectors and matrices
	Real*8  :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), &
		Param(*), JParam(*), CommonPar(*)
	
	!   =====   END ARGUMENTS  =======
     		
	Integer :: i,j,k,l,m,n,p,q, A, B, ApRow, BpCol , App,Bqp! counters 
	Integer :: NdofLag
	integer, parameter ::  NodG = 2, NodLag = 3
	Real*8 :: MatPar(maxMatPar)
	real*8 , allocatable ::  Xel(:) 
	Real*8 :: eps, pen 
	Real*8 :: afib(NdimE), Areafib, Lfib, signal
	integer :: iShiftFluc, iShiftLag , iMaterial, constLaw
	
	iShiftFluc = nint(CommonPar(1))
	iShiftLag = nint(CommonPar(2))
	iMaterial = nint(CommonPar(3))
	pen = commonPar(4)
	eps = commonPar(5)

	call getMaterial(constLaw, matPar, iMaterial)
	Areafib = matpar(3)
	
 	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	call get_afib(afib,Lfib,Xel)
	
	NdofLag = NdimE*NdimE
				
	B = NodLag
	BpCol = (B-1)*iDofT + iShiftLag
	
	AE = 0.0d0
	
	Do A=1,NodG 
		ApRow  = (A-1)*iDofT + iShiftFluc
		
		Do p = 1 , NdimE
			App = ApRow+p
												 
			Do q = 1 , NdofLag
				Bqp = BpCol + q
				
				call getIJfromK(i,j,q)
				
				signal = (-1.0d00)**A
			
				if(i == p) then
					AE(Bqp,App) = pen*signal*Areafib*afib(j)
					AE(App,Bqp) = AE(Bqp,App) 
				end if
				AE(Bqp,Bqp) = eps
				BE(Bqp) = eps*Sol1(Bqp)
			end do
		end do ! loop Ap
	Enddo !LoopRow

end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine minRestrictionFibresS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	use funcAux
	use globalVariables, only : NdimE
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)
	integer ::  p , q, App, ApRow, BpCol, Bqp, A,B, iShiftFluc, iShiftLag , NdofLag
	integer , parameter :: NodG = 2, NodLag = 3
	
	iShiftFluc = nint(CommonPar(1))
	iShiftLag = nint(CommonPar(2))
	
	NdofLag = NdimE*NdimE
	
	B = NodLag
	BpCol = (B-1)*iDofT + iShiftLag
	
	Do A=1,NodG
		ApRow  = (A-1)*iDofT + iShiftFluc
		
		Do p = 1 , NdimE
			App = ApRow+p
												 
			Do q = 1 , NdofLag
				Bqp = BpCol + q
			
				Coupling(App,Bqp) = 1
				Coupling(Bqp,App) = 1
				Coupling(Bqp,Bqp) = 1
			end do
		end do ! loop Ap
	Enddo !LoopRow

!~ 	
end subroutine





