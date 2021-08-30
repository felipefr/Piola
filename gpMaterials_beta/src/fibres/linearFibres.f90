Subroutine linearFibres(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
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
	Real*8 ::  Lfib, E, Area , afib(NdimE) , r
	real*8 , allocatable :: Xel(:) !! all have dimension NodG*NdimE
	integer :: NodG, iMaterial
	integer :: iShiftU, constLaw , iFemType

	iShiftU = nint(CommonPar(1))
	iFemType = nint(commonPar(2))
	iMaterial = nint(commonPar(3))
	     
    call setNodG(iFemtype, NodG)	
		
	call getMaterial(constLaw, matPar, iMaterial)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)

	afib(1) = Xel(3) - Xel(1)
	afib(2) = Xel(4) - Xel(2)  
	
	Lfib = norm2(afib)
	
	afib = afib/Lfib
	
	E = matPar(1)
	Area = matPar(2)
	r = E*Area/Lfib
	
		
	Do A=1, NodG !LoopRow
		ApRow  = (A-1)*iDofT + iShiftU
				
		Do p = 1 , NdimE
			App = ApRow+p
			
!~ 			BE(App) = (-1.0d0)**(A+1) * r*afib(p)
					
			do B= 1 , NodG ! LoopCol ! not considering the simmetry
				BpCol  = (B-1)*iDofT + iShiftU
				
				do q=1,NdimE
					Bqp = BpCol+q
					
					AE(App,Bqp) = (-1.0d0)**(A+B) * r *afib(p) *afib(q) 
			
				end do ! loop Bq
			Enddo !LoopCol
		end do ! loop Ap
	Enddo !LoopRow
	
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine linearFibresS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
	use ptsGaussLib2, only : setNodG
    
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
    
    Integer :: A,B,q,p, ApRow, BpCol
    integer , parameter :: NdimE = 2
    integer :: iShiftU, iFemType, NodG

	iShiftU = nint(CommonPar(1))
	iFemType = nint(commonPar(2))
	     
    call setNodG(iFemtype, NodG)

    do A = 1,NodG
        ApRow = (A-1) * iDofT + iShiftU
		
		do B = 1,NodG
			BpCol = (B-1) * iDofT + iShiftU

			do p=1,nDimE
			do q=1,nDimE

				Coupling (ApRow + p, BpCol + q) = 1

			enddo
            enddo
		enddo
	enddo


end Subroutine
