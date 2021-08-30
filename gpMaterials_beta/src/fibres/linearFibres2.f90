Subroutine linearFibres2(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
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
	Real*8 :: MatPar(maxMatPar), DSx, Sx
	Real*8 ::  Lfib, Ex, Area , afib(NdimE) , K(NdimE,NdimE), Res(NdimE), T(NdimE,NdimE) , Tt(NdimE,NdimE) , du(NdimE)
	real*8 , allocatable :: Xel(:) !! all have dimension NodG*NdimE
	integer :: NodG, iMaterial
	integer :: iShiftU, iShiftDeltaU, constLaw , iFemType

	iShiftU = nint(CommonPar(1))
	iShiftDeltaU = nint(CommonPar(2)) !!! now it's not being used
	iFemType = nint(commonPar(3))
	iMaterial = nint(commonPar(4))
	 
    call setNodG(iFemtype, NodG)	
		
	call getMaterial(constLaw, matPar, iMaterial)
	Area = matpar(3)
 	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)

	!! definition of local2global and global2local transformations
	afib(1) = Xel(3) - Xel(1)
	afib(2) = Xel(4) - Xel(2)  
	
	Lfib = norm2(afib)
	
	afib = afib/Lfib

	! T is global2local
	T(1,1) = afib(1) !! cos
	T(1,2) = afib(2) !! sin
	T(2,1) = -afib(2) !! -sin
	T(2,2) = afib(1) !! cos
	
	Tt = transpose(T) !! local2global
	
	!! compute K
	K(1,1) = 1.0d0
	K(1,2) = 0.0d0
	K(2,1) = 0.0d0
	K(2,2) = 0.0d0
	K = (Area*matpar(1)/Lfib)*K
	
	!! to global : matrices and vector blocks
	K = matmul(matmul(Tt,K),T)
	
!~ 	call numprint(K)

	!! proper global resulting matrix combining the blocks
	Do A=1, NodG !LoopRow
		ApRow  = (A-1)*iDofT + iShiftU
				
		Do p = 1 , NdimE
			App = ApRow+p
			
!~ 			BE(App) = (-1.0d0)**(A+1)*Res(p)
					
			do B= 1 , NodG ! LoopCol ! not considering the simmetry
				BpCol  = (B-1)*iDofT + iShiftU
				
				do q=1,NdimE
					Bqp = BpCol+q
					
					AE(App,Bqp) = (-1.0d0)**(A+B) * K(p,q)
			
				end do ! loop Bq
			Enddo !LoopCol
		end do ! loop Ap
	Enddo !LoopRow
	
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine linearFibres2S(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
	use ptsGaussLib2, only : setNodG
    
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
    
    Integer :: A,B,q,p, ApRow, BpCol
    integer , parameter :: NdimE = 2
    integer :: iShiftU, iShiftDeltaU, iFemType, NodG

	iShiftU = nint(CommonPar(1))
	iShiftDeltaU = nint(CommonPar(2)) !! not being used at the moment
	iFemType = nint(commonPar(3))
	     
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
