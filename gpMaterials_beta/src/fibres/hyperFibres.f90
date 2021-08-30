Subroutine HyperFibres(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
	use globalVariables, only : getF, NdimE, getMaterial, maxMatPar
	use fibresLib
	
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT, NodElt! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

!~ 	!   =====   END ARGUMENTS  =======
         
	Integer :: A, B, p, q, ApRow, BpCol, App, Bqp
	Real*8 :: MatPar(maxMatPar)
	Real*8 :: F(NdimE,NdimE), DeltaUf(NdimE), afib(NdimE), qfib(NdimE), Lfib
	Real*8 :: DSfib(NdimE,NdimE), Sfib(NdimE), signal
	real*8 , allocatable ::  SolUf(:) , Xel(:) !! all have dimension NodG*NdimE
	integer :: NodG, iFtype, iMaterial
	integer :: iShiftUf, iShiftDeltaU, constLaw 

	iShiftUf = nint(CommonPar(1))
	iShiftDeltaU = nint(CommonPar(2))
	iFtype = nint(commonPar(3)) 
	iMaterial = nint(commonPar(4))
	
	NodG = NodElt	
		
	call getMaterial(constLaw, matPar, iMaterial)

	call getSliceAllocate(SolUf,Sol1,1,NodG,iShiftUf + 1 ,iShiftUf + NdimE, iDofT)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	
	call getF(F,iFtype)
	
	DeltaUf(1) = SolUf(3) - SolUf(1)
	DeltaUf(2) = SolUf(4) - SolUf(2)
	
	afib(1) = Xel(3) - Xel(1)
	afib(2) = Xel(4) - Xel(2)  
	
	Lfib = norm2(afib)
	afib = afib/Lfib
	
	qfib = matmul(F,afib) + DeltaUf/Lfib
	
	call calcSfib(Sfib,qfib,matPar,constLaw,NDimE)
	call calcDSfib(DSfib,qfib,matPar,constLaw,NDimE)
		
	Do A=1, NodG !LoopRow
		ApRow  = (A-1)*iDofT + iShiftDeltaU
				
		Do p = 1 , NdimE
			App = ApRow+p
			
			signal = (-1.0d0)**A
			
			BE(App) = - signal*Sfib(p) 
			
			do B= 1 , NodG ! LoopCol ! not considering the simmetry
				BpCol  = (B-1)*iDofT + iShiftDeltaU
				
				signal = (-1.0d0)**(A + B)
				do q=1,NdimE
					Bqp = BpCol+q
					
					AE(App,Bqp) = signal*DSfib(p,q)/Lfib  
			
				end do ! loop Bq
			Enddo !LoopCol
		end do ! loop Ap
	Enddo !LoopRow
	
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine hyperFibresS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
    
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
    
    Integer :: A,B,q,p, ApRow, BpCol,iShiftU , iShiftDeltaU
    integer , parameter :: NdimE = 2, NodG = 2

    iShiftDeltaU = nint(commonPar(2))

    do A = 1,NodG
        ApRow = (A-1) * iDofT + iShiftDeltaU
		
		do B = 1,NodG
			BpCol = (B-1) * iDofT + iShiftDeltaU

			do p=1,nDimE
			do q=1,nDimE

				Coupling (ApRow + p, BpCol + q) = 1

			enddo
            enddo
		enddo
	enddo


end Subroutine
