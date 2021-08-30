Subroutine fibresTraction(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
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
	Real*8 :: MatPar(maxMatPar), damagePar(maxDamagePar)
	Real*8 :: F(NdimE,NdimE), afib(NdimE), qfib(NdimE), Lfib
	Real*8 ::  tfib(NdimE), signal
	real*8 , allocatable ::  Xel(:) !! all have dimension NodG*NdimE
	integer :: NodG, iFtype, iMaterial, iDamageParType
	integer :: iShiftTraction, constLaw 

	iShiftTraction = nint(commonPar(1))
	iFtype = nint(commonPar(2)) 
	iMaterial = nint(commonPar(3))
	iDamageParType = nint(commonPar(4))
	
	NodG = NodElt - 1	
	
	call getDamagePar(damagePar,iDamageParType)
	call getMaterial(constLaw, matPar, iMaterial)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	
	call getF(F,iFtype)
	
	afib(1) = Xel(3) - Xel(1)
	afib(2) = Xel(4) - Xel(2)  
	
	Lfib = norm2(afib)
	afib = afib/Lfib
	
	qfib = matmul(F,afib)
	
	call calcSfib(tfib,qfib,matPar,constLaw,NDimE)
	
	call damageUpdateFib(damagePar,Param,norm2(qfib))
	call damageModifySFib(tfib, Param)
		
	Do A=1, NodG !LoopRow
		ApRow  = (A-1)*iDofT + iShiftTraction
				
		if(A==1) then !! 1 or 2
			signal = -1.0d0
		else
			signal = 1.0d0
		end if
		
		Do p = 1 , NdimE
			App = ApRow+p
			
			AE(App,App) = 1.0d0 
			BE(App) = signal*tfib(p) 
			
		end do ! loop Ap
	Enddo !LoopRow
		
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine fibresTractionS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
    
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
    
    Integer :: A,B,q,p, ApRow, BpCol, iShiftTraction
    integer , parameter :: NdimE = 2, NodG = 2
    
    iShiftTraction = nint(commonPar(1))

    do A = 1,NodG
        ApRow = (A-1) * iDofT + iShiftTraction
		
		do p=1,nDimE
			Coupling (ApRow + p, ApRow + p) = 1
		enddo
	enddo

end Subroutine
