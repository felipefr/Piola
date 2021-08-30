!Last Modification : 11/03/2014
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Element for Hyperelasticity 
! Still in development
!     ------------------------------------------------------------------

subroutine damageGenericFibres(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	
	use funcAux
	use globalVariables, only : NdimE, getDamagePar, getMaterial, maxDamagePar, maxMatPar, getF
	use fibresLib
	use ptsGaussLib2, only : setNodG
	use damageFibLib

	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

	!   =====   END ARGUMENTS  =======
         
	Real*8 :: MatPar(maxMatPar),  damagePar(maxDamagePar)
	real*8 , allocatable ::  SolU(:) , Xel(:) !! all have dimension NodG*NdimE
	integer :: NodG, iShiftU, constLaw, iFemType, iMaterial, iDamageParType , iFtype
	real*8 :: lambda

	iShiftU = nint(commonPar(1))
	iFemType = nint(commonPar(2))
	iFtype = nint(commonPar(3))
	iMaterial = nint(commonPar(4))	
	iDamageParType = nint(commonPar(5))
	
    call setNodG(iFemtype, NodG)	
	call getDamagePar(damagePar,iDamageParType)
	call getMaterial(constLaw,matPar,iMaterial)	

	call getSliceAllocate(SolU,Sol1,1,NodG,iShiftU + 1 ,iShiftU + NdimE, iDofT)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	
	lambda = getLambdaFibre(SolU,Xel,iFtype,NdimE)
	call damageUpdateFib(damagePar,Param , lambda)		
		
	AE(1,1) = 1.0d0
	BE(1) = Sol1(1)
	
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine damageGenericFibresS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd, ipRow , i
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)

	Coupling=0 
	
	Coupling(1,1) = 1
	
	return
end subroutine
