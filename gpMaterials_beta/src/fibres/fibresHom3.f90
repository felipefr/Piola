Subroutine fibresHom3(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
	use globalVariables, only : getF, NdimE, getMaterial, getDamagePar, maxMatPar, maxDamagePar, addToPKhomGen
	use fibresLib
	use ptsGaussLib2, only : setNodG
	use damageFibLib
	use DETERMINANT , only : iElem
	
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT, NodElt! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

!~ 	!   =====   END ARGUMENTS  =======
	Real*8 :: MatPar(maxMatPar),  damagePar(maxDamagePar)
	Real*8 :: DeltaU(NdimE), afib(NdimE), qfib(NdimE), Lfib,  PK(NdimE,NdimE)
	Real*8 :: Sfib(NdimE), AreaFib, Vfib
	real*8 , allocatable ::  SolU(:) , Xel(:) !! all have dimension NodG*NdimE
	integer :: NodG, iMaterial, iDamageParType
	integer :: iShiftU, iShiftDeltaU, constLaw , iFemType, iFtype, iElemBegin, Nelem

	iShiftU = nint(CommonPar(1))
	iShiftDeltaU = nint(CommonPar(2)) !!! now it's not being used
	iFemType = nint(commonPar(3))
	iFType = nint(commonPar(4))
	iMaterial = nint(commonPar(5))
	iDamageParType = nint(commonPar(6))
	iElemBegin = nint(commonPar(7))
	Nelem = nint(commonPar(8))
	
    call setNodG(iFemtype, NodG)	
	
	call getDamagePar(damagePar,iDamageParType)
	call getMaterial(constLaw, matPar, iMaterial)
	AreaFib = matpar(3)
 	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	call getSliceAllocate(SolU,Sol1,1,NodG,iShiftU + 1 ,iShiftU + NdimE, iDofT)

	call setFibreKinematic(qfib,DeltaU,afib,Lfib,SolU,Xel,iFtype)

	call calcSfib(Sfib,qfib,matPar,constLaw,NDimE)

!~ 	call numprint(Sfib)
!~ 	pause

	call damageModifySfib(Sfib, Param)

	Vfib = AreaFib*Lfib
	
	PK(1,1) = Vfib*Sfib(1)*afib(1)
	PK(1,2) = Vfib*Sfib(1)*afib(2)
	PK(2,1) = Vfib*Sfib(2)*afib(1)
	PK(2,2) = Vfib*Sfib(2)*afib(2)	

	call addToPKhomGen(PK,Vfib,iElem,iElemBegin,Nelem)
	
	AE(1,1) = 1.0d0
	BE(1) = Sol1(1)

end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine fibresHom3S(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
    
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)

	Coupling (1, 1) = 1

end Subroutine
