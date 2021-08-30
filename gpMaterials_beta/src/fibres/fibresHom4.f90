!> Element implementation for homogenisation of P for hyperelastic networks
!!
!! @param iOption = nint(CommonPar(1)) , 1 - divide by volume and write in file , otherwise just acumulates
!! @param iShiftUf = nint(CommonPar(2))
!! @param iFemType  = nint(CommonPar(3)) 
!! @param iFType  = nint(commonPar(4))
!! @param iMaterial  = nint(commonPar(5))
!! @param iDamageParType  = nint(commonPar(6)) ( <0 if not defined )
!!
!! @author Rocha, Felipe Figueredo
    
Subroutine fibresHom4(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
	use globalVariables, only : getF, NdimE, getMaterial, getDamagePar, maxMatPar, &
							maxDamagePar, contributePKhom, writePKhom, allocatePKhom
	use fibresLib
	use ptsGaussLib2, only : setNodG
	use damageFibLib
	
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT, NodElt! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

!~ 	!   =====   END ARGUMENTS  =======
	Real*8 :: MatPar(maxMatPar),  damagePar(maxDamagePar)
	Real*8 :: DeltaU(NdimE), afib(NdimE), qfib(NdimE), Lfib,  PK(NdimE,NdimE)
	Real*8 :: Sfib(NdimE), AreaFib, Vfib
	real*8 , allocatable ::  SolUf(:) , Xel(:) !! all have dimension NodG*NdimE
	integer :: NodG, iMaterial, iDamageParType
	integer :: iShiftUf, constLaw , iFemType, iFtype, iOption

	iOption = nint(CommonPar(1))
	
	if(iOption ==1) then
		AE(1,1) = 1.0d0
		BE(1) = Sol1(1)
		call writePKhom()
		return
	end if
	
	iShiftUf = nint(CommonPar(2))
	iFemType = nint(commonPar(3))
	iFType = nint(commonPar(4))
	iMaterial = nint(commonPar(5))
	iDamageParType = nint(commonPar(6))
	
    call setNodG(iFemtype, NodG)	
	
	if(iDamageParType>-1) then
		call getDamagePar(damagePar,iDamageParType)
	end if
	call getMaterial(constLaw, matPar, iMaterial)
	AreaFib = matpar(3)
 	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	call getSliceAllocate(SolUf,Sol1,1,NodG,iShiftUf + 1 ,iShiftUf + NdimE, iDofT)

	call setFibreKinematic(qfib,DeltaU,afib,Lfib,SolUf,Xel,iFtype)

	call calcSfib(Sfib,qfib,matPar,constLaw,NDimE)

	if(iDamageParType>-1) then
		call damageModifySfib(Sfib, Param)
	end if

	Vfib = AreaFib*Lfib
	
	PK(1,1) = Vfib*Sfib(1)*afib(1)
	PK(1,2) = Vfib*Sfib(1)*afib(2)
	PK(2,1) = Vfib*Sfib(2)*afib(1)
	PK(2,2) = Vfib*Sfib(2)*afib(2)	
	
	call allocatePKhom() !! desconsidered if already allocated 

	call contributePKhom(PK,Vfib)
	
	AE(1,1) = 1.0d0
	BE(1) = Sol1(1)

end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine fibresHom4S(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
    
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)

	Coupling (1, 1) = 1

end Subroutine
