    !> Element implementation for hyperlastic fibres
    !!
    !! \f[
    !!   \sum_f \mathbf{s}_f \cdot \hat{\mathbf{q}}_f 
    !! \f]
    !!
	!! @param iShiftU = nint(CommonPar(1))
	!! @param iShiftDeltaU = nint(CommonPar(2)) !!! now it's not being used
	!! @param iFemType  = nint(CommonPar(3)) 
	!! @param iFType  = nint(commonPar(4))
	!! @param iMaterial  = nint(commonPar(5))
	!! @param iDamageParType  = nint(commonPar(6))
	!!
    !! @author Rocha, Felipe Figueredo
   
Subroutine nonlinearFibres_qfBased(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, &
									Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	
	
	use funcAux , only : setAEandBEfib, getSliceAllocate
	use globalVariables, only : getF, NdimE, getMaterial, getDamagePar, maxMatPar, maxDamagePar, addToPKhomGen
	use fibresLib
	use ptsGaussLib2, only : setNodG
	use damageFibLib
		
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT, NodElt! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)
!~ 	!   =====   END ARGUMENTS  =======
    
	! ======== COMMONPAR ARGUMENTS =====================
	integer iShiftU !> @param  = nint(CommonPar(1))
	integer iShiftDeltaU !> @var  = nint(CommonPar(2)) !!! now it's not being used
	integer iFemType !> @var  = nint(CommonPar(2)) = nint(commonPar(3))
	integer iFType !> @var  = nint(commonPar(4))
	integer iMaterial !> @var  = nint(commonPar(5))
	integer iDamageParType !> @var = nint(commonPar(6))
	
	!====== OTHER VARIABLES ==============================
	Real*8 :: MatPar(maxMatPar),  damagePar(maxDamagePar)
	Real*8 :: DeltaU(NdimE), afib(NdimE), qfib(NdimE), Lfib, Res(NdimE) , K(NdimE,NdimE)
	Real*8 :: DSfib(NdimE,NdimE), Sfib(NdimE), signal, Area
	real*8 , allocatable ::  SolU(:) , Xel(:) !! all have dimension NodG*NdimE
	integer :: NodG
	integer :: constLaw 

	iShiftU = nint(CommonPar(1))
	iShiftDeltaU = nint(CommonPar(2)) !!! now it's not being used
	iFemType = nint(commonPar(3))
	iFType = nint(commonPar(4))
	iMaterial = nint(commonPar(5))
	iDamageParType = nint(commonPar(6))
	
    call setNodG(iFemtype, NodG)	
	
	if(iDamageParType>-1) then
		call getDamagePar(damagePar,iDamageParType)
	end if
	call getMaterial(constLaw, matPar, iMaterial)
	Area = matpar(3)
 	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	call getSliceAllocate(SolU,Sol1,1,NodG,iShiftU + 1 ,iShiftU + NdimE, iDofT)

	call setFibreKinematic(qfib,DeltaU,afib,Lfib,SolU,Xel,iFtype)
	
	call calcSfib(Sfib,qfib,matPar,constLaw,NDimE)
	call calcDSfib(DSfib,qfib,matPar,constLaw,NDimE)

!~ 	call damageUpdateFib(damagePar,Param,norm2(qfib))
	if(iDamageParType>-1) then
		call damageModifySfib(Sfib, Param)
		call damageModifyDSfib(DSfib, Param)
	end if

	Res = -Area*Sfib	
	K = (Area/Lfib)*DSfib
	
	Res = Res + matmul(K,DeltaU)
	
	call setAEandBEfib(AE,BE,K,Res,NodG,idofT,iShiftU,NdimE)
	
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine nonlinearFibres_qfBasedS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux, only : setCoupling_pureDisplacement, numprint
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

	call setCoupling_pureDisplacement(Coupling,NodG,idofT,iShiftU,NdimE)

end Subroutine



