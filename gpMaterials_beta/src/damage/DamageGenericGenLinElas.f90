    !> Generic element for damage
    !!
	!! @param iShiftU = nint(CommonPar(1))
	!! @param iFemType  = nint(CommonPar(2)) 
	!! @param iMaterial  = nint(commonPar(3))
	!! @param iDamageParType  = nint(commonPar(4))
	!!
    !! @author Rocha, Felipe Figueredo
    
Subroutine damageGenericGenLinElas(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, Sol1, &
									CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	
	use funcAux
	use globalVariables, only : NdimE, getDamagePar, getMaterial, maxDamagePar, maxMatPar, getF
	use damageNewLib
	use ptsGaussLib2
	use linElasLib

	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

	!   =====   END ARGUMENTS  =======
         
	Integer :: nG, i,j
	Real*8 :: GradU(NdimE,NdimE) , epsMat(NdimE*NdimE-1)
	Real*8 ::energy, MatPar(maxMatPar),  damagePar(maxDamagePar)
	real*8 , allocatable ::  SolU(:) , Xel(:) !! all have dimension NodG*NdimE
	integer :: NodG, iShiftU, constLaw, iFemType, iMaterial, iDamageParType, iBubble, iSimplex, pOrder, NGP
	type(ptGaussClass) :: PtG

	iShiftU = nint(commonPar(1))
	iFemType = nint(commonPar(2))
	iMaterial = nint(commonPar(3))	
	iDamageParType = nint(commonPar(4))
	
	call getDamagePar(damagePar,iDamageParType)
	call getMaterial(constLaw,matPar,iMaterial)	
	call setFEMtype(iFEMtype,NodG,pOrder,NGP,iSimplex,iBubble)

	call getSliceAllocate(SolU,Sol1,1,NodG,iShiftU + 1 ,iShiftU + NdimE, iDofT)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	
	call PtG%init(Xel,NodG,NdimE,NGP,pOrder,iBubble, iSimplex) 		

	do nG = 1 , NGP
		call PtG%calcGradU(GradU,SolU,nG)	
		
		call strainEnergyLinElas(energy,GradU,matPar)
		
		call damageUpdate(damagePar,Param , energy , nG)		
		
		call writeDamageFromParam(Param,energy,nG)
		
	end do
	
	AE(1,1) = 1.0d0
	BE(1) = Sol1(1)
	
	
End

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine damageGenericGenLinElasS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd, ipRow , i
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)

	Coupling=0 
	
	Coupling(1,1) = 1
	
	return
end Subroutine
