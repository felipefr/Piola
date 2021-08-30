!Last Modification : 11/03/2014
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Element for Hyperelasticity 
! Still in development
!     ------------------------------------------------------------------

Subroutine damageGeneric(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
	use finiteStrainLib
	use damageLib
	use ptsGaussLib2
	
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

	!   =====   END ARGUMENTS  =======
         
	Integer :: nG, i,j, ipDim, ipDimE, ipDof,ipDofmod ! counters 
	Integer , Parameter :: NDimE = 2 
	Real*8 :: GradU(NdimE,NdimE) , F(NdimE,NdimE)
	Real*8 ::energy, MatPar(10)
	real*8 , allocatable ::  SolU(:) , Xel(:) !! all have dimension NodG*NdimE
	integer :: NodG, iShiftDamage, iShiftU, constLaw, iFemType, iBubble, iSimplex, pOrder, NGP
	Real*8 :: RPar(1), HPar(3), H 
	integer , parameter :: iShiftDamageCommonPar = 6 , iShiftDamageParam = 0
	type(ptGaussClass) :: PtG
!~ 
	matPar = 0.0d0
	F = 0.0d0
	
	iShiftU = nint(commonPar(1))
	iFemType = nint(commonPar(2))
	constLaw = nint(CommonPar(3)) 	
	matPar(1:3) = CommonPar(4:6)
	
	call setFEMtype(iFEMtype,NodG,pOrder,NGP,iSimplex,iBubble)

	call getSliceAllocate(SolU,Sol1,1,NodG,iShiftU + 1 ,iShiftU + NdimE, iDofT)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	
	call PtG%init(Xel,NodG,NdimE,NGP,pOrder,iBubble, iSimplex)
!~ 		

	do nG = 1 , NGP
		call PtG%calcGradU(GradU,SolU,nG)
!~ 		
		F = deltaKron(1:NdimE,1:NdimE) + GradU
!~ 		
		call strainEnergy(energy,F,matPar,constLaw)
!~ 		
		call damageUpdate(CommonPar,Param,energy, iShiftDamageCommonPar, iShiftDamageParam, nG)		
	end do
	
	AE(1,1) = 1.0d0
	BE(1) = Sol1(1)
	
End

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine damageGenericS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd, ipRow , i
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)

	Coupling=0 
	
	Coupling(1,1) = 1
	
	return
end Subroutine

