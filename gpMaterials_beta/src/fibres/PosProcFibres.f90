!Last Modification : 11/03/2014
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Element for Hyperelasticity 
! Still in development
!     ------------------------------------------------------------------

!~ 
Subroutine posProcFibres(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
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
    
    integer :: i
	Integer , parameter :: iShiftField = 6,  NFieldPar = 2
	Real*8 :: MatPar(maxMatPar),  damagePar(maxDamagePar)
	Real*8 :: DeltaU(NdimE), afib(NdimE), qfib(NdimE), Lfib
	Real*8 :: Sfib(NdimE), Area
	real*8 , allocatable ::  SolU(:) , Xel(:) !! all have dimension NodG*NdimE
	integer :: NodG, iMaterial, iDamageParType
	integer :: iShiftUf, constLaw , iFemType, iFtype
	real*8 :: VField(4)
	integer ::  NFields, CodeField(5) , PosField(5)
	
	matPar = 0.0d0

	iShiftUf = nint(CommonPar(1))
	iFemType = nint(commonPar(2))
	iFType = nint(commonPar(3))
	iMaterial = nint(commonPar(4))
	iDamageParType = nint(commonPar(5))
	NFields = nint(CommonPar(6))
	
	do i = 1, NFields
		CodeField(i) = nint(CommonPar( iShiftField + (i-1)*NFieldPar + 1) )
		PosField(i)  = nint(CommonPar( iShiftField + (i-1)*NFieldPar + 2) )
	end do	
	
	
	call setNodG(iFemtype, NodG)	
	
	call getDamagePar(damagePar,iDamageParType)
	call getMaterial(constLaw, matPar, iMaterial)
	Area = matpar(3)
 	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	call getSliceAllocate(SolU,Sol1,1,NodG,iShiftUf + 1 ,iShiftUf + NdimE, iDofT)

	call setFibreKinematic(qfib,DeltaU,afib,Lfib,SolU,Xel,iFtype)
	
	call calcSfib(Sfib,qfib,matPar,constLaw,NDimE)
	
	do i = 1 , NFields
		select case(CodeField(i))
			case(1) !! plot fibres strecht
				VField(i) = norm2(qfib)
			case(2) !! plot fibre tension
				VField(i) = norm2(Sfib)
			case default
				write(0,*) "Code Field not found"
				VField(i) = 0
		end select

		Param(LengthParam + PosField(i) ) = VField(i) 
		
	end do
	
!~ 	AE(1,1) = 1.0d0
!~ 	BE(1) = Sol1(1)	
	
End Subroutine



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine posProcFibresS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd, ipRow , i
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)

!~ 	Coupling=0 
	
!~ 	Coupling(1,1) = 1
	
	return
end Subroutine
