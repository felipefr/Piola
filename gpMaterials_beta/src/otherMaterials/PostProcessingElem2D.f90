!Last Modification : 11/03/2014
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Element for Hyperelasticity 
! Still in development
!     ------------------------------------------------------------------

!~ 
Subroutine postProcessingElem2D(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
	use finiteStrainLib !! New way
	use damageLib
	use multiscaleLib
	use ptsGaussLib2

	implicit none

	!	PARAMETERS
	
	integer :: NParamField = 2

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT,NodElt ! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

	!   =====   END ARGUMENTS  =======
         
	Integer :: nG,i,j,k,l,m,n,p,q, A, B, ApRow, BpCol , App,Bqp, ipDim, ipDimE, ipDof,ipDofmod, NodG ! counters 
	Integer , parameter :: NDimE = 2
	Real*8 :: GradU(NdimE,NdimE) , F(NdimE,NdimE)
	Real*8 :: Ppk(NdimE,NdimE), D(NdimE,NdimE,NdimE,NdimE), PpkMean(NdimE,NdimE), DMean(NdimE,NdimE,NdimE,NdimE)
	Real*8 :: MatPar(10) , VField(4) , energy, energyMean, dVng, dVtot
	real*8 , allocatable :: SolU(:) ,  Xel(:) 
	integer :: NparamElem, NGP,pOrder , IisPeriodic , iShiftU, iShiftParam, constLaw , NFields, CodeField(5) , PosField(5) !! maximum 5 at now
	logical :: isPeriodic, ifAnalytical
	integer, save :: iElem = 0
    integer , parameter :: OUnitP11hom = 48, OUnitP12hom = 49, OUnitP21hom = 50
    integer , parameter :: OUnitP22hom = 51, OUnitDamage = 52, OUnitVonMises = 53
    integer , parameter :: OUnitDisplacement = 54, OUnitRd = 55, OUnitQd = 56
	integer , parameter :: iShiftDamageParam = 0 
	integer :: ishiftDamageCommonPar , IdamageModify, iFemType,  iSimplex,  iBubble,  Nelem 
	Real*8 :: normal(2), tangent(2) ,  theta , alpha, beta(2)
	type(ptGaussClass) :: PtG

	matPar = 0.0d0


	iShiftU = nint(commonPar(1))
	iFemType = nint(commonPar(2))
	Nelem = nint(commonPar(3))	
	constLaw = nint(CommonPar(4)) 	
	matPar(1:3) = CommonPar(5:7)
	NFields = nint(CommonPar(8))
	do i = 1, NFields
		CodeField(i) = nint(CommonPar( 8 + (i-1)*NParamField + 1) )
		PosField(i)  = nint(CommonPar( 8 + (i-1)*NParamField + 2) )
	end do
	
	IdamageModify = nint(commonPar(8 + NFields*NParamField + 1))  
	
	!! commonPar for the damage starts after the commonPar of the actual function
	ishiftDamageCommonPar = 8 + NFields*NParamField + 1
	
	
	call setFEMtype(iFEMtype,NodG,pOrder,NGP,iSimplex,iBubble)
	
	call getSliceAllocate(SolU,Sol1,1,NodG,iShiftU + 1 ,iShiftU + NdimE, iDofT)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	
	
	F = deltaKron(1:NdimE,1:NdimE) + GradU
	
	
	call PtG%init(Xel,NodG,NdimE,NGP,pOrder,iBubble, iSimplex)
	
	PpkMean = 0.0d0
	DMean = 0.0d0
	energyMean = 0.0d0
	dVtot = 0.0d0
	
	do nG = 1 , NGP
		call PtG%calcGradU(GradU,SolU,nG)
		
		F = deltaKron(1:NdimE,1:NdimE) + GradU
		
		call calcPpk(Ppk,F,NdimE,matPar,constLaw)		
		call calcD(D,F,NdimE,matPar,constLaw)	
		call damageModifyTangentNew(D,Ppk,commonPar,Param,NdimE, &
			iShiftDamageCommonPar,iShiftDamageParam, nG,IdamageModify)
		
		call strainEnergy(energy,F,matPar,constLaw)
		
		dVng=ptG%dV(nG)
		
		dVtot = dVtot + dVng
		PpkMean = PpkMean + dVng*Ppk
		DMean = DMean + dVng*D
		energyMean = energyMean + dVng*energy 	
	end do
	
	PpkMean = PpkMean/dVtot
	DMean = DMean/dVtot
	energyMean = energyMean/dVtot
	
	do i = 1 , NFields
		select case(CodeField(i))
			case(1) !! energy 
				VField(i) = energy
			case(2) !! vol
				VField(i) = dVtot
			case(3) !! Pnn				
				theta = datan(0.3d0/1.0d0)
				normal(1) = -dsin(theta)
				normal(2) = dcos(theta)				
				VField(i) = dot_product(matmul(PpkMean,normal),normal)
			case(4) !! Pnt
				theta = datan(0.3d0/1.0d0)
				normal(1) = -dsin(theta)
				normal(2) = dcos(theta)
				tangent(1) = dcos(theta)
				tangent(2) = dsin(theta)
				VField(i) = dot_product(matmul(PpkMean,normal),tangent)
			
			case(5) !! Pnb
				alpha = 0.5d0*PI
				theta = datan(0.3d0/1.0d0)
				normal(1) = -dsin(theta)
				normal(2) = dcos(theta)
				beta(1) = dcos(theta + alpha)
				beta(2) = dsin(theta + alpha)
				
				VField(i) = dot_product(matmul(PpkMean,normal),beta)
!~ 				VField(i) = dot_product(matmul(PpkMean,normal),normal)
				
			case(6) !! Pnb update (in case of fibers for example)
				alpha = 0.5d0*PI
				theta = datan(0.3d0/1.0d0)
				normal(1) = -dsin(theta)
				normal(2) = dcos(theta)
				beta(1) = dcos(theta + alpha)
				beta(2) = dsin(theta + alpha)
				
				VField(i) = Param(LengthParam + PosField(i) ) + dot_product(matmul(PpkMean,normal),beta)
				
			case(7) !! Bifurcation Indicator
				alpha = 0.5d0*PI
				theta = datan(0.3d0/1.0d0)
				normal(1) = -dsin(theta)
				normal(2) = dcos(theta)
				beta(1) = dcos(theta + alpha)
				beta(2) = dsin(theta + alpha)
				
				VField(i) = 0.0d0
				do j = 1, NdimE
				do k = 1, NdimE
					VField(i) = VField(i) + GradU(j,k)*beta(j)*normal(k)
				end do
				end do
				write(*,*) VField(i)
!~ 				PAUSE

			case(11) !! P11
				VField(i) = PpkMean(1,1)
			case(12) !! P12
				VField(i) = PpkMean(1,2)
			case(13) !! P21
				VField(i) = PpkMean(2,1)
			case(14) !! P22
				VField(i) = PpkMean(2,2)
			
			case default
				write(0,*) "Code Field not found"
				VField(i) = 0
		end select
		
		Param(LengthParam + PosField(i) ) = VField(i)

	end do
	
	AE(1,1) = 1.0d0
	BE(1) = Sol1(1)	
	
End Subroutine



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine postProcessingElem2DS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd, ipRow , i
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)

	Coupling=0 
	
	Coupling(1,1) = 1
	
	return
end Subroutine
