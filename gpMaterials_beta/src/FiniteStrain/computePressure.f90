!Last Modification : 11/03/2014
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Element for Hyperelasticity 
! Still in development
!     ------------------------------------------------------------------

!~ 
Subroutine computePressure(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
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
	Real*8 :: GradU(NdimE,NdimE) , F(NdimE,NdimE), F0(NdimE,NdimE), detF, detF0, Fbar(NdimE,NdimE)
	Real*8 :: Ppk(NdimE,NdimE), D(NdimE,NdimE,NdimE,NdimE), sigma(NdimE,NdimE)
	Real*8 :: MatPar(10) , VField(4) , energy, energyMean, dVng, dVtot, pressure
	real*8 , allocatable :: SolU(:) ,  Xel(:), Psi(:,:) , W(:)
	integer :: NGP,pOrder, iShiftP, iShiftU, constLaw
    integer , parameter :: iShiftDamageParam = 0 
	integer :: ishiftDamageCommonPar , IdamageModify, iFemType,  iSimplex,  iBubble,  Nelem 
	type(PtGaussClass) :: PtG, PtG0
	
	matPar = 0.0d0

	iShiftP = nint(commonPar(1))
	iShiftU = nint(commonPar(2))
	iFemType = nint(commonPar(3))
	constLaw = nint(CommonPar(4)) 	
	matPar(1:3) = CommonPar(5:7)
	
	IdamageModify = nint(commonPar(8))  
	
	!! commonPar for the damage starts after the commonPar of the actual function
	ishiftDamageCommonPar = 8
	
	call setFEMtype(iFEMtype,NodG,pOrder,NGP,iSimplex,iBubble)
	
	call getSliceAllocate(SolU,Sol1,1,NodG,iShiftU + 1 ,iShiftU + NdimE, iDofT)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	
	allocate(Psi(NdimE,NGP))
	allocate(W(NGP))
	
	Psi(1,1) = -1.0d0; Psi(2,1) = -1.0d0
	Psi(1,2) = 1.0d0; Psi(2,2) = -1.0d0
	Psi(1,3) = -1.0d0; Psi(2,3) = 1.0d0
	Psi(1,4) = 1.0d0; Psi(2,4) = 1.0d0
	
	!! W is dummy
	W(1) = 0.25d0
	W(2) = 0.25d0
	W(3) = 0.25d0
	W(4) = 0.25d0
		
	call PtG%initWithPsi(Psi,W,Xel,NodG,NdimE,NGP,pOrder,iBubble, iSimplex)
	
	nG = 1
	call PtG0%init(Xel,NodG,NdimE,1,pOrder,iBubble,iSimplex)
	call PtG0%calcGradU(GradU,SolU,nG)
	F0 = deltaKron(1:NdimE,1:NdimE) + GradU
	call calcI3(detF0,F0)
!~ 	
	do nG = 1 , NodG
		call PtG%calcGradU(GradU,SolU,nG)

		F = deltaKron(1:NdimE,1:NdimE) + GradU
!~ 		
		call calcI3(detF,F) 
!~ 		Fbar = ((detF0/detF)**(1.0d0/2.0d0))*F !! for 2D
		Fbar = F !! for 2D

		call calcPpk(Ppk,Fbar,NdimE,matPar,constLaw)		
		call calcD(D,Fbar,NdimE,matPar,constLaw)	
		call damageModifyTangentNew(D,Ppk,commonPar,Param,NdimE, &
			iShiftDamageCommonPar,iShiftDamageParam, nG,IdamageModify)
 		
!~ 		sigma = (1.0d0/detF0)*matmul(Ppk,transpose(Fbar))
		sigma = (1.0d0/detF)*matmul(Ppk,transpose(Fbar))
	
	
		pressure = detF
		
!~ 		pressure = - (sigma(1,1) + sigma(2,2))/3.0d0
 
		ipDof = (nG-1)*idofT + iShiftP  + 1
!~ 		
		AE(ipDof,ipDof) =  1.0d0
		BE(ipDof) =  pressure
!~ 				 
	end do
	
!~ 	AE(1,1) = 1.0d0
!~ 	BE(1) = Sol1(1)
	
	deallocate(Psi)
	
End Subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine computePressureS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	use ptsGaussLib2
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd, ipRow , i
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)

	integer :: femType, NodG, iShiftP, ApRow, A


	Coupling = 0
    
    iShiftP = nint(commonPar(1))
    femType = nint(commonPar(3)) 
    
    call setNodG(femtype, NodG)

    do A = 1,NodG
        ApRow = (A-1) * iDofT + iShiftP + 1
		Coupling (ApRow,ApRow) = 1
	enddo
	
!~ 	
!~ 	Coupling=0 
!~ 	
!~ 	Coupling(1,1) = 1
	
	call numprint(Coupling)
	
	return
end Subroutine
