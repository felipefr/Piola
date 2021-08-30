    !> Generic element for linear elasticity with discontinuity through EFEM (with no static condesation for beta)
    !!
	!! @param iShiftU = nint(CommonPar(1))
	!! @param iShiftDeltaU = nint(CommonPar(2))
	!! @param iShiftBeta = nint(CommonPar(3))
	!! @param iShiftDeltaBeta = nint(CommonPar(4))
	!! @param iFemType  = nint(CommonPar(5)) 
	!! @param iFType  = nint(commonPar(6))
	!! @param iMaterial  = nint(commonPar(7))
	!! @param iDamageParType  = nint(commonPar(8))
	!!
    !! @author Rocha, Felipe Figueredo

Subroutine LinElas_EFEM_voigt(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
	use globalVariables, only : getF, NdimE, getDamagePar, getMaterial, maxDamagePar, maxMatPar
	use damageNewLib
	use ptsGaussLib2
	use linElasLib
	
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT, NodElt! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

!~ 	!   =====   END ARGUMENTS  =======
         
	Integer :: nG, i,j,k,l,m,n,p,q, A, B, ApRow, BpCol , App,Bqp, ApDim, BpDim, ipDim, ipDimE, ipDof,ipDofmod ! counters 
	Integer :: constLaw ! if NparamDamage = 0 , no damage
	Real*8 :: dV , Det , MatPar(maxMatPar), damagePar(maxDamagePar)
	Real*8 :: GradU(NdimE,NdimE)
	Real*8 :: sigmaMat(NdimE*NdimE-1), epsMat(NdimE*NdimE-1)
	real*8 , allocatable ::  SolU(:), SolU0(:), Xel(:)  !! all have dimension NodG*NdimE
	integer :: pOrder , NGP, NodG, iFEMtype , iFtype, iMaterial, iDamageParType, iSimplex, iBubble
	integer :: iShiftU, iShiftBeta, iShiftDeltaU, iShiftDeltaBeta, iShiftInfoCrack, SolitaryNode
	type(ptGaussClass) :: PtG
	real*8 :: beta(NdimE), Dts(NdimE,NdimE), betaOtimesBeta(NdimE,NdimE)
	real*8 :: Auu(3*NdimE,3*NdimE), Aub(3*NdimE,NdimE), Abu(NdimE,3*NdimE), Abb(NdimE,NdimE), Bu(3*NdimE), Bb(NdimE), normal(NdimE)
	real*8 :: cTs, aTs, damageTs, AreaS, alphaAux1 , alphaAux2
	real*8 :: Bmat(NdimE*NdimE-1,3*NdimE), Dmat(NdimE*NdimE-1,NdimE*NdimE-1), Dmat0(NdimE*NdimE-1,NdimE*NdimE-1), & 
			  BmatT(3*NdimE,NdimE*NdimE-1), BmatBeta(NdimE*NdimE-1,NdimE), BmatBetaT(NdimE,NdimE*NdimE-1) 
	integer , parameter :: NodBeta = 1
	integer :: flagBif, isEFEM
	logical :: isBif
	real*8 :: normBeta,  nu, E
	
	
	iShiftU = nint(CommonPar(1))
	iShiftDeltaU = nint(CommonPar(2))
	iShiftBeta = nint(CommonPar(3))
	iShiftDeltaBeta = nint(CommonPar(4))
	iFEMtype =  nint(commonPar(5)) 
	iMaterial = nint(commonPar(6))
	iDamageParType = nint(commonPar(7))
	isEFEM = nint(commonPar(8))
	aTs = commonPar(9)
	cTs = commonPar(10)
	iShiftInfoCrack = nint(commonPar(11))
	
	if(iDamageParType>0) call getDamagePar(damagePar,iDamageParType)
	call getMaterial(constLaw, matPar, iMaterial)
	call setFEMtype(iFEMtype,NodG,pOrder,NGP,iSimplex,iBubble)

	call getSliceAllocate(SolU,Sol1,1,NodG,iShiftU + 1 ,iShiftU + NdimE, iDofT)
	call getSliceAllocate(SolU0,Sol0,1,NodG,iShiftU + 1 ,iShiftU + NdimE, iDofT)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	
	call PtG%init(Xel,NodG,NdimE,NGP,pOrder,iBubble, iSimplex)

	AreaS = 0.5d0
	solitaryNode = nint(Param(iShiftInfoCrack + 3))
	normal(1) = 1.0d0
	normal(2) = 0.0d0

	flagBif = nint(Param(LengthParam + iShiftInfoCrack + 4))
	E = matpar(1)
	nu = matpar(2)
	
!~ 	!call buildElasticityTensor(D,E,nu,NdimE)
 	call buildElasticityTensorMatSym(Dmat0,E,nu,NdimE)
 	
	if(flagBif < 1 .and. isEFEM > 0  ) then
			
		call PtG%calcGradU(GradU,SolU0,1)

		call getEpsMatSym(epsMat,GradU)
		
		sigmaMat = matmul(Dmat,epsMat)
		
		if(sigmaMat(1) > cTs) then
			flagBif = 2  !!! it means that the bifurcation was detected at this step
			Param(LengthParam + iShiftInfoCrack + 4) = 1
		end if
	
	end if
 
	if(flagBif > 0 .and. isEFEM > 0) then
			
		beta(1:NdimE) = Sol1(iShiftBeta + 1 : iShiftBeta + NdimE) !!! in the first node, temporary
		normBeta = norm2(beta)
	
		damageTs = aTs*normBeta
		
		call VotimesW(betaOtimesBeta,beta,beta)
		
		alphaAux1 = -cTs*aTs/normBeta
		alphaAux2 = cTs*(1.0d0 - aTs*normBeta)
		
		Dts = alphaAux1*betaOtimesBeta + alphaAux2*deltaKron(1:NdimE,1:NdimE)
		
!~ 	!	Dts = -aTs*cTs*deltaKron(1:NdimE,1:NdimE)

	end if
		
	Param(LengthParam + iShiftInfoCrack + 1 : LengthParam + iShiftInfoCrack + NdimE) =  beta(1:NdimE)

	AE = 0.0d0
	BE = 0.0d0
	Bu = 0.0d0
	Bb = 0.0d0
	Auu = 0.0d0
	Aub = 0.0d0
	Abu = 0.0d0
	Abb = 0.0d0

	
	do nG = 1, NGP ! LoopGauss
		
		call PtG%calcGradU(GradU,SolU,nG)
						
		call PtG%getBmatSym(Bmat,nG)
		BmatT = transpose(Bmat)				
		
		BmatBeta = Bmat(:,5:6) !!! just the last node
		BmatBetaT = transpose(BmatBeta) !!! just the last node
		
		call damageModifyTangentLinElas(Dmat,Dmat0,damagePar,Param,nG)
!~ 		Dmat = Dmat0
		
		dV=ptG%dV(nG)
								
		Bu = 0.0d0
		Auu = matmul(BmatT,matmul(Dmat,Bmat))*dV
!~ 
		if(flagBif > 0 .and. isEFEM > 0) then
			Aub = -matmul(BmatT,matmul(Dmat,BmatBeta))*dV
			Abu = transpose(Aub)
			
			Bb = 0.0d0 
			Abb = matmul(BmatBetaT,matmul(Dmat,BmatBeta))*dV + Dts*AreaS
		else
			Abb = deltaKron(1:NdimE,1:NdimE)
		end if
 		
	EndDo !LoopGauss

	call setBlockToVector(BE,Bu,NodG,iDofT, NdimE, iShiftDeltaU ) !(A,B, NodG, iDofA, iDofB, iShiftA) , B in A
	call setBlockToMatrixSquare(AE,Auu, NodG,iDofT, NdimE, iShiftDeltaU ) !(A,B, NodG, iDofA, iDofB, iShiftA) , B in A
	
	if(isEFEM > 0 ) then
		call setBlockToMatrixSquare(AE,Abb, 1,iDofT, NdimE, iShiftDeltaBeta ) !(A,B, NodG, iDofA, iDofB, iShiftA) , B in A
		call setBlockToVector(BE,Bb,1,iDofT, NdimE, iShiftDeltaBeta ) !(A,B, NodG, iDofA, iDofB, iShiftA) , B in A
		call setBlockToMatrix(AE,Aub, NodG,1, iDofT,iDofT, NdimE,NdimE, iShiftDeltaU, iShiftDeltaBeta) ! (A,B, NodGi,NodGj, iDofAi,iDofAj, iDofBi,iDofBj, iShiftAi, iShiftAj)
		call setBlockToMatrix(AE,Abu, 1, NodG, iDofT, iDofT, NdimE, NdimE, iShiftDeltaBeta, iShiftDeltaU ) ! (A,B, NodGi,NodGj, iDofAi,iDofAj, iDofBi,iDofBj, iShiftAi, iShiftAj)
	end if
	
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine LinElas_EFEM_voigtS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
    use ptsGaussLib2
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
    
    Integer :: A,B,q,p, ApRow, BpCol,iShiftU , iFemtype, NodG, iShiftBeta, isEFEM, iShiftDeltaBeta, iShiftDeltaU
    integer , parameter :: NdimE = 2

    Coupling = 0
    
	iShiftU = nint(CommonPar(1))
	iShiftDeltaU = nint(CommonPar(2))
	iShiftBeta = nint(CommonPar(3))
	iShiftDeltaBeta = nint(CommonPar(4))
	iFEMtype =  nint(commonPar(5)) 
	isEFEM = nint(commonPar(8))
    
    call setNodG(iFemtype, NodG)
    
    call setSimbolicBlockToMatrixSquare(Coupling,NodG,iDofT, NdimE, iShiftDeltaU )
	
	if(isEFEM > 0 ) then
		call setSimbolicBlockToMatrixSquare(Coupling,1,iDofT, NdimE, iShiftDeltaBeta ) 
		call setSimbolicBlockToMatrix(Coupling,NodG,1, iDofT,iDofT, NdimE,NdimE, iShiftDeltaU, iShiftDeltaBeta) 
		call setSimbolicBlockToMatrix(Coupling,1, NodG, iDofT, iDofT, NdimE, NdimE, iShiftDeltaBeta, iShiftDeltaU ) 
	end if
	
end Subroutine


