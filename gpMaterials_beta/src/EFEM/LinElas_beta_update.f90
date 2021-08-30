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

Subroutine LinElas_beta_update(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, & 
								Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
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
	real*8 , allocatable ::  SolU(:), SolU0(:), Xel(:), SolDeltaU(:)  !! all have dimension NodG*NdimE
	integer :: pOrder , NGP, NodG, iFEMtype , iFtype, iMaterial, iDamageParType, iSimplex, iBubble
	integer :: iShiftU, iShiftBeta, iShiftDeltaU, iShiftDeltaBeta, iShiftInfoCrack, SolitaryNode
	type(ptGaussClass) :: PtG
	real*8 :: beta(NdimE), Dts(NdimE,NdimE), dbeta(NdimE)
	real*8 :: Auu(3*NdimE,3*NdimE), Aub(3*NdimE,NdimE), Abu(NdimE,3*NdimE), Abb(NdimE,NdimE), &
			Bu(3*NdimE), Bb(NdimE), AbbInv(NdimE,NdimE)
	real*8 :: cTs, aTs, damageTs, AreaS, normal(NdimE), detAbb !(dummy)
	real*8 :: Bmat(NdimE*NdimE-1,3*NdimE), Dmat(NdimE*NdimE-1,NdimE*NdimE-1), & 
			  BmatT(3*NdimE,NdimE*NdimE-1), BmatBeta(NdimE*NdimE-1,NdimE), BmatBetaT(NdimE,NdimE*NdimE-1) 
	integer , parameter :: NodBeta = 1
	integer :: flagBif, isEFEM
	logical :: isBif
	real*8 :: normBeta,  nu, E

	iShiftU = nint(CommonPar(1))
	iShiftDeltaU = nint(CommonPar(2))
	iFEMtype =  nint(commonPar(3)) 
	iMaterial = nint(commonPar(4))
	isEFEM = nint(commonPar(5))
	aTs = commonPar(6)
	cTs = commonPar(7)
	iShiftInfoCrack = nint(commonPar(8)) 
	
	flagBif = nint(Param(LengthParam + iShiftInfoCrack + 4))
	
	if(flagBif > 0 .and. isEFEM > 0) then !! update is trivial when not bifurcated
	
		call getMaterial(constLaw, matPar, iMaterial)
		call setFEMtype(iFEMtype,NodG,pOrder,NGP,iSimplex,iBubble)

		call getSliceAllocate(SolU,Sol1,1,NodG,iShiftU + 1 ,iShiftU + NdimE, iDofT)
		call getSliceAllocate(SolU0,Sol0,1,NodG,iShiftU + 1 ,iShiftU + NdimE, iDofT)
		call getSliceAllocate(SolDeltaU,Sol1,1,NodG,iShiftDeltaU + 1 ,iShiftDeltaU + NdimE, iDofT)
		call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
		
		call PtG%init(Xel,NodG,NdimE,NGP,pOrder,iBubble, iSimplex)

		AreaS = 0.5d0
		solitaryNode = nint(Param(iShiftInfoCrack + 3))
		normal(1) = 1.0d0
		normal(2) = 0.0d0

		E = matpar(1)
		nu = matpar(2)

		call buildElasticityTensorMatSym(Dmat,E,nu,NdimE)
				
		beta(1:NdimE) = Param(iShiftInfoCrack + 1 : iShiftInfoCrack + NdimE) !!! in the first node, temporary
		normBeta = norm2(beta)
	
		damageTs = aTs*normBeta

		Dts = -aTs*cTs*deltaKron(1:NdimE,1:NdimE)			

		AE = 0.0d0
		BE = 0.0d0
		Bu = 0.0d0
		Bb = 0.0d0
		Auu = 0.0d0
		Aub = 0.0d0
		Abu = 0.0d0
		Abb = 0.0d0

		
		do nG = 1, NGP ! LoopGauss
		
			call PtG%getBmatSym(Bmat,nG)
			BmatT = transpose(Bmat)				
			
			BmatBeta = Bmat(:,5:6) !!! just the last node
			BmatBetaT = transpose(BmatBeta) !!! just the last node
			
			dV=ptG%dV(nG)
												
			Aub = -matmul(BmatT,matmul(Dmat,BmatBeta))*dV
			Abu = transpose(Aub)
				
			Bb = 0.0d0 
			Abb = matmul(BmatBetaT,matmul(Dmat,BmatBeta))*dV + Dts*AreaS
			
		EndDo !LoopGauss
		
		call matInv(AbbInv,detAbb,Abb) 	

		dbeta = matmul(AbbInv,Bb - matmul(Abu,SolDeltaU))
		beta = beta + dbeta	
		Param(LengthParam + iShiftInfoCrack + 1 : LengthParam + iShiftInfoCrack + NdimE) =  beta(1:NdimE)

	end if
	
	AE(1,1) = 1.0d0
	BE(1) = Sol1(1)
	
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine LinElas_beta_updateS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
    use ptsGaussLib2
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)

    Coupling(1,1) = 1
    
end Subroutine


