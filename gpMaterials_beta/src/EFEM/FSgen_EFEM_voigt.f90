    !> Generic element for finite strains with discontinuity through EFEM (with no static condesation for beta)
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

Subroutine FSgen_EFEM_voigt(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
	use globalVariables, only : getF, NdimE, getDamagePar, getMaterial, maxDamagePar, maxMatPar
	use finiteStrainLib
	use damageNewLib
	use ptsGaussLib2
	
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT, NodElt! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

!~ 	!   =====   END ARGUMENTS  =======
         
	Integer :: nG, i,j,k,l,m,n,p,q, A, B, ApRow, BpCol , App,Bqp, ApDim, BpDim, ipDim, ipDimE, ipDof,ipDofmod ! counters 
	Integer :: constLaw ! if NparamDamage = 0 , no damage
	Real*8 :: PhiA,PhiB, dV , Det , MatPar(maxMatPar), damagePar(maxDamagePar)
	Real*8 :: dPhiA(NdimE), dPhiB(NdimE)
	Real*8 :: GradU(NdimE,NdimE) , DGradU(NdimE,NdimE), F(NdimE,NdimE)
	Real*8 :: D(NdimE,NdimE,NdimE,NdimE), Ppk(NdimE,NdimE)
	real*8 , allocatable ::  SolU(:), SolU0(:), Xel(:)  !! all have dimension NodG*NdimE
	integer :: pOrder , NGP, NodG, iFEMtype , iFtype, iMaterial, iDamageParType, iSimplex, iBubble
	integer :: iShiftU, iShiftDeltaU, iShiftBeta, iShiftDeltaBeta, iShiftInfoCrack, SolitaryNode
	type(ptGaussClass) :: PtG
	real*8 :: beta(NdimE), beta0(NdimE), dVarPhi(NdimE), beta_Ten_dVarPhi(NdimE,NdimE), Ts(NdimE), Dts(NdimE,NdimE)
	real*8 :: Auu(3*NdimE,3*NdimE), Aub(3*NdimE,NdimE), Abu(NdimE,3*NdimE), Abb(NdimE,NdimE), Bu (3*NdimE), Bb(NdimE), normal(NdimE)
	real*8 :: cTs, aTs, damageTs, AreaS
	real*8 :: Bmat(NdimE*NdimE,3*NdimE), Dmat(NdimE*NdimE,NdimE*NdimE) , Pmat(NdimE*NdimE), & 
			  BmatT(3*NdimE,NdimE*NdimE), BmatBeta(NdimE*NdimE,NdimE), BmatBetaT(NdimE,NdimE*NdimE) 
	integer , parameter :: NodBeta = 1
	integer :: flagBif
	logical :: isBif
	real*8 :: normal_beta(NdimE,NdimE), normBeta, der_damage
	
	
	iShiftU = nint(CommonPar(1))
	iShiftDeltaU = nint(CommonPar(2))
	iShiftBeta = nint(CommonPar(3))
	iShiftDeltaBeta = nint(CommonPar(4))
	iFEMtype =  nint(commonPar(5)) 
	iFtype = nint(commonPar(6)) 
	iMaterial = nint(commonPar(7))
	iDamageParType = nint(commonPar(8))
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
	
	if(flagBif < 1 ) then
			
		call PtG%calcGradU(GradU,SolU0,1)
		call getF(F,iFtype) !!! with identity summed
		F = F + GradU
		call calcPpk(Ppk,F,NdimE,matPar,constLaw)
		
		if(Ppk(1,1) > cTs) then
			flagBif = 2  !!! it means that the bifurcation was detected at this step
			Param(LengthParam + iShiftInfoCrack + 4) = 1
		end if
	
	end if

	if(flagBif > 0 ) then
			
		beta(1:NdimE) = Sol1(iShiftBeta + 1 : iShiftBeta + NdimE) !!! in the first node, temporary
		beta0(1:NdimE) = Sol0(iShiftBeta + 1 : iShiftBeta + NdimE) !!! in the first node, temporary
		
		normBeta = norm2(beta)
		
		damageTs = 1.0d0 - dexp(-aTs*normBeta)
		der_damage = aTs*dexp(-aTs*normBeta)
		
		call VotimesW(normal_beta,normal,beta)
		Ts = (1.0d0 - damageTs)*cTs*normal
		Dts = - (cTs*der_damage/normBeta)*normal_beta
		
	else
		beta = 0.0d0
		beta0 = 0.0d0
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
		
		call getF(F,iFtype) !!! with identity summed
		
		F = F + GradU 
		
		if(flagBif > 0) then
			dVarPhi = PtG%dPhi_G(:,solitaryNode,nG) !! solitary node is the last one, convention
			call VotimesW(beta_Ten_dVarPhi,beta,dVarPhi)  
			F = F - beta_Ten_dVarPhi
		end if
				
		call calcPpk(Ppk,F,NdimE,matPar,constLaw)		
		call calcD(D,F,NdimE,matPar,constLaw)	
		
		if(iDamageParType>0) then
			call damageUpdateEquilibrium(F,matPar,damagePar,Param,constLaw,nG) !~ 		
			call damageModifyTangent(D,Ppk,damagePar,Param, nG, 1)
		end if						
						
		call PtG%getBmat(Bmat,nG)
		BmatT = transpose(Bmat)				
		call voigtTen4toTen2(Dmat,D)
		call voigtTen2toVec(Pmat,Ppk)
		
		BmatBeta = Bmat(:,5:6) !!! just the last node
		BmatBetaT = transpose(BmatBeta) !!! just the last node
		
		dV=ptG%dV(nG)
								
		Bu = - matmul(BmatT,Pmat)*dV
		Auu = matmul(BmatT,matmul(Dmat,Bmat))*dV

		if(flagBif > 0) then
			Aub = -matmul(BmatT,matmul(Dmat,BmatBeta))*dV
			Abu = transpose(Aub)
			
			Bb = matmul(BmatBetaT,Pmat)*dV - Ts*AreaS 
			Abb = matmul(BmatBetaT,matmul(Dmat,BmatBeta))*dV + Dts*AreaS
		else
			Abb = deltaKron(1:NdimE,1:NdimE)
		end if
 		
	EndDo !LoopGauss
 	
 	if(iShiftU == iShiftDeltaU .and. iShiftBeta == iShiftDeltaBeta) then
		Bu = Bu + matmul(Auu,SolU) + matmul(Aub,Beta)
		Bb = Bb + matmul(Abu,SolU) + matmul(Abb,Beta) 
	end if

	call setBlockToVector(BE,Bu,NodG,iDofT, NdimE, iShiftDeltaU ) !(A,B, NodG, iDofA, iDofB, iShiftA) , B in A
	call setBlockToVector(BE,Bb,1,iDofT, NdimE, iShiftDeltaBeta ) !(A,B, NodG, iDofA, iDofB, iShiftA) , B in A
	
	call setBlockToMatrixSquare(AE,Auu, NodG,iDofT, NdimE, iShiftDeltaU ) !(A,B, NodG, iDofA, iDofB, iShiftA) , B in A
	call setBlockToMatrixSquare(AE,Abb, 1,iDofT, NdimE, iShiftDeltaBeta ) !(A,B, NodG, iDofA, iDofB, iShiftA) , B in A
	
	call setBlockToMatrix(AE,Aub, NodG,1, iDofT,iDofT, NdimE,NdimE, iShiftDeltaU, iShiftDeltaBeta) ! (A,B, NodGi,NodGj, iDofAi,iDofAj, iDofBi,iDofBj, iShiftAi, iShiftAj)
	call setBlockToMatrix(AE,Abu, 1, NodG, iDofT, iDofT, NdimE, NdimE, iShiftDeltaBeta, iShiftDeltaU ) ! (A,B, NodGi,NodGj, iDofAi,iDofAj, iDofBi,iDofBj, iShiftAi, iShiftAj)
	
	if(iShiftU /= iShiftDeltaU) then
		call setBlockDiagonalConstantToMatrix(AE,1.0d0,NodG, iDofT, NdimE, iShiftU, iShiftU)
		call setBlockDiagonalConstantToMatrix(AE,-1.0d0,NodG, iDofT, NdimE, iShiftU, iShiftDeltaU)
		call setBlockToVector(BE,SolU,NodG,iDofT, NdimE, iShiftU ) !(A,B, NodG, iDofA, iDofB, iShiftA) , B in A
		call setBlockDiagonalConstantToMatrix(AE,1.0d0,1, iDofT, NdimE, iShiftBeta, iShiftBeta)
		call setBlockDiagonalConstantToMatrix(AE,-1.0d0,1, iDofT, NdimE, iShiftBeta, iShiftDeltaBeta)
		call setBlockToVector(BE,Beta,1,iDofT, NdimE, iShiftBeta ) !(A,B, NodG, iDofA, iDofB, iShiftA) , B in A
	end if
	
	
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine FSgen_EFEM_voigtS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
    use ptsGaussLib2
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
    
    Integer :: A,B,q,p, ApRow, BpCol,iShiftU , iShiftDeltaU, Femtype, NodG, iShiftDeltaBeta, iShiftBeta
    integer , parameter :: NdimE = 2

    Coupling = 0
    
    iShiftU = nint(commonPar(1))
    iShiftDeltaU = nint(commonPar(2))
    iShiftBeta = nint(commonPar(3))
    iShiftDeltaBeta = nint(commonPar(4))
    femType = nint(commonPar(5)) 
    
    call setNodG(femtype, NodG)
    
    call setSimbolicBlockToMatrixSquare(Coupling,NodG,iDofT, NdimE, iShiftDeltaU )
	call setSimbolicBlockToMatrixSquare(Coupling,1,iDofT, NdimE, iShiftDeltaBeta ) 
	
	call setSimbolicBlockToMatrix(Coupling,NodG,1, iDofT,iDofT, NdimE,NdimE, iShiftDeltaU, iShiftDeltaBeta) 
	call setSimbolicBlockToMatrix(Coupling,1, NodG, iDofT, iDofT, NdimE, NdimE, iShiftDeltaBeta, iShiftDeltaU ) 
	
	
	if(iShiftU /= iShiftDeltaU) then
		call setSimbolicBlockDiagonalToMatrix(Coupling,NodG, iDofT, NdimE, iShiftU, iShiftU)
		call setSimbolicBlockDiagonalToMatrix(Coupling,NodG, iDofT, NdimE, iShiftU, iShiftDeltaU)
		call setSimbolicBlockDiagonalToMatrix(Coupling, 1, iDofT, NdimE, iShiftBeta, iShiftBeta)
		call setSimbolicBlockDiagonalToMatrix(Coupling,1, iDofT, NdimE, iShiftBeta, iShiftDeltaBeta)
	end if
	
end Subroutine
