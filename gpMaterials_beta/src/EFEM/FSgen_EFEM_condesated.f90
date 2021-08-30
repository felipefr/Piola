    !> Generic element for finite strains with discontinuity through EFEM (with no static condesation for beta)
    !!
	!! @param iShiftU = nint(CommonPar(1))
	!! @param iShiftDeltaU = nint(CommonPar(2))
	!! @param iFemType  = nint(CommonPar(3)) 
	!! @param iFType  = nint(commonPar(4))
	!! @param iMaterial  = nint(commonPar(5))
	!! @param iDamageParType  = nint(commonPar(6))
	!!
    !! @author Rocha, Felipe Figueredo

Subroutine FSgen_EFEM_condesated(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
	use globalVariables, only : getF, NdimE, getDamagePar, getMaterial, maxDamagePar, maxMatPar
	use finiteStrainLib
	use damageNewLib
	use ptsGaussLib2
	use efemLib
	
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT, NodElt! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

!~ 	!   =====   END ARGUMENTS  =======
         
	Integer :: A,B,q,p, ApRow, BpCol, App, nG ! counters 
	Integer :: constLaw ! if NparamDamage = 0 , no damage
	Real*8 :: dV , Det , MatPar(maxMatPar), damagePar(maxDamagePar)
	Real*8 :: GradU(NdimE,NdimE), F(NdimE,NdimE)
	Real*8 :: D(NdimE,NdimE,NdimE,NdimE), Ppk(NdimE,NdimE)
	real*8 , allocatable ::  SolU(:), SolDeltaU(:), Xel(:), SolU0(:)  !! all have dimension NodG*NdimE
	integer :: pOrder , NGP, NodG, iFEMtype , iFtype, iMaterial, iDamageParType, iSimplex, iBubble
	integer :: iShiftU, iShiftDeltaU, iShiftInfoCrack, SolitaryNode
	type(ptGaussClass) :: PtG
	real*8 :: beta(NdimE), beta0(NdimE), dVarPhi(NdimE), beta_Ten_dVarPhi(NdimE,NdimE), Ts(NdimE), Dts(NdimE,NdimE)
	real*8 :: Auu(3*NdimE,3*NdimE), Aub(3*NdimE,NdimE), Abu(NdimE,3*NdimE), & 
				Abb(NdimE,NdimE), Bu (3*NdimE), Bb(NdimE), AbbInv(NdimE,NdimE), detAbb
	real*8 :: cTs, aTs, damageTs, AreaS
	real*8 :: Bmat(NdimE*NdimE,3*NdimE), Dmat(NdimE*NdimE,NdimE*NdimE) , Pmat(NdimE*NdimE), & 
			  BmatT(3*NdimE,NdimE*NdimE), BmatBeta(NdimE*NdimE,NdimE), BmatBetaT(NdimE,NdimE*NdimE) , normal(NdimE)
	integer :: flagBif
	logical :: isBif
	
	iShiftU = nint(CommonPar(1))
	iShiftDeltaU = nint(CommonPar(2))
	iFEMtype =  nint(commonPar(3)) 
	iFtype = nint(commonPar(4)) 
	iMaterial = nint(commonPar(5))
	iDamageParType = nint(commonPar(6))
	aTs = nint(commonPar(7))
	cTs = nint(commonPar(8))
	iShiftInfoCrack = nint(CommonPar(9))

	if(iDamageParType>0) call getDamagePar(damagePar,iDamageParType)
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
			
		beta(1:NdimE) = Param(LengthParam + iShiftInfoCrack + 1 : LengthParam + iShiftInfoCrack + NdimE) 
		beta0(1:NdimE) = Param(iShiftInfoCrack + 2 : iShiftInfoCrack + NdimE + 1)
		damageTs = 1.0d0 - dexp(-aTs*norm2(beta0))
! 		damageTs = aTs*norm2(beta0)

		if(flagBif == 1) then !! that is, we have a previous step of bifurcation
			call posProcessBeta(Beta,PtG,Param,matPar,damagePar,SolU,SolDeltaU, cTs, aTs, damageTs, & 
						AreaS, constLaw, iFtype,SolitaryNode, NodG, NGP, iDamageParType)
		end if
		
		Param(LengthParam + iShiftInfoCrack + 1 : LengthParam + iShiftInfoCrack + NdimE) =  beta(1:NdimE)
		Ts = (1.0d0 - damageTs)*cTs*beta
		Dts = (1.0d0 - damageTs)*cTs*deltaKron(1:NdimE,1:NdimE)
	
	end if

	
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
!~ 		
		if(iDamageParType>0) then
			call damageUpdateEquilibrium(F,matPar,damagePar,Param,constLaw,nG) !~ 		
			call damageModifyTangent(D,Ppk,damagePar,Param, nG, 1)
		end if						
						
		call PtG%getBmat(Bmat,nG)
		BmatT = transpose(Bmat)				
		call voigtTen4toTen2(Dmat,D)
		call voigtTen2toVec(Pmat,Ppk)

		dV=ptG%dV(nG)
								
		Bu = - matmul(BmatT,Pmat)*dV
		Auu = matmul(BmatT,matmul(Dmat,Bmat))*dV

		if(flagBif > 0) then
					
			BmatBeta = Bmat(:,(solitaryNode-1)*NdimE + 1:(solitaryNode-1)*NdimE + 2) 
			BmatBetaT = transpose(BmatBeta) 
						
			Aub = -matmul(BmatT,matmul(Dmat,BmatBeta))*dV
			Abu = transpose(Aub)
			
			Bb = matmul(BmatBetaT,Pmat)*dV - Ts*AreaS 
			Abb = matmul(BmatBetaT,matmul(Dmat,BmatBeta))*dV + Dts*AreaS
		end if
 		
	EndDo !LoopGauss
 	
	if(iShiftU == iShiftDeltaU) then !!! I dont know if it is before or after the if below
		Bu = Bu + matmul(Auu,SolU) 
!~ 		+ matmul(Aub,Beta)
	end if

 	if(flagBif > 0) then
		call matInv(AbbInv,detAbb,Abb) 	
		Auu = Auu - matmul(Aub,matmul(AbbInv,Abu))
		Bu = Bu - matmul(Aub,matmul(AbbInv,Bb))		
	end if
	
	call setBlockToVector(BE,Bu,NodG,iDofT, NdimE, iShiftDeltaU ) !(A,B, NodG, iDofA, iDofB, iShiftA) , B in A
	
	call setBlockToMatrixSquare(AE,Auu, NodG,iDofT, NdimE, iShiftDeltaU ) !(A,B, NodG, iDofA, iDofB, iShiftA) , B in A
	
	if(iShiftU /= iShiftDeltaU) then
		call setBlockDiagonalConstantToMatrix(AE,1.0d0,NodG, iDofT, NdimE, iShiftU, iShiftU)
		call setBlockDiagonalConstantToMatrix(AE,-1.0d0,NodG, iDofT, NdimE, iShiftU, iShiftDeltaU)
		call setBlockToVector(BE,SolU,NodG,iDofT, NdimE, iShiftU ) !(A,B, NodG, iDofA, iDofB, iShiftA) , B in A
	end if

end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine FSgen_EFEM_condesatedS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
    use ptsGaussLib2
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
    
    Integer :: A,B,q,p, ApRow, BpCol,iShiftU , iShiftDeltaU, Femtype, NodG
    integer , parameter :: NdimE = 2

    Coupling = 0
    
    iShiftU = nint(commonPar(1))
    iShiftDeltaU = nint(commonPar(2))
    femType = nint(commonPar(3)) 
    
    call setNodG(femtype, NodG)
        
	call setSimbolicBlockToMatrixSquare(Coupling,NodG,iDofT, NdimE, iShiftDeltaU )

	if(iShiftU /= iShiftDeltaU) then
		call setSimbolicBlockDiagonalToMatrix(Coupling,NodG, iDofT, NdimE, iShiftU, iShiftU)
		call setSimbolicBlockDiagonalToMatrix(Coupling,NodG, iDofT, NdimE, iShiftU, iShiftDeltaU)
	end if

end Subroutine



