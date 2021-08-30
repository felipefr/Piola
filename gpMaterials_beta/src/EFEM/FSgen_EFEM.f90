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

Subroutine FSgen_EFEM(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
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
	real*8 , allocatable ::  SolU(:) , Xel(:)  !! all have dimension NodG*NdimE
	integer :: pOrder , NGP, NodG, iFEMtype , iFtype, iMaterial, iDamageParType, iSimplex, iBubble
	integer :: iShiftU, iShiftDeltaU, iShiftBeta, iShiftDeltaBeta
	type(ptGaussClass) :: PtG
	real*8 :: beta(NdimE), beta0(NdimE), dVarPhi(NdimE), beta_Ten_dVarPhi(NdimE,NdimE), Ts(NdimE), Dts(NdimE,NdimE)
	real*8 :: Auu(3*NdimE,3*NdimE), Aub(3*NdimE,NdimE), Abu(NdimE,3*NdimE), Abb(NdimE,NdimE), Bu (3*NdimE), Bb(NdimE)
	real*8 :: cTs, aTs, damageTs, AreaS
	integer , parameter :: NodBeta = 1
	logical :: isBif
	
	iShiftU = nint(CommonPar(1))
	iShiftDeltaU = nint(CommonPar(2))
	iShiftBeta = nint(CommonPar(3))
	iShiftDeltaBeta = nint(CommonPar(4))
	iFEMtype =  nint(commonPar(5)) 
	iFtype = nint(commonPar(6)) 
	iMaterial = nint(commonPar(7))
	iDamageParType = nint(commonPar(8))

	if(iDamageParType>0) call getDamagePar(damagePar,iDamageParType)
	call getMaterial(constLaw, matPar, iMaterial)
	call setFEMtype(iFEMtype,NodG,pOrder,NGP,iSimplex,iBubble)

	call getSliceAllocate(SolU,Sol1,1,NodG,iShiftU + 1 ,iShiftU + NdimE, iDofT)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	
	aTs = 1.0d0
	cTs = 80.0d0
	AreaS = 0.5d0
	beta(1:NdimE) = Sol1(iShiftBeta + 1 : iShiftBeta + NdimE) !!! in the first node, temporary
	beta0(1:NdimE) = Sol0(iShiftBeta + 1 : iShiftBeta + NdimE) !!! in the first node, temporary
!~ 	damageTs = 1.0d0 - dexp(-aTs*norm2(beta0))
	damageTs = aTs*norm2(beta0)
	
	Ts = (1.0d0 - damageTs)*cTs*beta
	Dts = (1.0d0 - damageTs)*cTs*deltaKron(1:NdimE,1:NdimE)
	

	if(Time > 15.0d0) then
!~ 		write(0,*) "=================== >Time = " , Time
		isBif = .True.
!~ 		pause
	else
		isBif = .False.
	end if
!~ 	
	call PtG%init(Xel,NodG,NdimE,NGP,pOrder,iBubble, iSimplex)
	
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
		
		dVarPhi = PtG%dPhi_G(:,NodG,nG) !! solitary node is the last one, convention
		
		call VotimesW(beta_Ten_dVarPhi,beta,dVarPhi)  
		
		F = F + GradU - beta_Ten_dVarPhi
				
		call calcPpk(Ppk,F,NdimE,matPar,constLaw)		
		call calcD(D,F,NdimE,matPar,constLaw)	
		
		if(iDamageParType>0) then
			call damageUpdateEquilibrium(F,matPar,damagePar,Param,constLaw,nG) !~ 		
			call damageModifyTangent(D,Ppk,damagePar,Param, nG, 1)
		end if						
						
		dV=ptG%dV(nG)
				
		! assembly Auu, and Bu
		Do A=1, NodG !LoopRow
			ApRow  = (A-1)*NdimE
			PhiA = PtG%Phi(A,nG)
			dPhiA= PtG%dPhi_G(:,A,nG)
		
			Do p = 1 , NdimE
				App = ApRow+p
				
				Bu(App) = Bu(App) - dot_product(Ppk(p,:),dPhiA)*DV  
				
				do B= 1 , NodG ! LoopCol ! not considering the simmetry
					BpCol  = (B-1)*NdimE
					PhiB  = PtG%Phi(B,nG)
					dPhiB = PtG%dPhi_G(:,B,nG) 
 					
					do q=1,NdimE
						Bqp = BpCol+q
						
						Auu(App,Bqp) = Auu(App,Bqp) + DGradGrad(D,dPhiA,dPhiB,p,q,NdimE)*DV  
				
					end do ! loop Bq
				Enddo !LoopCol
			end do ! loop Ap
		Enddo !LoopRow


		if(isBif) then
		
			! assembly Abu, Aub and Bb 
			Do p = 1 , NdimE
				App = p
				
				Bb(p) = Bb(p) + dot_product(Ppk(p,:),dVarPhi)*DV  - Ts(p)*AreaS
				
				do B= 1 , NodG ! LoopCol ! not considering the simmetry
					BpCol  = (B-1)*NdimE
					dPhiB = PtG%dPhi_G(:,B,nG) 
					
					do q=1,NdimE
						Bqp = BpCol+q
						
						Abu(p,Bqp) = Abu(p,Bqp) - DGradGrad(D,dVarPhi,dPhiB,p,q,NdimE)*DV  
						
						Aub(Bqp,p) = Abu(p,Bqp)
				
					end do ! loop Bq
				Enddo !LoopCol
			end do ! loop Ap
			
			! assembly Abb
			Do p = 1 , NdimE
				do q=1,NdimE					
					Abb(p,q) = Abb(p,q) + DGradGrad(D,dVarPhi,dVarPhi,p,q,NdimE)*DV + Dts(p,q)*AreaS 
				end do ! q
			end do ! p
			
		else
			Abb = deltaKron(1:NdimE,1:NdimE)
		end if
 		
	EndDo !LoopGauss
 	
	Bu = Bu + matmul(Auu,SolU) + matmul(Aub,Beta)
	Bb = Bb + matmul(Abu,SolU) + matmul(Abb,Beta) 

	call setBlockToVector(BE,Bu,NodG,iDofT, NdimE, iShiftDeltaU ) !(A,B, NodG, iDofA, iDofB, iShiftA) , B in A
	call setBlockToVector(BE,Bb,1,iDofT, NdimE, iShiftDeltaBeta ) !(A,B, NodG, iDofA, iDofB, iShiftA) , B in A
	
	call setBlockToMatrixSquare(AE,Auu, NodG,iDofT, NdimE, iShiftDeltaU ) !(A,B, NodG, iDofA, iDofB, iShiftA) , B in A
	call setBlockToMatrixSquare(AE,Abb, 1,iDofT, NdimE, iShiftDeltaBeta ) !(A,B, NodG, iDofA, iDofB, iShiftA) , B in A
	
	call setBlockToMatrix(AE,Aub, NodG,1, iDofT,iDofT, NdimE,NdimE, iShiftDeltaU, iShiftDeltaBeta) ! (A,B, NodGi,NodGj, iDofAi,iDofAj, iDofBi,iDofBj, iShiftAi, iShiftAj)
	call setBlockToMatrix(AE,Abu, 1, NodG, iDofT, iDofT, NdimE, NdimE, iShiftDeltaBeta, iShiftDeltaU ) ! (A,B, NodGi,NodGj, iDofAi,iDofAj, iDofBi,iDofBj, iShiftAi, iShiftAj)
		
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine FSgen_EFEMS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
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
    
!~     NodG = 6
!~     write(0,*) NodG
!~     pause
    
    do A = 1,NodG
        ApRow = (A-1) * iDofT + iShiftDeltaU
		
		do B = 1,NodG
			BpCol = (B-1) * iDofT + iShiftDeltaU

			do p=1,nDimE
			do q=1,nDimE

				Coupling (ApRow + p, BpCol + q) = 1

			enddo
            enddo
		enddo
	enddo

	 do A = 1,1
			ApRow = (A-1) * iDofT + iShiftDeltaBeta
						
		do B = 1,NodG
			BpCol = (B-1) * iDofT + iShiftDeltaU

			do p=1,nDimE
			do q=1,nDimE

				Coupling (ApRow + p, BpCol + q) = 1

			enddo
			enddo
		enddo
	enddo

    do A = 1,NodG
        ApRow = (A-1) * iDofT + iShiftDeltaU
		
		do B = 1,1
			BpCol = (B-1) * iDofT + iShiftDeltaBeta

			do p=1,nDimE
			do q=1,nDimE

				Coupling (ApRow + p, BpCol + q) = 1

			enddo
            enddo
		enddo
	enddo

    do A = 1,1
        ApRow = (A-1) * iDofT + iShiftDeltaBeta
		
		do B = 1,1
			BpCol = (B-1) * iDofT + iShiftDeltaBeta

			do p=1,nDimE
			do q=1,nDimE

				Coupling (ApRow + p, BpCol + q) = 1

			enddo
            enddo
		enddo
	enddo

	call numPrint(Coupling)

end Subroutine
