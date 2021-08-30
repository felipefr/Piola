Subroutine finiteStrainSpatial(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
	use loadingLib
	use finiteStrainLib
	use damageLib
	use ptsGaussLib2
	
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT, NodElt! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

!~ 	!   =====   END ARGUMENTS  =======
         
	Integer :: nG, i,j,k,l,m,n,p,q, A, B, ApRow, BpCol , App,Bqp, ipDim, ipDimE, ipDof,ipDofmod ! counters 
	Integer :: constLaw ! if NparamDamage = 0 , no damage
	integer , parameter ::  NdimE = 2! this code only makes sense in 2D and for triangles 
	Real*8 :: PhiA,PhiB, dV , Det , bm(NdimE), MatPar(9)
	Real*8 :: dPhiA(NdimE), dPhiB(NdimE), detFdummy, detFmu
	Real*8 , dimension(NdimE,NdimE) :: GradU , DGradU, F, Finv, Fmu, FmuInv
	Real*8 :: Ds(NdimE,NdimE,NdimE,NdimE), sigma(NdimE,NdimE)
	real*8 , allocatable ::  SolU(:) , Xel(:) !! all have dimension NodG*NdimE
	integer :: pOrder , NGP, NodG, iFEMtype , iSimplex, iBubble
	Real*8 :: LoadPar(6) , he, energy
	integer :: iShiftU, iShiftDeltaU , LoadType , LoadTypeProg
	!! commonPar for the damage starts after the commonPar of the actual function
	integer , parameter :: iShiftDamageParam = 0 , ishiftDamageCommonPar = 16
	integer :: IdamageModify
	type(ptGaussClass) :: PtG
	
	iShiftU = nint(CommonPar(1))
	iShiftDeltaU = nint(CommonPar(2))
	iFEMtype =  nint(commonPar(3)) 
	constLaw = nint(CommonPar(4))
	matPar(1:3) = commonPar(5:7)
	LoadType = nint(commonPar(8))
	LoadTypeProg = nint(commonPar(9))
	LoadPar(:) = CommonPar(10:15)
	IdamageModify = nint(commonPar(16))
	
	call setFEMtype(iFEMtype,NodG,pOrder,NGP,iSimplex,iBubble)

	call getSliceAllocate(SolU,Sol1,1,NodG,iShiftU + 1 ,iShiftU + NdimE, iDofT)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	
	Xel = Xel + SolU !! Xs
	
	call PtG%init(Xel,NodG,NdimE,NGP,pOrder,iBubble, iSimplex)
	
	Do nG = 1, NGP ! LoopGauss
		
		call PtG%calcGradU(GradU,SolU,nG)
		
!~ 		gradU(1,1) = 1.1d0 
!~ 		gradU(1,2) = 0.1d0
!~ 		gradU(2,1) = 0.3d0
!~ 		gradU(2,2) = 0.99d0
		
		call setF(F,LoadPar,Time,DelT,LoadTypeProg,LoadType)
		
		F = F + deltaKron(1:NdimE,1:NdimE)  
		
		call matInv(Finv,detFdummy,F)
		
		FmuInv = Finv - GradU

		call matInv(Fmu,detFmu,FmuInv)
		detFmu = 1.0d0/detFmu
		
		call calcPpk(sigma,Fmu,NdimE,matPar,constLaw)		
		call calcD(Ds,Fmu,NdimE,matPar,constLaw)	
		call modifyToSpatial(Ds,sigma,Fmu,detFmu,NdimE)
		
!~ 		write(0,*) "Q1"
!~ 		call numprint(Ds(:,:,1,1))
!~ 		call numprint(Ds(:,:,1,2))
!~ 		call numprint(Ds(:,:,2,1))
!~ 		call numprint(Ds(:,:,2,2))
!~  		pause
									
		dV=ptG%dV(nG)
		
		if(iShiftU == iShiftDeltaU) call T4xT2(DGradU,Ds,GradU)
		
		Do A=1, NodG !LoopRow
			ApRow  = (A-1)*iDofT + iShiftDeltaU
			PhiA = PtG%Phi(A,nG)
			dPhiA= PtG%dPhi_G(:,A,nG)
			
			
			Do p = 1 , NdimE
				App = ApRow+p
				
				BE(App) = BE(App) - dot_product(sigma(p,:),dPhiA)*DV  
				if(iShiftU == iShiftDeltaU) BE(App) = BE(App) + dot_product(DGradU(p,:),dphiA)*DV 
				
				do B= 1 , NodG ! LoopCol ! not considering the simmetry
					BpCol  = (B-1)*iDofT + iShiftDeltaU
					PhiB  = PtG%Phi(B,nG)
					dPhiB = PtG%dPhi_G(:,B,nG) 
 					
					do q=1,NdimE
						Bqp = BpCol+q
						
						AE(App,Bqp) = AE(App,Bqp) + DGradPhiBdotGradPhiA(Ds,dPhiB,dPhiA,p,q,NdimE)*DV  
!~ 						write(0,*) AE(App,Bqp) 
					end do ! loop Bq
				Enddo !LoopCol
			end do ! loop Ap
		Enddo !LoopRow
	EndDo !LoopGauss

!~ 	pause

 	if(iShiftU /= iShiftDeltaU) then
		Do A=1, NodG !LoopRow
			ApRow  = (A-1)*iDofT + iShiftU
			
			Do p = 1 , NdimE
				App = ApRow+p
				
				BE(App) = SolU((A-1)*NdimE + p)  				  
				
				AE(App,App) = 1.0d0
				AE(App,App + iShiftDeltaU - iShiftU) = -1.0d0
				
			end do ! loop Ap
		end do !LoopRow
	end if
	
!~ 			pause

	
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine finiteStrainSpatialS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
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

 	if(iShiftU /= iShiftDeltaU) then
		do A = 1,NodG
			ApRow = (A-1) * iDofT + iShiftU
			do p=1,nDimE
				Coupling (ApRow + p, ApRow + p) = 1
				Coupling (ApRow + p, ApRow + p + iShiftDeltaU - iShiftU) = 1
			enddo
		enddo
	end if

 	
	call numPrint(Coupling)

end Subroutine
