Subroutine FSFbar(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
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
         
	Integer :: nG, i,j,k,l,m,n,p,q, A, B, ApRow, BpCol , App,Bqp, ipDim, ipDimE, ipDof,ipDofmod ! counters 
	Integer :: constLaw ! if NparamDamage = 0 , no damage
	Real*8 :: PhiA,PhiB, dV , Det , bm(NdimE), MatPar(maxMatPar), damagePar(maxDamagePar)
	Real*8 :: dPhiA(NdimE), dPhiB(NdimE), dPhiB0(NdimE), detF0, detF
	Real*8 , dimension(NdimE,NdimE) :: GradU , DGradU, Ppk, F, F0, Fbar
	Real*8 , dimension(NdimE,NdimE,NdimE,NdimE) :: D ,  Qstar, Q0star
	real*8 , allocatable ::  SolU(:) , Xel(:) !! all have dimension NodG*NdimE
	integer :: pOrder , NGP, NodG, iFEMtype , iFtype, iMaterial, iDamageParType, iSimplex, iBubble
	integer :: iShiftU, iShiftDeltaU
	type(ptGaussClass) :: PtG, PtG0
	
	iShiftU = nint(CommonPar(1))
	iShiftDeltaU = nint(CommonPar(2))
	iFEMtype =  nint(commonPar(3)) 
	iFtype = nint(commonPar(4)) 
	iMaterial = nint(CommonPar(5))
	iDamageParType = nint(commonPar(6))
	
	if(iDamageParType>0) call getDamagePar(damagePar,iDamageParType)
	call getMaterial(constLaw, matPar, iMaterial)
	call setFEMtype(iFEMtype,NodG,pOrder,NGP,iSimplex,iBubble)

	call getSliceAllocate(SolU,Sol1,1,NodG,iShiftU + 1 ,iShiftU + NdimE, iDofT)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	
	call PtG%init(Xel,NodG,NdimE,NGP,pOrder,iBubble, iSimplex)
		
	nG = 1
	call PtG0%init(Xel,NodG,NdimE,1,pOrder,iBubble,iSimplex)
	call PtG0%calcGradU(GradU,SolU,nG)
	call getF(F0,iFtype) !!! with identity summed
	F0 = F0 + GradU
	call calcI3(detF0,F0)
	
	Do nG = 1, NGP ! LoopGauss
	
		GradU = 0.0d0 
		
		call PtG%calcGradU(GradU,SolU,nG)
		
		call getF(F,iFtype) !!! with identity summed
		
		F = F + GradU
		
		call calcI3(detF,F)
		
		Fbar = ((detF0/detF)**(1.0d0/2.0d0))*F !! for 2D
		
		call calcPpk(Ppk,Fbar,NdimE,matPar,constLaw)		
		call calcD(D,Fbar,NdimE,matPar,constLaw)	

		!! if necessary depending on the integration scheme
!~ 		call damageUpdateEquilibrium(Fbar,matPar,damagePar,Param,constLaw,nG) 
!~ 		
!~ 		call damageModifyTangent(D,Ppk,damagePar,Param, nG,1)
		
		call buildFbarTensors(D,Ppk,Qstar,Q0star,F,F0, NdimE)
		
		dV=ptG%dV(nG)
		
		Do A=1, NodG !LoopRow
			ApRow  = (A-1)*iDofT + iShiftDeltaU
			PhiA = PtG%Phi(A,nG)
			dPhiA= PtG%dPhi_G(:,A,nG)
			
			Do p = 1 , NdimE
				App = ApRow+p
				
				BE(App) = BE(App) - dot_product(Ppk(p,:),dPhiA)*DV  
				
				do B= 1 , NodG ! LoopCol ! not considering the simmetry
					BpCol  = (B-1)*iDofT + iShiftDeltaU
					PhiB  = PtG%Phi(B,nG)
					dPhiB = PtG%dPhi_G(:,B,nG) 
					dPhiB0 = PtG0%dPhi_G(:,B,1)
			
					do q=1,NdimE
						Bqp = BpCol+q
						
						AE(App,Bqp) = AE(App,Bqp) + DGradPhiBdotGradPhiA(D,dPhiB,dPhiA,p,q,NdimE)*DV 
						AE(App,Bqp) = AE(App,Bqp) + DGradPhiBdotGradPhiA(Q0star,dPhiB0,dPhiA,p,q,NdimE)*DV
						AE(App,Bqp) = AE(App,Bqp) - DGradPhiBdotGradPhiA(Qstar,dPhiB,dPhiA,p,q,NdimE)*DV
					end do ! loop Bq
				Enddo !LoopCol
			end do ! loop Ap
		Enddo !LoopRow
	EndDo !LoopGauss
	
 	
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

end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine FSFbarS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
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


