    !> Adjoint Problem for the damage sensivity analysis
    !!
	!! @param iShiftU = nint(CommonPar(1))
	!! @param iShiftUd = nint(CommonPar(2))
	!! @param iShiftL = nint(CommonPar(3))
	!! @param iFEMtype =  nint(commonPar(4)) 
	!! @param iMaterial = nint(commonPar(5))
	!! @param iDamageParType = nint(commonPar(6))
	!! @param weight1 = CommonPar(7)
	!! @param weight2 = CommonPar(8)
	!! @param normL2_ref = CommonPar(9)
	!! @param normEnergy_ref = CommonPar(10)
    !! @author Rocha, Felipe Figueredo
    

Subroutine AdjointProblem(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
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
	Real*8 :: PhiA,PhiB, dV , Det , MatPar(maxMatPar), damagePar(maxDamagePar)
	Real*8 :: dPhiA(NdimE), dPhiB(NdimE)
	Real*8 :: GradU(NdimE,NdimE) , DGradU(NdimE,NdimE), F(NdimE,NdimE)
	Real*8 :: D(NdimE,NdimE,NdimE,NdimE), PpkDif(NdimE,NdimE), SolDifU(6), Udif(NdimE), U(NdimE)
	real*8 , allocatable ::  SolU(:), SolUd(:) , Xel(:) !! all have dimension NodG*NdimE
	integer :: pOrder , NGP, NodG, iFEMtype , iMaterial, iDamageParType, iSimplex, iBubble
	integer :: iShiftU, iShiftUd, iShiftL
	type(ptGaussClass) :: PtG
	
	real*8 :: weight1, weight2, normL2_ref, normEnergy_ref, weight3
	
	iShiftU = nint(CommonPar(1))
	iShiftUd = nint(CommonPar(2))
	iShiftL = nint(CommonPar(3))
	iFEMtype =  nint(commonPar(4)) 
	iMaterial = nint(commonPar(5))
	iDamageParType = nint(commonPar(6))
	weight1 = CommonPar(7)
	weight2 = CommonPar(8)
	normL2_ref = CommonPar(9)
	normEnergy_ref = CommonPar(10)
	weight3 = CommonPar(11)
	
	weight1 = weight1/normL2_ref
	weight2 = weight2/normEnergy_ref

	if(iDamageParType>0) call getDamagePar(damagePar,iDamageParType)
	call getMaterial(constLaw, matPar, iMaterial)
	call setFEMtype(iFEMtype,NodG,pOrder,NGP,iSimplex,iBubble)

	call getSliceAllocate(SolU,Sol1,1,NodG,iShiftU + 1 ,iShiftU + NdimE, iDofT)
	call getSliceAllocate(SolUd,Sol1,1,NodG,iShiftUd + 1 ,iShiftUd + NdimE, iDofT)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	
	SolDifU = SolUd - SolU !! ud - u formulation
	

	
	call PtG%init(Xel,NodG,NdimE,NGP,pOrder,iBubble, iSimplex)
	
	Do nG = 1, NGP ! LoopGauss
		
		!!! SolU
		GradU = 0.0d0
		call PtG%calcGradU(GradU,SolU,nG)
		call PtG%calcU(U,SolU,nG)
		
		F = deltaKron(1:NdimE,1:NdimE) + GradU
		
		call calcD(D,F,NdimE,matPar,constLaw)	
		
		call damageModifyTangent(D,PpkDif,damagePar,Param, nG, 1) !!! PpkDif is dummy here
		
		!!! SolDifU
		GradU = 0.0d0
		call PtG%calcGradU(GradU,SolDifU,nG)
		call PtG%calcU(Udif,SolDifU,nG)
		F = deltaKron(1:NdimE,1:NdimE) + GradU
		
		call calcPpk(PpkDif,F,NdimE,matPar,constLaw)
											
		dV=ptG%dV(nG)
		
		Do A=1, NodG !LoopRow
			ApRow  = (A-1)*iDofT + iShiftL
			PhiA = PtG%Phi(A,nG)
			dPhiA= PtG%dPhi_G(:,A,nG)
			ipDim = (A-1)*NdimE
			
			Do p = 1 , NdimE
				App = ApRow+p
				  
				BE(App) = BE(App) + 2.0d0 * weight1 * Udif(p) * PhiA * DV
				BE(App) = BE(App) + weight2 * dot_product(PpkDif(p,:),dPhiA) * DV
				BE(App) = BE(App) - 2.0d0 * weight3 * U(p) * PhiA * DV
				
				
!~ 				write(0,*) BE(App)

				do B= 1 , NodG ! LoopCol ! not considering the simmetry
					BpCol  = (B-1)*iDofT + iShiftL
					PhiB  = PtG%Phi(B,nG)
					dPhiB = PtG%dPhi_G(:,B,nG) 
 					
					do q=1,NdimE
						Bqp = BpCol+q
						
!~ 						AE(App,Bqp) = AE(App,Bqp) + DGradPhiBdotGradPhiA(D,dPhiB,dPhiA,p,q,NdimE)*DV  
						AE(App,Bqp) = AE(App,Bqp) + DGradPhiBdotGradPhiA(D,dPhiA,dPhiB,q,p,NdimE)*DV !! have to be changed
				
					end do ! loop Bq
				Enddo !LoopCol
			end do ! loop Ap
		Enddo !LoopRow
	EndDo !LoopGauss
	
	deallocate(SolU,Xel,SolUd)
	
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine AdjointProblemS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
    use ptsGaussLib2
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
    
    Integer :: A,B,q,p, ApRow, BpCol,iShiftL , Femtype, NodG
    integer , parameter :: NdimE = 2

    Coupling = 0
    
    iShiftL = nint(commonPar(3))
    femType = nint(commonPar(4)) 
    
    call setNodG(femtype, NodG)
    
    do A = 1,NodG
        ApRow = (A-1) * iDofT + iShiftL
		
		do B = 1,NodG
			BpCol = (B-1) * iDofT + iShiftL

			do p=1,nDimE
			do q=1,nDimE

				Coupling (ApRow + p, BpCol + q) = 1

			enddo
            enddo
		enddo
	enddo
 	
end Subroutine
