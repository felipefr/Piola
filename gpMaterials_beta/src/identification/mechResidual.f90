    !> Compute Mechanical Equilibrium Residual
    !!
	!! @param iShiftU = nint(CommonPar(1))
	!! @param iShiftRes = nint(CommonPar(2))
	!! @param iFemType  = nint(CommonPar(3)) 
	!! @param iFType  = nint(commonPar(4))
	!! @param iMaterial  = nint(commonPar(5))
	!! @param iDamageParType  = nint(commonPar(6))
	!!
    !! @author Rocha, Felipe Figueredo
Subroutine mechResidual(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
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

 	!   =====   END ARGUMENTS  =======
         
	Integer :: nG, i,j,k,l,m,n,p,q, A, B, ApRow, BpCol , App,Bqp, ipDim, ipDimE, ipDof,ipDofmod ! counters 
	Integer :: constLaw ! if NparamDamage = 0 , no damage
	Real*8 :: PhiA,PhiB, dV , Det , MatPar(maxMatPar), damagePar(maxDamagePar)
	Real*8 :: dPhiA(NdimE), dPhiB(NdimE)
	Real*8 :: GradU(NdimE,NdimE) , DGradU(NdimE,NdimE), F(NdimE,NdimE)
	Real*8 :: D(NdimE,NdimE,NdimE,NdimE), Ppk(NdimE,NdimE)
	real*8 , allocatable ::  SolU(:) , Xel(:) !! all have dimension NodG*NdimE
	integer :: pOrder , NGP, NodG, iFEMtype , iFtype, iMaterial, iDamageParType, iSimplex, iBubble
	integer :: iShiftU, iShiftRes
	type(ptGaussClass) :: PtG
	
	iShiftU = nint(CommonPar(1))
	iShiftRes = nint(CommonPar(2))
	iFEMtype =  nint(commonPar(3)) 
	iFtype = nint(commonPar(4)) 
	iMaterial = nint(commonPar(5))
	iDamageParType = nint(commonPar(6))

	if(iDamageParType>0) call getDamagePar(damagePar,iDamageParType)
	call getMaterial(constLaw, matPar, iMaterial)
	call setFEMtype(iFEMtype,NodG,pOrder,NGP,iSimplex,iBubble)

	call getSliceAllocate(SolU,Sol1,1,NodG,iShiftU + 1 ,iShiftU + NdimE, iDofT)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	
	call PtG%init(Xel,NodG,NdimE,NGP,pOrder,iBubble, iSimplex)
	
	Do nG = 1, NGP ! LoopGauss
		
		call PtG%calcGradU(GradU,SolU,nG)
		
		call getF(F,iFtype) !!! with identity summed
		
		F = F + GradU
		
		call calcPpk(Ppk,F,NdimE,matPar,constLaw)		
		call calcD(D,F,NdimE,matPar,constLaw)	
		
		if(iDamageParType>0) then
			call damageUpdateEquilibrium(F,matPar,damagePar,Param,constLaw,nG) !~ 		
			call damageModifyTangent(D,Ppk,damagePar,Param, nG, 1)
		end if						
		
						
		dV=ptG%dV(nG)
		
		Do A=1, NodG !LoopRow
			ApRow  = (A-1)*iDofT + iShiftRes
			PhiA = PtG%Phi(A,nG)
			dPhiA= PtG%dPhi_G(:,A,nG)
			
			Do p = 1 , NdimE
				App = ApRow+p
				
				BE(App) = BE(App) - dot_product(Ppk(p,:),dPhiA)*DV
				
				AE(App,App) = 1.0d0
								  
			end do ! loop App
		Enddo !LoopRow
	EndDo !LoopGauss
	
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine mechResidualS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
    use ptsGaussLib2
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
    
    Integer :: A, p, ApRow, iShiftRes, Femtype, NodG
    integer , parameter :: NdimE = 2

    Coupling = 0
   
    iShiftRes = nint(commonPar(2))
    femType = nint(commonPar(3)) 
    
    call setNodG(femtype, NodG)
    
    do A = 1,NodG
        ApRow = (A-1) * iDofT + iShiftRes
		
		do p=1,nDimE
			Coupling (ApRow + p, ApRow + p) = 1
		enddo
	enddo

 	
	call numPrint(Coupling)

end Subroutine
