    !> computes the increment of the damage for identification problem
    !!
	!! @param	iShiftU = nint(CommonPar(1))
	!! @param	iShiftL = nint(CommonPar(2))
	!! @param	iShiftD = nint(CommonPar(3))
	!! @param	iShiftDsens = nint(CommonPar(4))
	!! @param	iFEMtype =  nint(commonPar(5)) 
	!! @param	iMaterial = nint(commonPar(6))
	!! @param	alpha = commonPar(7)
    !! @author Rocha, Felipe Figueredo
    

Subroutine computeSens(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
	
	use funcAux
	use globalVariables, only : getF, NdimE, getDamagePar, getMaterial, maxDamagePar, maxMatPar
	use finiteStrainLib
	use ptsGaussLib2
	
    implicit none

    Real*8 :: dummy
    integer :: lengthJParam,lengthParam

    COMMON dummy,LengthParam,LengthJParam
    
    !   ===== SUBROUTINE ARGUMENTS  =======
    integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
    Real*8 :: DelT, DTm,Time ! all reals
    Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

!~     !   =====   END ARGUMENTS  =======
	Integer :: nG, i,j,k,l,m,n,p,q, A, B, ApRow, BpCol , App,Bqp, ipDim, ipDimE, ipDof,ipDofmod ! counters 
	Real*8 :: PhiA,PhiB, dV
	Real*8 :: dPhiA(NdimE), dPhiB(NdimE)
    Real*8 :: MatPar(maxMatPar), deltaDamage, damage , damageNew, alpha
    integer :: damageEvol,constLaw, Idamage, IdamageNew, iShiftU, iShiftL, iFEMtype, iMaterial, iShiftD, iShiftDsens
    integer :: pOrder, iBu, iSimplex, NGP, NodG
    real*8 , allocatable ::  SolU(:), SolL(:) , Xel(:) !! all have dimension NodG*NdimE
    real*8 :: PK(NdimE,NdimE) , GradU(NdimE,NdimE) , F(NdimE,NdimE) , GradL(NdimE,NdimE)
	type(ptGaussClass) :: PtG
    
	iShiftU = nint(CommonPar(1))
	iShiftL = nint(CommonPar(2))
	iShiftD = nint(CommonPar(3))
	iShiftDsens = nint(CommonPar(4))
	iFEMtype =  nint(commonPar(5)) 
	iMaterial = nint(commonPar(6))
	alpha = commonPar(7)
	
	call getMaterial(constLaw, matPar, iMaterial)
	call setFEMtype(iFEMtype,NodG,pOrder,NGP,iSimplex,iBu)

	call getSliceAllocate(SolU,Sol1,1,NodG,iShiftU + 1 ,iShiftU + NdimE, iDofT)
	call getSliceAllocate(SolL,Sol1,1,NodG,iShiftL + 1 ,iShiftL + NdimE, iDofT)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	
	call PtG%init(Xel,NodG,NdimE,NGP,pOrder,iBu, iSimplex)
	
	deltaDamage = 0.0d0
	
	Do nG = 1, NGP ! LoopGauss
		
		!!! SolU
		call PtG%calcGradU(GradU,SolU,nG)

		F = deltaKron(1:NdimE,1:NdimE) + GradU
		
		call calcPpk(PK,F,NdimE,matPar,constLaw)	
		
		!!! SolL
		call PtG%calcGradU(GradL,SolL,nG)

		dV=ptG%dV(nG)
			
		deltaDamage = deltaDamage + dot_product2(PK,GradL)*dV					
	EndDo !LoopGauss
	
	Do nG = 1, NGP ! LoopGauss
		
		!!! SolU
		call PtG%calcGradU(GradU,SolU,nG)

		F = deltaKron(1:NdimE,1:NdimE) + GradU
		
		call calcPpk(PK,F,NdimE,matPar,constLaw)	
		
		!!! SolL
		call PtG%calcGradU(GradL,SolL,nG)

		dV=ptG%dV(nG)
		
		Do A=1, NodG !LoopRow
			ApRow  = (A-1)*iDofT + iShiftDsens + 1
			PhiA = PtG%Phi(A,nG)
			dPhiA= PtG%dPhi_G(:,A,nG)
			ipDim = (A-1)*NdimE

			AE(ApRow,ApRow) = 1.0d0
			BE(ApRow) = BE(ApRow) - PhiA*dot_product2(PK,GradL)*dV
	
			do B= 1 , NodG ! LoopCol ! not considering the simmetry
				BpCol  = (B-1)*iDofT + iShiftD + 1
				PhiB  = PtG%Phi(B,nG)
				dPhiB = PtG%dPhi_G(:,B,nG) 
				
				BE(ApRow) = BE(ApRow) + 2.0d0*alpha*Sol1(BpCol)*dot_product(dPhiA,dPhiB)*dV
				
			Enddo !LoopCol
		Enddo !LoopRow
	EndDo !LoopGauss

	deallocate(SolU,Xel,SolL)
	
end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine computeSensS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use ptsGaussLib2
	IMPLICIT NONE

    Integer :: MaxLRows,Ndim,iDofT, iAdd, ipRow
    Integer ::  Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
    
    Integer :: A, ApRow, iShiftDsens , iFemtype, NodG
	
	iShiftDsens = nint(CommonPar(4))
	iFEMtype =  nint(commonPar(5)) 
	    
    call setNodG(ifemtype, NodG)
	
	Do A=1, NodG !LoopRow
		ApRow  = (A-1)*iDofT + iShiftDsens + 1
		Coupling(ApRow,ApRow) = 1
	Enddo !LoopRow

    return
end Subroutine


