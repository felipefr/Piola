    !> computes the increment of the damage for identification problem
    !!
	!! @param	iShiftU = nint(CommonPar(1))
	!! @param	iShiftRes = nint(CommonPar(2))
	!! @param	iFEMtype =  nint(commonPar(3)) 
	!! @param	iMaterial = nint(commonPar(4))
    !! @author Rocha, Felipe Figueredo
    

Subroutine fillParamSensibility_mechRes(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, &
										Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
	
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
	Integer :: nG, i,ip, p, A  ! counters 
	Real*8 :: dV , dPhiA(NdimE)
    Real*8 :: MatPar(maxMatPar)
    integer :: constLaw, Idamage, IdamageNew, iShiftU, iShiftRes, iFEMtype, iMaterial
    integer :: pOrder, iBu, iSimplex, NGP, NodG
    real*8 , allocatable ::  SolU(:), SolRes(:), Xel(:) !! all have dimension NodG*NdimE
    real*8 :: PK(NdimE,NdimE) , GradU(NdimE,NdimE) , F(NdimE,NdimE) 
	type(ptGaussClass) :: PtG
	
	real*8 :: derCost, cost
	integer :: Icost, IderCost, IcostNew, IderCostNew
    
	iShiftU = nint(CommonPar(1))
	iShiftRes = nint(CommonPar(2))
	iFEMtype =  nint(commonPar(3)) 
	iMaterial = nint(commonPar(4))
	
	call getMaterial(constLaw, matPar, iMaterial)
	call setFEMtype(iFEMtype,NodG,pOrder,NGP,iSimplex,iBu)

	call getSliceAllocate(SolU,Sol1,1,NodG,iShiftU + 1 ,iShiftU + NdimE, iDofT)
	call getSliceAllocate(SolRes,Sol1,1,NodG,iShiftRes + 1 ,iShiftRes + NdimE, iDofT)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	
	call PtG%init(Xel,NodG,NdimE,NGP,pOrder,iBu, iSimplex)
		
	derCost = 0.0d0
	Cost = 0.0d0
	Do nG = 1, NGP ! LoopGauss
		
		dV=ptG%dV(nG)
		
		!!! SolU
		call PtG%calcGradU(GradU,SolU,nG)
		F = deltaKron(1:NdimE,1:NdimE) + GradU
		call calcPpk(PK,F,NdimE,matPar,constLaw)	
	
		Do A=1, NodG !LoopRow
			ip  = (A-1)*NdimE 
			dPhiA= PtG%dPhi_G(:,A,nG)
			
			Do p = 1 , NdimE
				derCost = derCost - dot_product(PK(p,:),dPhiA)*DV*SolRes(ip + p)			  
			end do ! loop App
		Enddo !LoopRow
		
	end do 
	
	!! Estimative for the cost inside the element (the real cost is given globally by the norm of residual vector evaluated in each node)
	
	do i = 1 , NodG
		ip = (i-1)*NdimE
		Cost = Cost + 0.5d0*( SolRes(ip+1)**2.0d0 + SolRes(ip+2)**2.0d0 ) 
	end do 
	 
	Cost = Cost/3.0d0 
	
	ICost = 3
	IderCost = 2
	ICostNew = lengthParam + ICost
	IderCostNew = lengthParam + IderCost

	Param(ICostNew) = Cost
    Param(IderCostNew) = derCost		
	
	AE(1,1) = 1.0d0
	BE(1) = Sol1(1)

	deallocate(SolU,SolRes,Xel)
	
end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine fillParamSensibility_mechResS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    IMPLICIT NONE

    Integer :: MaxLRows,Ndim,iDofT, iAdd, ipRow
    Integer ::  Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
    
	Coupling(1,1) = 1
   
end Subroutine


