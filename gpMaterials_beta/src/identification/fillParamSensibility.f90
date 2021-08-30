    !> computes the increment of the damage for identification problem
    !!
	!! @param	iShiftU = nint(CommonPar(1))
	!! @param	iShiftUd = nint(CommonPar(2))
	!! @param	iShiftL = nint(CommonPar(3))
	!! @param	iFEMtype =  nint(commonPar(4)) 
	!! @param	iMaterial = nint(commonPar(5))
	!! @param	weight1 = CommonPar(6)
	!! @param	weight2 = CommonPar(7)
	!! @param	normL2_ref = CommonPar(8)
	!! @param	normEnergy_ref = CommonPar(9)
    !! @author Rocha, Felipe Figueredo
    

Subroutine fillParamSensibility(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
	
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
    integer :: damageEvol,constLaw, Idamage, IdamageNew, iShiftU, iShiftUd, iShiftL, iFEMtype, iMaterial, iShiftD, iShiftDsens
    integer :: pOrder, iBu, iSimplex, NGP, NodG
    real*8 , allocatable ::  SolU(:), SolUd(:), SolL(:) , Xel(:) !! all have dimension NodG*NdimE
    real*8 :: PK(NdimE,NdimE) , GradU(NdimE,NdimE) , F(NdimE,NdimE) , GradL(NdimE,NdimE), SolDif (6), Udif(NdimE), U(NdimE)
	type(ptGaussClass) :: PtG
	
	real*8 :: derCost, cost, weight1, weight2, normL2_ref, normEnergy_ref, energy, weight3
	integer :: Icost, IderCost, IcostNew, IderCostNew
    
	iShiftU = nint(CommonPar(1))
	iShiftUd = nint(CommonPar(2))
	iShiftL = nint(CommonPar(3))
	iFEMtype =  nint(commonPar(4)) 
	iMaterial = nint(commonPar(5))
	weight1 = CommonPar(6)
	weight2 = CommonPar(7)
	normL2_ref = CommonPar(8)
	normEnergy_ref = CommonPar(9)
	weight3 = CommonPar(10)
	
	weight1 = weight1/normL2_ref
	weight2 = weight2/normEnergy_ref
	
	call getMaterial(constLaw, matPar, iMaterial)
	call setFEMtype(iFEMtype,NodG,pOrder,NGP,iSimplex,iBu)

	call getSliceAllocate(SolU,Sol1,1,NodG,iShiftU + 1 ,iShiftU + NdimE, iDofT)
	call getSliceAllocate(SolUd,Sol1,1,NodG,iShiftUd + 1 ,iShiftUd + NdimE, iDofT)
	call getSliceAllocate(SolL,Sol1,1,NodG,iShiftL + 1 ,iShiftL + NdimE, iDofT)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	
	call PtG%init(Xel,NodG,NdimE,NGP,pOrder,iBu, iSimplex)
	
	SolDif = SolUd - SolU
	
	derCost = 0.0d0
	Cost = 0.0d0
	Do nG = 1, NGP ! LoopGauss
		
		dV=ptG%dV(nG)
		
		!!! SolU
		call PtG%calcGradU(GradU,SolU,nG)
		call PtG%calcU(U,SolU,nG)
		F = deltaKron(1:NdimE,1:NdimE) + GradU
		call calcPpk(PK,F,NdimE,matPar,constLaw)	
		
		!!! SolL
		call PtG%calcGradU(GradL,SolL,nG)
		
		!!! SolDifU
		GradU = 0.0d0
		call PtG%calcGradU(GradU,SolDif,nG)
		F = deltaKron(1:NdimE,1:NdimE) + GradU
		call strainEnergy(energy,F,matPar,constLaw)
		
		call PtG%calcU(Udif,SolDif,nG)
		
	
		derCost = derCost + dot_product2(PK,GradL)*dV					
		
		Cost = Cost +  weight1 * dot_product(Udif,Udif)*DV 
		Cost = Cost + weight2 * energy * DV
		Cost = Cost + weight3 * dot_product(U,U)*DV 

	end do 
	
	ICost = 3
	IderCost = 2
	ICostNew = lengthParam + ICost
	IderCostNew = lengthParam + IderCost

	Param(ICostNew) = Cost
    Param(IderCostNew) = derCost		
	
	AE(1,1) = 1.0d0
	BE(1) = Sol1(1)

	deallocate(SolU,SolUd,Xel,SolL)
	
end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine fillParamSensibilityS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    IMPLICIT NONE

    Integer :: MaxLRows,Ndim,iDofT, iAdd, ipRow
    Integer ::  Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
    
	Coupling(1,1) = 1
   
end Subroutine


