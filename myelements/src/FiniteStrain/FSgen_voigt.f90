    !> Generic element for finite strains
    !!
    !!
	!! @param iShiftU = nint(CommonPar(1))
	!! @param iShiftDeltaU = nint(CommonPar(2))
	!! @param iFemType  = nint(CommonPar(3)) 
	!! @param iFType  = nint(commonPar(4))
	!! @param iMaterial  = nint(commonPar(5))
	!!
    !! @author Rocha, Felipe Figueredo

Subroutine FSgen_voigt(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
	use globalVariables, only : NdimE, getMaterial, maxMatPar
	use finiteStrainLib
	use ptsGaussLib2
	
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT, NodElt! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

!~ 	!   =====   END ARGUMENTS  =======
         
	Integer :: nG, i,j,k,l,m,n,p,q, A, B, ApRow, BpCol , App,Bqp, ipDim, ipDimE, ipDof,ipDofmod ! counters 
	Integer :: constLaw 
	Real*8 :: dV , Det , MatPar(maxMatPar)
	Real*8 :: GradU(NdimE,NdimE) , F(NdimE,NdimE)
	Real*8 :: D(NdimE,NdimE,NdimE,NdimE), Ppk(NdimE,NdimE), Dmat(NdimE*NdimE,NdimE*NdimE) , Pmat(NdimE*NdimE)
	Real*8 :: Bmat(NdimE*NdimE,NodElt*NdimE), BmatT(NodElt*NdimE,NdimE*NdimE), Auu(NodElt*NdimE,NodElt*NdimE), Bu(NodElt*NdimE)
	real*8 , allocatable ::  SolU(:) , Xel(:) !! all have dimension NodG*NdimE
	integer :: pOrder , NGP, NodG, iFEMtype , iFtype, iMaterial, iSimplex, iBubble
	integer :: iShiftU, iShiftDeltaU
	type(ptGaussClass) :: PtG
	
	iShiftU = nint(CommonPar(1))
	iShiftDeltaU = nint(CommonPar(2))
	iFEMtype =  nint(commonPar(3)) 
	iFtype = nint(commonPar(4)) 
	iMaterial = nint(commonPar(5))

	call getMaterial(constLaw, matPar, iMaterial)
	call setFEMtype(iFEMtype,NodG,pOrder,NGP,iSimplex,iBubble)

	call getSliceAllocate(SolU,Sol1,1,NodG,iShiftU + 1 ,iShiftU + NdimE, iDofT)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	
	call PtG%init(Xel,NodG,NdimE,NGP,pOrder,iBubble, iSimplex)
	
	AE =  0.0d0
	BE =  0.0d0
	
	Auu = 0.0d0
	Bu = 0.0d0
	
	Do nG = 1, NGP ! LoopGauss
		
		call PtG%calcGradU(GradU,SolU,nG)
		
		F = deltaKron(1:NdimE,1:NdimE) + GradU
		
		call calcPpk(Ppk,F,NdimE,matPar,constLaw)		
		call calcD(D,F,NdimE,matPar,constLaw)						
		
		call PtG%getBmat(Bmat,nG)
		BmatT = transpose(Bmat)				
		call voigtTen4toTen2(Dmat,D)
		call voigtTen2toVec(Pmat,Ppk)
		
		dV=ptG%dV(nG)
		
		Bu = Bu - matmul(BmatT,Pmat)*dV
		Auu = Auu + matmul(BmatT,matmul(Dmat,Bmat))*dV
		
	EndDo !LoopGauss
	
	if(iShiftU == iShiftDeltaU ) Bu = Bu + matmul(Auu,SolU)
	
	call setBlockToVector(BE,Bu,NodG,iDofT, NdimE, iShiftDeltaU ) !(A,B, NodG, iDofA, iDofB, iShiftA) , B in A
	call setBlockToMatrixSquare(AE,Auu, NodG,iDofT, NdimE, iShiftDeltaU ) !(A,B, NodG, iDofA, iDofB, iShiftA) , B in A
	
 	if(iShiftU /= iShiftDeltaU) then
		call setBlockDiagonalConstantToMatrix(AE,1.0d0,NodG, iDofT, NdimE, iShiftU, iShiftU)
		call setBlockDiagonalConstantToMatrix(AE,-1.0d0,NodG, iDofT, NdimE, iShiftU, iShiftDeltaU)
		call setBlockToVector(BE,SolU,NodG,iDofT, NdimE, iShiftU ) !(A,B, NodG, iDofA, iDofB, iShiftA) , B in A
	end if

end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine FSgen_voigtS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
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
