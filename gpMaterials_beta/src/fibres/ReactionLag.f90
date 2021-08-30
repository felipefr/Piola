Subroutine ReactionLag &
    (AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)

	use funcAux
	use globalVariables, only : getF, NdimE
	use ptsGaussLib2, only : setNodG
    implicit none
    
    !   ===== SUBROUTINE ARGUMENTS  =======
    integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
    Real*8 :: DelT, DTm,Time ! all reals
    Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

    !   =====   END ARGUMENTS  =======
   
	integer :: ApRowU, ApRowLag, A, AppU, AppLag, ApDim, p
    Real*8 :: Ud(NdimE) , pen, LoadPar(6), eps, F(NdimE,NdimE)
    real*8 , allocatable :: Xel(:) !! all have dimension NodG*NdimE
    integer :: iShiftU, iShiftLag, iFtype , NodG, iFemType
	
	iShiftU = nint(commonPar(1))
	iShiftLag = nint(commonPar(2))
	iFemType = nint(commonPar(3))
	iFtype = nint(commonPar(4))
	pen = commonPar(5)
	eps = commonPar(6)
	
	call setNodG(iFemtype, NodG)	
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	call getF(F,iFtype)
		
	Ud = 0.0d0
	
	Do A=1, NodG !LoopRow
		ApRowU  = (A-1)*iDofT + iShiftU
		ApRowLag = (A-1)*iDofT + iShiftLag
		ApDim = (A-1)*NdimE
		
		Ud = matmul(F,Xel(ApDim + 1: ApDim + NdimE)) - Xel(ApDim + 1: ApDim + NdimE) 
		
		Do p = 1 , NdimE
			AppU = ApRowU+p
			AppLag = ApRowLag+p
			
			BE(AppLag) = pen*Ud(p)
!~ 			BE(AppU) = pen*Ud(p)
!~ 			AE(AppU,AppU) = eps
			AE(AppU,AppLag) = pen
			AE(AppLag, AppU) = pen
			AE(AppLag, AppLag) = eps
			
		end do
	
	end do
		
end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine ReactionLagS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	use funcAux
	use ptsGaussLib2, only : setNodG
    
    implicit none
    
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
    
    integer :: ApRowU, ApRowLag, A, AppU, AppLag , p
    integer :: iShiftU, iShiftLag, iFtype , NodG, NdimE, iFemType
    
    iShiftU = nint(commonPar(1))
	iShiftLag = nint(commonPar(2))
	iFemType = nint(commonPar(3))
	
	call setNodG(iFemtype, NodG)
	NdimE = 2

	Do A=1, NodG !LoopRow
		ApRowU  = (A-1)*iDofT + iShiftU
		ApRowLag = (A-1)*iDofT + iShiftLag
		
		Do p = 1 , NdimE
			AppU = ApRowU+p
			AppLag = ApRowLag+p
			
!~ 			Coupling(AppU,AppU) = 1
			Coupling(AppU,AppLag) = 1
			Coupling(AppLag, AppU) = 1
			Coupling(AppLag, AppLag) = 1
			
		end do
	
	end do
	
	call numprint(Coupling)
	
end Subroutine
