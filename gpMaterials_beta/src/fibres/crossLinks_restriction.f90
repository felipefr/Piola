    !> Incorporates restriction in cross-links for mikado network : \n
    !! 4   2 \n
	!!  \ /  \n
    !!   5   \n
	!!  / \  \n
    !! 1   3 \n
	!! @param iShiftFluc = nint(CommonPar(1))
	!! @param iShiftLag = nint(CommonPar(2))
	!! @param pen  = commonPar(3)
	!! @param eps  = commonPar(4)

    !! @author Rocha, Felipe Figueredo

!     ------------------------------------------------------------------
Subroutine crossLinks_restriction(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, &
						Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux , only : numprint, setAEandBE_lagrangeMultipliers, getSliceAllocate
	use globalVariables, only : NdimE, getMaterial, maxMatPar, getAnisoTensorInv
	use fibresLib
		
	IMPLICIT NONE
	
	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	! Reals Vectors and matrices
	Real*8  :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), &
		Param(*), JParam(*), CommonPar(*)
	
	!   =====   END ARGUMENTS  =======
     		
	Integer :: i,j,k,l,m,n,p,q, A, B, ApRow, BpCol , App,Bqp, ApDim, ApDim2, AppDim! counters 
	integer, parameter ::  NodG = 5, NodLag = 5 , NdofLag = 2*NdimE
	Real*8 :: MatPar(maxMatPar)
	real*8 , allocatable ::  Xel(:) 
	Real*8 :: eps, pen 
	Real*8 :: deltaX(NdimE), Lfib(NodG-1), MatB(NdofLag,NodG*NdimE), MatC(NdofLag,NdofLag), t1,t2
	integer :: iShiftFluc, iShiftLag 
	
	iShiftFluc = nint(CommonPar(1))
	iShiftLag = nint(CommonPar(2))
	pen = commonPar(3)
	eps = commonPar(4)

	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	
	ApDim2 = (NodLag-1)*NdimE
	do i = 1, NodG-1
		ApDim = (i-1)*NdimE
		deltaX = Xel(ApDim:ApDim+1) - Xel(ApDim2:ApDim2+1)
		Lfib(i) = dsqrt(dot_product(deltaX,deltaX))  
	end do 
	
	t1 = Lfib(1)/(Lfib(1) + Lfib(2))
	t2 = Lfib(3)/(Lfib(3) + Lfib(4))  
	
	write(0,*) t1,t2
	
	MatB = 0.0d0
	
	MatB(1,1) = t1
	MatB(2,2) = t1
	MatB(1,3) = 1.0d0 - t1
	MatB(2,4) = 1.0d0 - t1
	MatB(3,5) = t2
	MatB(4,6) = t2
	MatB(3,7) = 1.0d0 - t2
	MatB(4,8) = 1.0d0 - t2
	MatB(1,9) = -1.0d0
	MatB(2,10) = -1.0d0
	MatB(3,9) = -1.0d0
	MatB(4,10) = -1.0d0
	
	MatC = 0.0d0
	do i = 1,NdofLag
		MatC(i,i) = eps
	end do
	
	call setAEandBE_lagrangeMultipliers(AE,BE,matB,matC,Sol1,eps,pen,NodG,NodLag,idofT,NdofLag,iShiftFluc,iShiftLag,NdimE)

end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine crossLinks_restrictionS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	use funcAux , only : numprint, setCoupling_lagrangeMultipliers
	
	use globalVariables, only : NdimE
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)
	integer ::  p , q, App, ApRow, BpCol, Bqp, A,B, iShiftFluc, iShiftLag 
	integer , parameter :: NodG = 5, NodLag = 5,  NdofLag = 2*NdimE
	
	iShiftFluc = nint(CommonPar(1))
	iShiftLag = nint(CommonPar(2))
 	
	call setCoupling_lagrangeMultipliers(Coupling,NodG,NodLag,idofT,NdofLag,iShiftFluc,iShiftLag,NdimE)
	
end subroutine





