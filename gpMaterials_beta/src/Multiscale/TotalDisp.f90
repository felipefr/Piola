    !> Implements periodic boundary condition for 2D continuum problems 
    !!
	!! @param iShiftFluc = nint(CommonPar(1))
	!! @param iShiftU = nint(CommonPar(2))
	!! @param iShiftDeltaU  = nint(CommonPar(3))
	!! @param iFemType  = nint(commonPar(4))
	!! @param iFtype  = nint(CommonPar(3))
	!! @param y0  = commonPar(6:7)
	!!
    !! @author Rocha, Felipe Figueredo
Subroutine TotalDisp &
    (AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)

	use globalVariables , only : NdimE, getF
	use ptsGaussLib2
    
    implicit none
    
    !   ===== SUBROUTINE ARGUMENTS  =======
!~     integer , save :: ElmNumber = 0
    integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
    Real*8 :: DelT, DTm,Time ! all reals
    Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

    !   =====   END ARGUMENTS  =======
	
    integer ::  i ,j, iShiftU, iShiftDeltaU, iShiftFluc ,ipDim, &
				ipU, ipDeltaU, ipp , ipf, iFType, iFemtype, NodG
    Real*8 :: F(NdimE,NdimE) , y0(NdimE), y(NdimE) 
    Real*8 , allocatable :: Xel(:)
	
	iShiftFluc = nint(commonPar(1))
	iShiftU = nint(commonPar(2))
	iShiftDeltaU = nint(commonPar(3))	
	iFemtype = nint(commonPar(4))	
	iFtype = nint(commonPar(5))
	y0(:) = CommonPar(6:7)
	y0 = 0.0d0
	call setNodG(iFemtype, NodG)
	
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	
	call getF(F,iFtype) !!! with identity summed
	
!~ 	F(1:NdimE,1:NdimE) = F(1:NdimE,1:NdimE) - deltaKron(1:NdimE,1:NdimE) !!! really we need F-I
	F(1,1) = F(1,1) - 1.0d0
	F(2,2) = F(2,2) - 1.0d0
	 
    do i=1,NodG
		
		ipDim = (i - 1)* NdimE !! for dimension  !!! paradoxally, taking i := i+NodG , becomes better
        ipf = (i-1) * iDofT + iShiftFluc !! for fluctuation 
        ipU = (i-1) * iDofT + iShiftU !! for total micro displacement = F (y-yo) + fluctuation 
		ipDeltaU = (i-1) * iDofT + iShiftDeltaU !! total increment   
       
        y(:) = Xel(ipDim + 1 : ipDim+NdimE)
        
        y = matmul(F,y-y0)
        
        do j=1,NdimE            
            AE (ipU + j , ipU + j)  = 1.0d0
            BE (ipU + j ) = Sol1 (ipf + j) + y(j)
            if(iShiftDeltaU >= 0 .and. (iShiftDeltaU /= iShiftFluc)) then
				AE (ipDeltaU + j , ipDeltaU + j)  = 1.0d0
!~ 				BE (ipDeltaU + j ) = Sol1 (ipf + j) - Sol0(ipf + j)
				BE (ipDeltaU + j ) = 0.0d0
		    end if
        enddo
         
    enddo
!~ 	   
!~ 	   AE(1,1) = 1.0d0
!~ 	   BE(1) = Sol1(1)
end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine TotalDispS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
    use ptsGaussLib2
    
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)

    integer :: i, j, ipU , ipDeltaU, iShiftFluc, iShiftU, iShiftDeltaU, NodG, iFemType
    integer , parameter :: NdimE = 2
	
!~ 	Coupling = 0

	iShiftFluc = nint(commonPar(1)) 
	iShiftU = nint(commonPar(2)) 
	iShiftDeltaU = nint(commonPar(3)) 
	iFemtype = nint(commonPar(4)) 
	
	call setNodG(iFemtype, NodG)
	
    do i=1,NodG
        ipU = (i-1) * iDofT + iShiftU !! for total micro displacement = F (y-yo) + fluctuation 
        ipDeltaU = (i-1) * iDofT + iShiftDeltaU !! total increment 
        
        do j=1,NdimE
            Coupling (ipU + j , ipU + j)  = 1	
            if(iShiftDeltaU >= 0 .and. (iShiftDeltaU /= iShiftFluc)) then
				Coupling (ipDeltaU + j , ipDeltaU + j)  = 1	
			end if
        enddo
    enddo
!~ 	
!~ 	write(*,*) "computeTotalDisp2DS Symbolic"
!~ 	call IprintMat(Coupling)
!~ 	

end Subroutine

