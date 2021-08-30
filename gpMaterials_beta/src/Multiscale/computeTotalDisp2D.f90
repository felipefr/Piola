! send to cluster
Subroutine computeTotalDisp2D &
    (AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)

	use loadingLib
	use ptsGaussLib2
    
    implicit none
    
    !   ===== SUBROUTINE ARGUMENTS  =======
!~     integer , save :: ElmNumber = 0
    integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
    Real*8 :: DelT, DTm,Time ! all reals
    Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

    !   =====   END ARGUMENTS  =======
	
	integer, parameter :: NdimE = 2
    integer ::  i ,j, iShiftU, iShiftFluc ,ipDim, &
				ip, ipp , ipf, LoadType , LoadTypeProg , iFemtype, NodG
    Real*8 :: LoadPar(6) , F(NdimE,NdimE) , y0(NdimE), y(NdimE) 
    Real*8 , allocatable :: Xel(:)
	
	iShiftFluc = nint(commonPar(1))
	iShiftU = nint(commonPar(2))	
	iFemtype = nint(commonPar(3))	
	LoadType = nint(commonPar(4))
	LoadTypeProg = nint(commonPar(5))
	LoadPar(:) = CommonPar(6:11)
	y0(:) = CommonPar(12:13)
	
	call setNodG(iFemtype, NodG)
	
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	
	call setF(F,LoadPar,Time,DelT,LoadTypeProg,LoadType) !! really F-I
	 
    do i=1,NodG
		
		ipDim = (i - 1)* NdimE !! for dimension  !!! paradoxally, taking i := i+NodG , becomes better
        ipf = (i-1) * iDofT + iShiftFluc !! for fluctuation 
        ip = (i-1) * iDofT + iShiftU !! for total micro displacement = F (y-yo) + fluctuation 
       
        y(:) = Xel(ipDim + 1 : ipDim+NdimE)
        
        y = matmul(F,y-y0)
        
        do j=1,NdimE            
            AE (ip + j , ip + j)  = 1.0d0
            BE (ip + j ) = Sol1 (ipf + j) + y(j) 
        enddo
         
    enddo
!~ 	   
!~ 	   AE(1,1) = 1.0d0
!~ 	   BE(1) = Sol1(1)
end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine computeTotalDisp2DS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
    use ptsGaussLib2
    
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)

    integer :: i, j, ip , iShiftU, NodG, iFemType
    integer , parameter :: NdimE = 2
	
	Coupling = 0

	
	iShiftU = nint(commonPar(2)) 
	iFemtype = nint(commonPar(3)) 
	
	call setNodG(iFemtype, NodG)
	
    do i=1,NodG
		
        ip = (i-1) * iDofT + iShiftU !! for total micro displacement = F (y-yo) + fluctuation 
        
        do j=1,NdimE
            Coupling (ip + j , ip + j)  = 1	
        enddo
    enddo
!~ 	
!~ 	write(*,*) "computeTotalDisp2DS Symbolic"
!~ 	call IprintMat(Coupling)
!~ 	

end Subroutine
