Subroutine insertDeformation &
    (AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)

	use globalVariables
	use funcAux
    
    implicit none
    
    !   ===== SUBROUTINE ARGUMENTS  =======
!~     integer , save :: ElmNumber = 0
    integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
    Real*8 :: DelT, DTm,Time ! all reals
    Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

    !   =====   END ARGUMENTS  =======
	
	integer, parameter :: nLoadParam = 9, iShiftLoadPar = 3
    integer ::  i ,j, ip, LoadType , LoadTypeProg , nbFtype, iShiftMinDetQ, evolStyle
    Real*8 :: LoadPar(6) 
    logical :: isIncremental , stopAfterBif
	
	iShiftMinDetQ = nint(commonPar(1))
	evolStyle = nint(commonPar(2))
	nbFtype = nint(commonPar(3))

	if(.not. isAllocated_Fsave() ) then
		call allocateFsave(NdimE,nbFtype)
	end if
	
	do i = 1 , nbFtype
		ip = (i-1)*nLoadParam + iShiftLoadPar
		
		stopAfterBif = .false.
		if(nint(commonPar(ip+1))>0) stopAfterBif = .true.
		
		LoadType = nint(commonPar(ip + 2))
		LoadTypeProg = nint(commonPar(ip + 3))
		LoadPar(:) = CommonPar( ip + 4 : ip + 9)
		
		if(.not. (stopAfterBif .and. Sol1( iShiftMinDetQ + 1) < 0.0d0 ) ) then
			call evolutionF(i,LoadPar,Time,delT,LoadTypeProg,LoadType,NdimE,evolStyle)
		end if
	end do

!~ 	call numprint(Fsave(:,:,1))
!~ 	call numprint(Fsave(:,:,2))
!~ 	pause
 	
	AE(1,1) = 1.0d0
	BE(1) = Sol1(1)

end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine insertDeformationS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)

	Coupling(1,1) = 1

end Subroutine
