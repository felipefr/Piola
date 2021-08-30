!Last Modification : 11/03/2014
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Element for Hyperelasticity 
! Still in development
!     ------------------------------------------------------------------

!~ 
Subroutine postProcessingRVEboundary(AE, BE, MaxLRows, XLL, NDim, iDofT, &
 NodElt, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux

	implicit none

	!	PARAMETERS

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT,NodElt ! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

	!   =====   END ARGUMENTS  =======
         
	integer, save :: iElem = 0
    integer , parameter :: OUnitPKhom = 48, nodG = 2, NdimE = 2
	integer :: iShiftPKhom, sizePKhom, Nelem, i, j, ip
	real*8 :: PKhom(NdimE,NdimE)

	iShiftPKhom = nint(commonPar(1))
	Nelem = nint(commonPar(2))
	
	sizePKhom = NdimE*NdimE	
	
	do i = 1 , NdimE
		ip = NodG*idofT + iShiftPKhom + (i-1)*NdimE
		do j = 1 , NdimE
			PKhom(i,j) = -Sol1(ip + j) 
		end do
	end do
	
	AE(1,1) = 1.0d0
	BE(1) = Sol1(1)
	
	open (OUnitPKhom, file='PKhom_Gamma.txt', Access = 'append')  
	
!~ 	!!! ======== For RVE =========
	if(mod(iElem,Nelem) == 0) then 
		write(OUnitPKhom,*) PKhom(1,1)
		write(OUnitPKhom,*) PKhom(1,2)
		write(OUnitPKhom,*) PKhom(2,1)
		write(OUnitPKhom,*) PKhom(2,2)
	end if
	iElem = iElem + 1
	
	close (OUnitPKhom)  

End Subroutine



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine postProcessingRVEboundaryS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd, ipRow , i
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)

	Coupling=0 
	
	Coupling(1,1) = 1
	
	return
end Subroutine
