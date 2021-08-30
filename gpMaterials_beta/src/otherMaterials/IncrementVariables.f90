	!> Element that solves the trivial diagonal problem of incrementing variables
	!! @param CommonPar = ( a , a*[bi,ci,di,ei] , dummy ) where : 
	!! @param Nvar = CommonPar(1) ==> a
	!! @param iShiftUi = CommonPar(see structure) ===> bi
	!! @param iShiftDeltaUi = CommonPar(see structure) ===> ci
	!! @param nodei (for each material = see structure) ===> ei 
	!! @author Rocha, Felipe Figueredo

Subroutine IncrementVariables(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------

	use funcAux
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

	!   =====   END ARGUMENTS  =======
    
    integer , parameter :: NauxParam = 3  	
    integer :: Nvar, iShiftDeltaU, iShiftU, lenU, node, NodG , node1, node2
    integer :: i, ip, jp, jpDelta, j, k
         	
	Nvar = nint(commonPar(1))
	
	NodG = 3
	
	do i = 1 , Nvar
		ip = (i-1)*NauxParam + 1
		iShiftU = nint(commonPar(ip + 1))
		iShiftDeltaU = nint(commonPar(ip + 2))
		node = nint(commonPar(ip + 3))
		
		lenU = iShiftU - iShiftDeltaU
				
		if(node == 0) then
			node1 = 1
			node2 = NodG
		else
			node1 = node
			node2 = node		
		end if
		
		do j = node1,node2
			jp = (j - 1)*iDofT + iShiftU
			jpDelta = (j - 1)*iDofT + iShiftDeltaU
			
			do k = 1,lenU
				AE(jp + k, jp + k) = 1.0d0
				BE(jp + k) = Sol1(jpDelta + k) + Sol0(jp + k)
			end do	
		end do
		
	end do
	
End

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Subroutine IncrementVariablesS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	use funcAux
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd, ipRow
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)

    integer , parameter :: NauxParam = 3  	
    integer :: Nvar, iShiftDeltaU, iShiftU, lenU, node, NodG , node1, node2
    integer :: i, ip, jp, jpDelta, j, k
         	
	Nvar = nint(commonPar(1))
	
	NodG = 3
	
	do i = 1 , Nvar
		ip = (i-1)*NauxParam + 1
		iShiftU = nint(commonPar(ip + 1))
		iShiftDeltaU = nint(commonPar(ip + 2))
		node = nint(commonPar(ip + 3))
		
		lenU = iShiftU - iShiftDeltaU
				
		if(node == 0) then
			node1 = 1
			node2 = NodG
		else
			node1 = node
			node2 = node		
		end if
		
		do j = node1,node2
			jp = (j - 1)*iDofT + iShiftU
			jpDelta = (j - 1)*iDofT + iShiftDeltaU
			
			do k = 1,lenU
				Coupling(jp + k, jp + k) = 1
			end do	
		end do
		
	end do

!~ 	call numprint(Coupling)
!~ 	pause
end Subroutine

