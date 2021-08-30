module globalVariables

implicit none

private constLaw_save, matPar_save
public setMaterial, getMaterial, isAllocatedMaterial, allocateMaterial, maxMatPar


integer , allocatable :: constLaw_save(:)
real*8 , allocatable :: matPar_save(:,:)
integer , parameter :: maxMatPar = 9, NdimE = 2

contains

subroutine getMaterial(constLaw,matPar,iMatType)
	integer , intent(out) :: constLaw
	real*8 , intent(out) :: matPar(:)
	integer , intent(in) :: iMatType
	integer n

	constLaw = constLaw_save(iMatType)
	matPar(:) = matPar_save(:,iMatType)
	
end subroutine

subroutine setMaterial(constLaw,matPar,iMatType)
	integer , intent(in) :: constLaw
	real*8 , intent(in) :: matPar(:)
	integer , intent(in) :: iMatType

	constLaw_save(iMatType) = constLaw
	matPar_save(:,iMatType) = matPar(:)
	
end subroutine

logical function isAllocatedMaterial() result(flag)
	flag = allocated(constLaw_save) 
end function

subroutine allocateMaterial(n)
	integer, intent(in) ::  n
	
	if(allocated(constLaw_save)) then
		write(*,*) "trying to allocate Material already allocated"
		stop
	else
		allocate(constLaw_save(n), matPar_save(maxMatPar,n))
	end if
	
end subroutine

end module
