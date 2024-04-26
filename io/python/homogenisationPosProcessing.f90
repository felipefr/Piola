module converterFuncs

use utils
implicit none

contains

subroutine homogenisationFibres_(PK, ParamFile, Nelem, NTimeSteps, NParamCell, NdimE, iShift_volfib , iShift_afib, iShift_sfib )    

real*8 , intent(out) :: PK(:,:,:) 
integer, intent(in) :: Nelem, NTimeSteps, NParamCell, NdimE, iShift_volfib , iShift_afib, iShift_sfib 
character(len=40), intent(in) :: ParamFile

integer :: i, j,n,k, Iafib, Isfib, Ivolfib, IUnitParam, OUnitHomog
			
Real*8 , allocatable :: Param(:), sfib(:), afib(:)
Real*8 :: dummy, VolTot , Volfib

IUnitParam = 12
OUnitHomog = 22

open (IUnitParam, file=trim(ParamFile), status= 'old' )

allocate(Param(NElem*NParamCell), sfib(NdimE) , afib(NdimE))

do n = 1 , NTimeSteps 
	call FindKeyword(IUnitParam,"*Real Parameters")
	read(IUnitParam,*) ( dummy , j=1,NElem + 1)
	read(iUnitParam,*) Param

	volTot = 0.0d0
	do k = 1 , Nelem		
		Ivolfib = (k-1)*NParamCell + iShift_volfib + 1  
		volfib = Param(Ivolfib)
		volTot = volTot + volfib
			
		Iafib =  (k-1)*NParamCell + iShift_afib + 1
		afib = Param(Iafib:Iafib+NdimE-1)
		
		Isfib =  (k-1)*NParamCell + iShift_sfib + 1
		sfib = Param(Isfib:Isfib+NdimE-1)
		
		do j = 1 , NdimE
		do i = 1 , NdimE
			PK(i,j,n) =  PK(i,j,n) + volfib*sfib(i)*afib(j)
		end do 
		end do 
				
	end do
	
	PK(:,:,n) = (1.0d0/VolTot) * PK(:,:,n)
	
	
end do

close(IUnitParam)

end subroutine 


end module
