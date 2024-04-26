subroutine homogenisationFibres()    
use utils

implicit none

character(len=40) :: arg, ParamFile
integer :: i, j,n,k,  Nelem, NTimeSteps, NParamCell, NdimE, iShift_volfib , iShift_afib, iShift_sfib , &
		   Iafib, Isfib, Ivolfib, IUnitParam, OUnitHomog
			
Real*8 , allocatable :: Param(:), sfib(:), afib(:), Phom(:,:)
Real*8 :: dummy, VolTot , Volfib

call getarg(2,ParamFile)
write(0,*) ParamFile

call getarg(3, arg)
write(0,*) arg
read(arg,*) Nelem, NTimeSteps, NParamCell, NdimE, iShift_volfib , iShift_afib, iShift_sfib 


IUnitParam = 12
OUnitHomog = 22

open (IUnitParam, file=trim(ParamFile), status= 'old' )
open (OUnitHomog, file='Homog.txt')

allocate(Param(NElem*NParamCell), sfib(NdimE) , afib(NdimE), Phom(NdimE,NdimE))

do n = 1 , NTimeSteps 
	call FindKeyword(IUnitParam,"*Real Parameters")
	read(IUnitParam,*) ( dummy , j=1,NElem + 1)
	read(iUnitParam,*) Param

	Phom = 0.0d0
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
			Phom(i,j) =  Phom(i,j) + volfib*sfib(i)*afib(j)
		end do 
		end do 
				
	end do
	
	Phom = (1.0d0/VolTot) * Phom
	
	write(OUnitHomog,*) Phom
	
end do

close(IUnitParam)
close(OUnitHomog)

end subroutine 

