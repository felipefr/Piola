module fibresModelsLib
use funcAux
implicit none
private  fibreExponential, fibreLinear
!~ public StrainEnergy, StrainEnergyC

contains

subroutine FibreExponential(Psi,I4,MatPar) ! see my dissertation
	implicit none
	real*8 , intent(out) :: Psi
	real*8 , intent(in) :: I4, MatPar(2) ! MatPar is only the material coeficient, a0 was splitted
 	real*8 :: c1,c2
	c1 = MatPar(1) ! k1
	c2 = MatPar(2) ! k2
	
!~ 	write(0,*) "c1 = " , c1 , "c2 = " , c2  
	
	Psi = 0.5d0*(c1/c2)*( Dexp(c2*(I4-1.0d0)**2.0d0)-1.0d0 ) 
	
end subroutine 

subroutine dFibreExponential(dPsi,I4,MatPar) ! see my dissertation
	implicit none
	real*8 , intent(out) :: dPsi
	real*8 , intent(in) :: I4, MatPar(2) ! MatPar is only the material coeficient, a0 was splitted
 	real*8 :: c1,c2
	c1 = MatPar(1) ! k1
	c2 = MatPar(2) ! k2
	
!~ 	write(0,*) "c1 = " , c1 , "c2 = " , c2  
	
	dPsi = c1*(I4-1.0d0)*Dexp(c2*(I4-1.0d0)**2.0d0) 
	
end subroutine

subroutine d2FibreExponential(d2Psi,I4,MatPar) ! see my dissertation
	implicit none
	real*8 , intent(out) :: d2Psi
	real*8 , intent(in) :: I4, MatPar(2) ! MatPar is only the material coeficient, a0 was splitted
 	real*8 :: c1,c2
	c1 = MatPar(1) ! k1
	c2 = MatPar(2) ! k2
	
	d2Psi = c1*Dexp(c2*(I4-1.0d0)**2.0d0)*(1.0d0 + 2.0d0*(I4-1.0d0)**2.0d0) 
	
end subroutine

subroutine FibreLinear(Psi,I4,MatPar) ! see my Wriggers  Sx = E Ex
	implicit none
	real*8 , intent(out) :: Psi
	real*8 , intent(in) :: I4, MatPar ! MatPar is only the material coeficient, a0 was splitted
 	real*8 :: c1
	c1 = MatPar ! k1
	
!~ 	write(0,*) "c1 = " , c1
	
	Psi = (c1/8.0d0)*(I4-1.0d0)**2.0d0 
	
end subroutine 

subroutine dFibreLinear(dPsi,I4,MatPar) ! see my Wriggers  Sx = E Ex
	implicit none
	real*8 , intent(out) :: dPsi
	real*8 , intent(in) :: I4, MatPar ! MatPar is only the material coeficient, a0 was splitted
 	real*8 :: c1
	c1 = MatPar ! k1
	
!~ 	write(0,*) "c1 = " , c1   
	
	dPsi = (c1/4.0d0)*(I4-1.0d0) 
	
end subroutine 

subroutine d2FibreLinear(d2Psi,I4,MatPar) ! see my Wriggers  Sx = E Ex
	implicit none
	real*8 , intent(out) :: d2Psi
	real*8 , intent(in) :: I4, MatPar ! MatPar is only the material coeficient, a0 was splitted
 	real*8 :: c1
	c1 = MatPar ! k1
	
	d2Psi = c1/4.0d0 
	
end subroutine 

subroutine strainEnergy(energy,qfib,matPar,constLaw)
	real*8, intent(in) :: qfib(:) , matPar(:)
	real*8 , intent(out) :: energy
	integer , intent(in) :: constLaw
	real*8 :: I4
	
	energy = 0.0d0
	
	I4 = dot_product(qfib,qfib)
	
	if(I4<1.0d0) then 
		return
	end if
	
	selectcase(constLaw) 
		case(1) 
			call FibreExponential(energy,I4,matPar(1:2))
		case(2) 
			call FibreLinear(energy,I4,matpar(1))
		case default
			write(0,*) constLaw, matpar(1:3), "Strain Energy case unknown"
			PAUSE
	end select

end subroutine

subroutine SfibAnalytic(Sfib,qfib,matPar,constLaw)
	real*8, intent(in) :: qfib(:) , matPar(:)
	real*8 , intent(out) :: Sfib(:)
	integer , intent(in) :: constLaw
	real*8 :: I4, dPsi
	
	Sfib = 0.0d0
	
	I4 = dot_product(qfib,qfib)

	if(I4<1.0d0) then 
		return
	end if
	
	select case(constLaw) 
		case(1) 
			call dFibreExponential(dPsi,I4,matPar(1:2))
		case(2) 
			call dFibreLinear(dPsi,I4,matpar(1))
		case default
			write(0,*) constLaw, matpar(1:3), "Strain Energy case unknown"
			PAUSE
	end select
	
	Sfib = 2.0d0*dPsi*qfib

end subroutine

subroutine DSfibAnalytic(DSfib,qfib,matPar,constLaw)
	real*8, intent(in) :: qfib(:) , matPar(:)
	real*8 , intent(out) :: DSfib(:,:)
	integer , intent(in) :: constLaw
	real*8 :: I4, dPsi, d2Psi
	
	DSfib = 0.0d0
	
	I4 = dot_product(qfib,qfib)
	
	if(I4<1.0d0) then 
		return
	end if
	
	select case(constLaw) 
		case(1) 
			call dFibreExponential(dPsi,I4,matPar(1:2))
			call d2FibreExponential(d2Psi,I4,matPar(1:2))
		case(2) 
			call dFibreExponential(dPsi,I4,matPar(1:2))
			call d2FibreLinear(d2Psi,I4,matpar(1))
		case default
			write(0,*) constLaw, matpar(1:3), "Strain Energy case unknown"
			PAUSE
	end select
	
	DSfib(1,1) = 4.0d0*d2Psi*qfib(1)*qfib(1) + 2.0d0*dPsi
	DSfib(1,2) = 4.0d0*d2Psi*qfib(1)*qfib(2)
	DSfib(2,1) = 4.0d0*d2Psi*qfib(2)*qfib(1)
	DSfib(2,2) = 4.0d0*d2Psi*qfib(2)*qfib(2) + 2.0d0*dPsi

end subroutine



subroutine strainEnergyEx(energy,Ex,matPar,constLaw)
	real*8, intent(in) :: Ex , matPar(:)
	real*8 , intent(out) :: energy
	integer , intent(in) :: constLaw
	real*8 :: I4
	
	energy = 0.0d0
	
	I4 = 2.0d0*Ex + 1.0d0
	
	if(I4<1.0d0) then 
		return
	end if
	
	select case(constLaw) 
		case(1) 
			call FibreExponential(energy,I4,matPar(1:2))
		case(2) 
			call FibreLinear(energy,I4,matpar(1))
		case default
			write(0,*) constLaw, matpar(1:3), "Strain Energy case unknown"
			PAUSE
	end select

end subroutine


subroutine SxAnalytic(Sx,Ex,matPar,constLaw)
	real*8, intent(in) :: Ex , matPar(:)
	real*8 , intent(out) :: Sx
	integer , intent(in) :: constLaw
	real*8 :: I4, dPsi
	
	Sx = 0.0d0
	
	I4 = 2.0d0*Ex + 1.0d0
	
	if(I4<1.0d0) then 
		return
	end if
	
	select case(constLaw) 
		case(1) 
			call dFibreExponential(dPsi,I4,matPar(1:2))
		case(2) 
			call dFibreLinear(dPsi,I4,matpar(1))
		case default
			write(0,*) constLaw, matpar(1:3), "Strain Energy case unknown"
			PAUSE
	end select
	
	Sx = 2.0d0*dPsi

end subroutine

subroutine DSxAnalytic(DSx,Ex,matPar,constLaw)
	real*8, intent(in) :: Ex , matPar(:)
	real*8 , intent(out) :: DSx
	integer , intent(in) :: constLaw
	real*8 :: I4, d2Psi
		
	DSx = 0.0d0
	
	I4 = 2.0d0*Ex + 1.0d0
	
	if(I4<1.0d0) then 
		return
	end if
	
	select case(constLaw) 
		case(1) 
			call d2FibreExponential(d2Psi,I4,matPar(1:2))
		case(2) 
			call d2FibreLinear(d2Psi,I4,matpar(1))
		case default
			write(0,*) constLaw, matpar(1:3), "Strain Energy case unknown"
			PAUSE
	end select
	
	DSx = 4.0d0*d2Psi
	
end subroutine

end module

