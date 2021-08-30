
module materialsLib
use funcAux
implicit none
private  Delfino, MooneyRivlin, Fiber, vol, calcInvCfromF, calcInvCfromC, calcInvCMfromF, calcInvCMfromC
public StrainEnergy, StrainEnergyC

contains

!~ subroutine vol(Psi,InvC,MatPar)
!~ 	implicit none
!~ 	real*8 , intent(out) :: Psi
!~ 	real*8 , intent(in) :: InvC(:),MatPar
!~ 	real*8 :: c1,I3, J
!~ 	c1 = MatPar ! kappa
!~ 
!~ 	J = dsqrt(InvC(3)) ! J
!~ 	
!~ 	Psi = c1*(J - 1.0d0)**2.0d0
!~ 
!~ end subroutine

!~ subroutine vol(Psi,invC,MatPar)
!~ 	implicit none
!~ 	real*8 , intent(out) :: Psi
!~ 	real*8 , intent(in) :: InvC(:),MatPar
!~ 	real*8 :: c1,J
!~ 	c1 = MatPar ! kappa
!~ 	
!~ 	J = dsqrt(invC(3)) ! J
!~ 	
!~ 	Psi = 0.5d0*c1*( J*J - 1.0d0 - 2.0d0*dlog(J) )
!~ 	return
!~ end subroutine

subroutine vol(Psi,invC,MatPar)
	implicit none
	real*8 , intent(out) :: Psi
	real*8 , intent(in) :: InvC(:),MatPar
	real*8 :: c1,J
	c1 = MatPar ! kappa
	
	J = dsqrt(invC(3)) ! J
	
	Psi = 0.5d0*c1*( 0.5d0*(J*J - 1.0d0) - dlog(J) )
	return
end subroutine

!~ 
!~ subroutine vol(Psi,invC,MatPar)
!~ 	implicit none
!~ 	real*8 , intent(out) :: Psi
!~ 	real*8 , intent(in) :: InvC(:),MatPar
!~ 	real*8 :: c1, J
!~ 	c1 = MatPar ! kappa
!~ 	
!~ 	J = dsqrt(InvC(3)) ! J
	
!~ 	Psi = c1*dlog(J)**2.0d0
!~ 	return
!~ end subroutine


!~ subroutine vol2(Psi,C,MatPar)
!~ 	implicit none
!~ 	real*8 , intent(out) :: Psi
!~ 	real*8 , intent(in) :: C(:,:),MatPar(2)
!~ 	real*8 :: c1,c2,I3
!~ 	
!~ 	c1 = 2.0d0*MatPar(1) ! kappa
!~ 	c2 = MatPar(2) ! kappa
!~ 	
!~ 	write(0,*) "New volumetric penalty function"
!~ 	
!~ 	call calcI3(I3,C)
!~ 	I3 = dlog(dsqrt(I3)) ! J
!~ 	
!~ 	Psi = -c1*I3 + c2*I3*I3
!~ 	return
!~ end subroutine

subroutine Delfino(Psi,invC,MatPar)
! Delfino et al. [10] proposed an (isotropic) rubber-like potential for carotid arteries
          ! [10] A. Delfino, N. Stergiopulos, J.E. Moore and J.-J. Meister, 
          ! Residual strain effects on the stress field in a thick wall
          ! finite element model of the human carotid bifurcation. 
          ! J. Biomech. 30 (1997) 777786.
	implicit none
	real*8 , intent(out) :: Psi
	real*8 , intent(in) :: invC(:),MatPar(2)
	real*8 :: c1,c2,c3,I1bar,J,rI33
	c1 = MatPar(1) ! a
	c2 = MatPar(2)*0.5d0 ! b/2
	
	J = dsqrt(InvC(3)) ! J
	rI33 = J**(-2.0D0/3.0D0)
	I1bar = rI33*InvC(1) ! I1bar
	
	Psi = 0.5d0*(c1/c2)*( Dexp(c2*(I1bar-3.0d0))-1.0d0 )
	
end subroutine

!~ subroutine MooneyRivlin(Psi,invC,MatPar)
!~  !Mooney Rivlin Bathe, K.J., Finite Element Procedures, pp.592; see also 6.4 pp.561
!~  !Crisfield M.A. Vol II, Chap 13, pp 62. See also, Section 10.3 pp 07 
!~ 	implicit none
!~ 	real*8 , intent(out) :: Psi
!~ 	real*8 , intent(in) :: invC(:),MatPar(2)
!~ 	real*8 :: c1,c2,I1bar,I2bar,J,rI33
!~ 	c1 = MatPar(1)
!~ 	c2 = MatPar(2)
!~ 
!~ 	J = dsqrt(InvC(3)) ! J
!~ 	rI33 = J**(-2.0D0/3.0D0)
!~ 	I1bar = rI33*InvC(1) ! I1bar
!~ 	I2bar = rI33*rI33*invC(2) ! I2bar
!~ 	Psi = c1*(I1bar-3.0d0) + c2*(I2bar-3.0d0)
!~ 	return
!~ end subroutine


subroutine MooneyRivlin(Psi,invC,MatPar)
 !Mooney Rivlin Bathe, K.J., Finite Element Procedures, pp.592; see also 6.4 pp.561
 !Crisfield M.A. Vol II, Chap 13, pp 62. See also, Section 10.3 pp 07 
	implicit none
	real*8 , intent(out) :: Psi
	real*8 , intent(in) :: invC(:),MatPar(2)
	real*8 :: c1,c2,I1bar,I2bar,J,rI33
	c1 = MatPar(1)
	c2 = MatPar(2)

	J = dsqrt(InvC(3)) ! J
	rI33 = J**(-2.0D0/3.0D0)
	I1bar = rI33*InvC(1) ! I1bar
	I2bar = rI33*rI33*invC(2) ! I2bar
	Psi = c1*(I1bar-3.0d0) + c2*(I2bar-3.0d0)
	return
end subroutine

subroutine Fiber(Psi,invCM,MatPar) ! see my dissertation
	implicit none
	real*8 , intent(out) :: Psi
	real*8 , intent(in) :: invCM(:),MatPar(2) ! MatPar is only the material coeficient, a0 was splitted
 	real*8 :: c1,c2,I4
	c1 = MatPar(1) ! k1
	c2 = MatPar(2) ! k2
	
	I4 = invCM(2) !! already incompressible
	
	Psi = 0.5d0*(c1/c2)*( Dexp(c2*(I4-1.0d0)*2.0d0)-1.0d0 ) 
	
end subroutine 

subroutine calcInvCMfromF(invCM,F,beta)
	real*8 , intent(in) :: F(:,:), beta
	real*8 , intent(out) :: InvCM(:)
	real*8 :: v(size(F,1)), J , a0(size(F,1))
		
	a0(1) = dcos(beta)
	a0(2) = dsin(beta)
	
	v = matmul(F,a0)
	InvCM(1) = dot_product(v,v)
	
	call calcI3(J,F)
	
	InvCM(2) = (J**(-2.0D0/3.0D0))*InvCM(1) ! incompressible
	
end subroutine


subroutine calcInvCMfromC(invCM,C,beta)
	real*8 , intent(in) :: C(:,:), beta
	real*8 , intent(out) :: InvCM(:)
	real*8 :: v(size(C,1)), J, a0(size(C,1))
	
	a0(1) = dcos(beta)
	a0(2) = dsin(beta)
	
	v = matmul(C,a0)	
	InvCM(1) = dot_product(v,a0)
	
	call calcI3(J,C)
	J = dsqrt(J)
	
	InvCM(2) = (J**(-2.0D0/3.0D0))*InvCM(1) ! incompressible
	
end subroutine

	
subroutine calcInvCfromC(invC,C)
	real*8 , intent(in) :: C(:,:)
	real*8 , intent(out) :: InvC(:)
		
	call calcI1(InvC(1),C)
	call calcI2(InvC(2),C)
	call calcI3(InvC(3),C)
	
end subroutine

	
subroutine calcInvCfromF(invC,F)
	real*8 , intent(in) :: F(:,:)
	real*8 , intent(out) :: InvC(:)
	real*8 :: C(size(F,1),size(F,2))
		
	C = matmul(transpose(F),F)
	
	call calcI1(InvC(1),C)
	call calcI2(InvC(2),C)
	call calcI3(InvC(3),C)
	
end subroutine

subroutine strainEnergy(energy,F,matPar,constLaw)
	real*8, intent(in) :: F(:,:) , matPar(:)
	real*8 , intent(out) :: energy
	integer , intent(in) :: constLaw
	real*8 :: rAux, energyTemp, invC(3)
	
	if(constLaw<10) then 
		call calcInvCfromF(invC,F)
	else
		call calcInvCMfromF(invC,F,matPar(3))
	end if
	
	energy = 0.0d0
	energyTemp = 0.0d0

	selectcase(constLaw) 
		case(3) 
			call Delfino(energy,invC,matPar(1:2))
			call vol(energyTemp,invC,matPar(3))
			energy = energy + energyTemp
		case(4)
			call MooneyRivlin(energy,invC,matPar(1:2))
            call vol(energyTemp,invC,matPar(3))
            energy = energy + energyTemp
		case(13)
			call Fiber(energy,invC,matPar(1:2))
		case default
			write(0,*) constLaw, matpar(1:3), "Strain Energy case unknown"
			PAUSE
	end select

end subroutine

subroutine strainEnergyC(energy,C,matPar,constLaw)
	real*8, intent(in) :: C(:,:) , matPar(:)
	real*8 , intent(out) :: energy
	integer , intent(in) :: constLaw
	real*8 :: rAux, energyTemp, invC(3)
	
	if(constLaw<10) then 
		call calcInvCfromC(invC,C)
	else
		call calcInvCMfromC(invC,C,matPar(3))
	end if
	
	energy = 0.0d0
	energyTemp = 0.0d0

	selectcase(constLaw) 
		case(3) 
			call Delfino(energy,invC,matPar(1:2))
			call vol(energyTemp,invC,matPar(3))
			energy = energy + energyTemp
		case(4)
			call MooneyRivlin(energy,invC,matPar(1:2))
            call vol(energyTemp,invC,matPar(3))
            energy = energy + energyTemp
		case(13)
			call Fiber(energy,invC,matPar(1:2))
		case default
			write(0,*) constLaw, matpar(1:3), "Strain Energy case unknown"
			PAUSE
	end select

end subroutine


end module

