module damageModelsLib
implicit none

private  Hexp, HexpReg, HexpRegSqrt, RmaxEnergy, RmaxEnergySqrt, HBlanco, HSanchez, derDamageSanchez
public chooseHfunc, chooseRfunc

abstract interface
	subroutine HfuncLike(H,HPar,r,q) 
		implicit none
		real*8 , intent(in) :: r,q,HPar(:)
		real*8 , intent(out) :: H
	end subroutine
	
	subroutine RfuncLike(r,RPar,rOld,energy) 
		implicit none
		real*8 , intent(in) :: rOld,energy,RPar(:)
		real*8 , intent(out) :: r
	end subroutine

	subroutine damageFuncLike(d,HPar,r,q) 
		implicit none
		real*8 , intent(in) :: r,q,HPar(:)
		real*8 , intent(out) :: d
	end subroutine
	
	subroutine derDamageFuncLike(derD,HPar,r,q) 
		implicit none
		real*8 , intent(in) :: r,q,HPar(:)
		real*8 , intent(out) :: derD
	end subroutine
end interface

contains

subroutine chooseHfunc(Hfunc,damageFunc,derDamageFunc,HType)
	procedure (HfuncLike), pointer :: Hfunc 
	procedure (damageFuncLike), pointer :: damageFunc 
	procedure (derDamageFuncLike), pointer :: derDamageFunc 
	integer, intent(in) :: HType

	selectcase(HType) 
		case(1) 
			Hfunc => Hexp
		case(2) 
			Hfunc => HexpReg
		case(3) 
			Hfunc => HexpRegSqrt
		case(4)		
			Hfunc => HBlanco
		case(5)		
			Hfunc => HSanchez  
		case(11) 
			damageFunc => DamageExp
			derDamageFunc => derDamageExp
		case(12) 
			damageFunc => DamageExpReg
			derDamageFunc => derDamageExpReg
		case(13) 
			damageFunc => DamageExpRegSqrt
			derDamageFunc => derDamageExpRegSqrt
		case default 
		
		write(*,*) "Damage Model Not found Htype = " , Htype
		pause
	end select

end subroutine

subroutine chooseRfunc(Rfunc,RType)
	procedure (RfuncLike), pointer :: Rfunc 
	integer, intent(in) :: RType

	select case(RType) 
		case(1) 
			Rfunc => RmaxEnergy
		case(2) 
			Rfunc => RmaxEnergySqrt
		case default 
			write(*,*) "Damage Model Not found Rtype = " , Rtype 
			pause
	end select

end subroutine

subroutine RmaxEnergy(r,RPar,rOld,energy) 
	real*8 , intent(in) :: rOld,energy,RPar(:)
	real*8 , intent(out) :: r
	real*8 :: r0 , rAux
	
	r0 = RPar(1)
	rAux = max(r0,rOld)
	r = max(rAux,energy)
	
end subroutine

subroutine RmaxEnergySqrt(r,RPar,rOld,energy) 
	real*8 , intent(in) :: rOld,energy,RPar(:)
	real*8 , intent(out) :: r
	real*8 :: r0 , rAux
	
	r0 = RPar(1)

	rAux = max(r0,rOld)
	r = max(rAux,dsqrt(2.0d0*energy))
	
end subroutine

subroutine Hexp(H,HPar,r,q) 
	real*8 , intent(in) :: r, q, HPar(:)
	real*8 , intent(out) :: H 
	real*8 :: dinf, beta, expmrbeta
	
	dinf = HPar(1)
	beta = HPar(2)
	
	expmrbeta = dexp(-r/beta)
	H = 1.0d0 - dinf*(1.0d0 - expmrbeta) - r*dinf*expmrbeta/beta

end subroutine

subroutine damageExp(d,HPar,r,q) 
	implicit none
	real*8 , intent(in) ::  r, q, HPar(:)
	real*8 , intent(out) :: d
	real*8 :: beta, dinf
	
	dinf = HPar(1)
	beta = HPar(2)
	
	d = dinf*(1.0d0 - dexp(-r/beta))

end subroutine

subroutine derDamageExp(derD,HPar,r,q) 
	implicit none
	real*8 , intent(in) ::  r,  q, HPar(:)
	real*8 , intent(out) :: derD
	real*8 :: beta, dinf
	
	dinf = HPar(1)
	beta = HPar(2)
	
	derD = dinf*dexp(-r/beta)/beta

end subroutine

subroutine HexpReg(H,HPar,r,q) 
	real*8 , intent(in) :: r, q, HPar(:)
	real*8 , intent(out) :: H 
	real*8 :: Gf, r0, beta, he0
	
	Gf = HPar(1)
	r0 = HPar(2)
	he0 = HPar(3)
	
	beta = (Gf/he0) - r0
	
	H = dexp(-(r-r0)/beta)*(1.0d0 - r/beta)

end subroutine

subroutine damageExpReg(d,HPar,r,q) 
	implicit none
	real*8 , intent(in) ::  r, q, HPar(:)
	real*8 , intent(out) :: d
	real*8 :: Gf, r0, beta, he0
	
	Gf = HPar(1)
	r0 = HPar(2)
	he0 = HPar(3)

	beta = (Gf/he0) - r0
	
	d = 1.0d0 - dexp(-(r-r0)/beta)

end subroutine

subroutine derDamageExpReg(derD,HPar,r,q) 
	implicit none
	real*8 , intent(in) ::  r, q, HPar(:)
	real*8 , intent(out) :: derD
	real*8 :: Gf, r0, beta, he0
	
	Gf = HPar(1)
	r0 = HPar(2)
	he0 = HPar(3)
	
	beta = (Gf/he0) - r0
	
	derD = dexp(-(r-r0)/beta)/beta

end subroutine


subroutine HexpRegSqrt(H,HPar,r,q) 
	implicit none
	real*8 , intent(in) :: r, q, HPar(:)
	real*8 , intent(out) :: H
	real*8 :: Gf, beta,r0, he0
	
	Gf = HPar(1)
	r0 = HPar(2)
	he0 = HPar(3) 

	beta = 0.5d0*(-r0 + dsqrt(4.0d0*Gf/he0  - r0*r0)) 	
 
	H = dexp(-(r-r0)/beta)*(1.0d0 - r/beta)

end subroutine

subroutine damageExpRegSqrt(d,HPar,r,q) 
	implicit none
	real*8 , intent(in) ::  r, q, HPar(:)
	real*8 , intent(out) :: d
	real*8 :: Gf, r0, beta, he0
	
	Gf = HPar(1)
	r0 = HPar(2)
	he0 = HPar(3)
	
!~ 	beta = 0.5d0*(-r0 + dsqrt(4.0d0*Gf/he0  - r0*r0)) 	
	beta = Gf
	
	d = 1.0d0 - dexp(-(r-r0)/beta)

end subroutine

subroutine derDamageExpRegSqrt(derD,HPar,r,q) 
	implicit none
	real*8 , intent(in) ::  r, q, HPar(:)
	real*8 , intent(out) :: derD
	real*8 :: Gf, r0, beta, he0
	
	Gf = HPar(1)
	r0 = HPar(2)
	he0 = HPar(2)
	 	
!~ 	beta = 0.5d0*(-r0 + dsqrt(4.0d0*Gf/he0  - r0*r0)) 	
	beta = Gf 	
	
	derD = dexp(-(r-r0)/beta)/beta

end subroutine

subroutine HBlanco(H,HPar,r,q)  
	implicit none
	real*8 , intent(in) :: r, q, HPar(:) 
	real*8 , intent(out) :: H
	real*8 :: chi, Gf, q0, he0
	
	Gf = HPar(1)
	q0 = HPar(2)
	he0 = HPar(3)
!~ 	chi = HPar(3)
	chi = 1.0d0
	
	q0 = (q0**(2.0d0-chi))/(2.0d0-chi)
	
	H = -q0*(r**chi)*he0/Gf

end subroutine

subroutine HSanchez(H,HPar,r,q) 
	implicit none
	real*8 , intent(in) :: r, q, HPar(:)
	real*8 , intent(out) :: H
	real*8 :: Gf, r0, he0, H0
	
	Gf = HPar(1)
	r0 = HPar(2)
	he0 = HPar(3) 
	
!~ 	H0 = -(he0*r0**2.0d0)/Gf
	
	H0 = (0.5d0 - Gf/(he0*r0**2.0d0))**(-1.0d0)
		 	
	H = H0*dexp(H0*(r-r0)/r0)

end subroutine

subroutine derDamageSanchez(derD,HPar,r,q) 
	implicit none
	real*8 , intent(in) :: r, q, HPar(:)
	real*8 , intent(out) :: derD
	real*8 :: H
	
	call HSanchez(H,Hpar,r,q)
	
	derD = (-H*r+q)/(r*r)

end subroutine

end module
