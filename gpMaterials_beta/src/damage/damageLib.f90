module damageLib
use funcAux

implicit none

! for compatibility is Real*8, integer, integer. Det to me is local, not shared in common
real*8 :: dummy
integer:: LengthParam,LengthJParam
COMMON dummy,LengthParam,LengthJParam

private  Hexp, HexpReg, HexpRegSqrt, RmaxEnergy, RmaxEnergySqrt, HBlanco, HSanchez, derDamageSanchez, damageImplex
public Hfunc, Rfunc, damageFunc, derDamageFunc, chooseHfunc, chooseRfunc, &
	   loadDamage, putDamage,  damageModifyTangentF, damageUpdate, damageModifyTangentNew, damageModifyTangentDouble, copyDamage

integer , parameter :: Idamage = 1 , Ird = 2 , Iqd = 3,  IrdDot = 4, IderDamage = 5, NdamageParam = 5
real*8 :: derDamage,damage,rd,qd, rdDot
real*8 :: derDamageNew,damageNew,rdNew,qdNew, rdDotNew

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

procedure (HfuncLike), pointer :: Hfunc => null ()
procedure (RfuncLike), pointer :: Rfunc => null ()
procedure (damageFuncLike), pointer :: damageFunc => null ()
procedure (derDamageFuncLike), pointer :: derDamageFunc => null ()

contains

subroutine damageUpdate(commonPar,Param,energy,iShiftDamageCommonPar,iShiftDamageParam, nG)
	integer, intent(in) :: iShiftDamageParam,  iShiftDamageCommonPar, nG
	real*8, intent(in) :: commonPar(*), energy 
	real*8, intent(inout) :: Param(*)
	integer :: Htype, Rtype
	real*8 :: Rpar(1) , Hpar(3) , H
	 		
	Rtype =  nint(commonPar(1 + iShiftDamageCommonPar))
	Rpar(1) = CommonPar(2 + iShiftDamageCommonPar)
	Htype =  nint(commonPar(3 + iShiftDamageCommonPar))
	HPar(1:3) = CommonPar(4 + iShiftDamageCommonPar: 6 + iShiftDamageCommonPar)
	
	call chooseRfunc(Rtype)
	call chooseHfunc(Htype) !! DerDamageFunc damageFunc as well if available

	call loadDamage(Param,iShiftDamageParam,nG,.false.)

	if(qd == 0.0d0 .and. rd == 0.0d0) then 
		rd = Rpar(1)
		qd = Rpar(1)
	end if
	
	call Rfunc(rdNew,RPar, rd, energy)
!~ 	call Hfunc(H, HPar, 0.5d0*(rd+rdNew) , qd, he) !! integration errors are bigger otherwise
	call Hfunc(H, HPar, rd , qd) !! integration errors are bigger otherwise

	rdDotNew = rdNew - rd
	qdNew = qd + H*rdDotNew
	
	if(qdNew < 0.0d0) qdNew = 0.000001d0
	
	damageNew = 1.0d0 - qdNew/rdNew !! integrating
		
	call Hfunc(H, HPar, rdNew , qdNew)
	
	derDamageNew = (-H*rdNew + qdNew)/(rdNew*rdNew)
!~ 	call derDamageFunc(derDamageNew, HPar, rdNew , qdNew)
	
	if(Rtype==2) derDamageNew = derDamageNew/rdNew
	
	call putDamage(Param,iShiftDamageParam,nG)
	
!~ 	if( LengthParam>60 .and. iShiftDamageParam == 0) then
!~ 		write(*,*) "========== damage copied ============="
!~ 		call copyDamage(Param,iShiftDamageParam,10)
!~ 		call copyDamage(Param,iShiftDamageParam,17)
!~ 		call copyDamage(Param,iShiftDamageParam,24)
!~ 		call copyDamage(Param,iShiftDamageParam,31)
!~ 	end if
	
!~ 	write(*,*) "Debug damage = "
!~ 	write(*,*) Rtype,RPar,Htype,HPar, energy
!~ 	write(*,*) damage, rd, qd, damageDot, rdDot, qdDot
!~ 	write(*,*) damageNew, rdNew, qdNew, damageDotNew, rdDotNew, qdDotNew
!~ 	PAUSE

end subroutine

function damageImplex(commonPar,iShiftDamageCommonPar) result(damageTilde)
	Real*8 , intent(in) :: commonPar(*)
	integer , intent(in) :: iShiftDamageCommonPar
	
	real*8 :: damageTilde, qTilde, rTilde, H, HPar(3), RPar(1)
	integer :: Rtype, Htype

	Rtype =  nint(commonPar(1 + iShiftDamageCommonPar))
	Rpar(1) = CommonPar(2 + iShiftDamageCommonPar)
	Htype =  nint(commonPar(3 + iShiftDamageCommonPar))
	HPar(1:3) = CommonPar(4 + iShiftDamageCommonPar: 6 + iShiftDamageCommonPar)
		
	call chooseRfunc(Rtype)
	call chooseHfunc(Htype)

	rTilde = rd + rdDot
	call Hfunc(H, HPar, rd , qd)
	qTilde = qd + H*(rTilde - rd)
	
	if(rTilde < 0.00000001d0) then
		rTilde = Rpar(1)
		qTilde = Rpar(1)
	end if  
	
	damageTilde = 1.0d0 - qTilde/rTilde
		
end function

subroutine damageModifyTangentDouble(D,DB,Ppk,commonPar,Param,NdimE, &
							iShiftDamageCommonPar,iShiftDamageParam,nG,IdamageModify)
			
	Real*8 , intent(inout):: DB(NdimE,NdimE,NdimE,NdimE),D(NdimE,NdimE,NdimE,NdimE),Ppk(NdimE,NdimE)
	Real*8 , intent(in) :: commonPar(*),Param(*)
	integer, intent(in) :: NdimE, iShiftDamageParam, iShiftDamageCommonPar, nG
	integer , intent(in) :: IdamageModify
	logical :: flagDamageDer, flagDamageNew
	integer :: i,j,k,l
	real*8 :: alpha , derdamageOld
	
	select case(IdamageModify)
		case(-1)
			return !!! Damage is not defined
		case(0,4)  !! 4 is for simplex
			flagDamageDer = .false.
			flagDamageNew = .false.
		case(1)	
			flagDamageDer = .true.
			flagDamageNew = .true.
			alpha = 0.5d0
		case(2) !! it could be useful
			flagDamageDer = .true.
			flagDamageNew = .false.
			alpha = 0.5d0
		case(3,5) !! 5 is the case where update inside finiteStrain but no tangent modification 
			flagDamageDer = .false.
			flagDamageNew = .true.			
		case default
			write(*,*) "IdamageModify should be -1,0,1,2,3, 4 or 5"
			STOP
	end select
	 
	
	call loadDamage(Param,iShiftDamageParam,nG,.false.)
	
	derDamageOld = derDamage
	
	call loadDamage(Param,iShiftDamageParam,nG,flagDamageNew)
		
	if(IdamageModify == 4) damage = damageImplex(CommonPar,iShiftDamageCommonPar)
		
	D = (1.0d0-damage)*D 
	DB = 0.0d0
	
	if(rdDot>0.0d0 .and. derDamage>0.0d0 .and. flagDamageDer) then 
			
		do i = 1,NdimE
		do j = 1,NdimE
		do k = 1,NdimE
		do l = 1,NdimE
				D(i,j,k,l) = D(i,j,k,l) - alpha*derDamage*Ppk(i,j)*Ppk(k,l)
				DB(i,j,k,l) = - (1.0d0 - alpha)*derDamage*Ppk(i,j)*Ppk(k,l)
		end do
		end do
		end do
		end do
	end if
	
	Ppk = (1.0d0-damage)*Ppk
	
end subroutine



subroutine damageModifyTangentNew(D,Ppk,commonPar,Param,NdimE,&
			iShiftDamageCommonPar,iShiftDamageParam,nG,IdamageModify)
			
	Real*8 , intent(inout):: D(NdimE,NdimE,NdimE,NdimE),Ppk(NdimE,NdimE)
	Real*8 , intent(in) :: commonPar(*),Param(*)
	integer, intent(in) :: NdimE, iShiftDamageParam, iShiftDamageCommonPar, IdamageModify, nG
	logical :: flagDamageDer, flagDamageNew
	integer :: i,j,k,l
	
	select case(IdamageModify)
		case(-1)
			return !!! Damage is not defined
		case(0,4)  !! 4 is for simplex
			flagDamageDer = .false.
			flagDamageNew = .false.
		case(1)	
			flagDamageDer = .true.
			flagDamageNew = .true.
		case(2) !! it could be useful
			flagDamageDer = .true.
			flagDamageNew = .false.
		case(3,5) !! 5 is the case where update inside finiteStrain but no tangent modification 
			flagDamageDer = .false.
			flagDamageNew = .true.			
		case default
			write(*,*) "IdamageModify should be -1,0,1,2,3, 4 or 5"
			STOP
	end select
	 
	call loadDamage(Param,iShiftDamageParam,nG,flagDamageNew)
		
	if(IdamageModify == 4) damage = damageImplex(CommonPar,iShiftDamageCommonPar)
		
	D = (1.0d0-damage)*D 
	
	if(rdDot>0.0d0 .and. derDamage>0.0d0 .and. flagDamageDer) then 
			
		do i = 1,NdimE
		do j = 1,NdimE
		do k = 1,NdimE
		do l = 1,NdimE
				D(i,j,k,l) = D(i,j,k,l) - derDamage*Ppk(i,j)*Ppk(k,l)
		end do
		end do
		end do
		end do
	end if
	
	Ppk = (1.0d0-damage)*Ppk
	
end subroutine

!~ 
subroutine damageModifyTangentF(D,Ppk,commonPar,Param,NdimE,&
			iShiftDamageCommonPar,iShiftDamageParam,nG,flagDamageDer,flagDamageNew)
			
	Real*8 , intent(inout):: D(NdimE,NdimE,NdimE,NdimE),Ppk(NdimE,NdimE)
	Real*8 , intent(in) :: commonPar(*),Param(*)
	integer, intent(in) :: NdimE, iShiftDamageParam, iShiftDamageCommonPar, nG
	logical , intent(in) :: flagDamageDer, flagDamageNew
	real*8 :: derDamageLocal, HPar(3), RPar(1), rTilde, qTilde, damageTilde,H
	integer :: i,j,k,l, Htype, Rtype
	logical :: isImplex = .false.
	 
	call loadDamage(Param,iShiftDamageParam,nG,flagDamageNew)
	
	if(isImplex .and. (.not. flagDamageNew) .and. (.not. flagDamageDer) ) then
		Rtype =  nint(commonPar(1 + iShiftDamageCommonPar))
		Rpar(1) = CommonPar(2 + iShiftDamageCommonPar)
		Htype =  nint(commonPar(3 + iShiftDamageCommonPar))
		HPar(1:3) = CommonPar(4 + iShiftDamageCommonPar: 6 + iShiftDamageCommonPar)
		
		call chooseRfunc(Rtype)
		call chooseHfunc(Htype) !! DerDamageFunc damageFunc as well if available
		
!~ 		!!call loadDamage(Param,iShiftDamageParam,.false.)
		rTilde = rd + rdDot
		call Hfunc(H, HPar, rd , qd)
		qTilde = qd + H*(rTilde - rd)
		damageTilde = 1.0d0 - qTilde/rTilde
		
		D = (1.0d0-damageTilde)*D 	
	else
		D = (1.0d0-damage)*D 	
	end if
	
	if(rdDot>0.0d0 .and. flagDamageDer) then 
		
		Rtype =  nint(commonPar(1 + iShiftDamageCommonPar))
		Rpar(1) = CommonPar(2 + iShiftDamageCommonPar)
		Htype =  nint(commonPar(3 + iShiftDamageCommonPar))
		HPar(1:3) = CommonPar(4 + iShiftDamageCommonPar: 6 + iShiftDamageCommonPar)
		
		call chooseRfunc(Rtype)
		call chooseHfunc(Htype) !! DerDamageFunc damageFunc as well if available
	
		call derDamageFunc(derDamageLocal,HPar,rd,qd)
		
		write(*,*) "comparison = " , derDamage, derDamageLocal, derDamageLocal/rd
!~ 		PAUSE
		
		do i = 1,NdimE
		do j = 1,NdimE
		do k = 1,NdimE
		do l = 1,NdimE
			if(Rtype == 1) then
				D(i,j,k,l) = D(i,j,k,l) - derDamageLocal*Ppk(i,j)*Ppk(k,l)
			else if(Rtype == 2) then
				D(i,j,k,l) = D(i,j,k,l) - (derDamageLocal/rd)*Ppk(i,j)*Ppk(k,l)
			else
				write(*,*) "Rtype should be 1 or 2"
				STOP
			end if
		end do
		end do
		end do
		end do
	end if
		
	if(isImplex .and. (.not. flagDamageNew) .and. (.not. flagDamageDer) ) then
		Ppk = (1.0d0-damageTilde)*Ppk
	else 
!~ 		!!write(*,*) "damage uncoupled = " , damage
		Ppk = (1.0d0-damage)*Ppk
	end if

end subroutine

subroutine chooseHfunc(HType)
	integer, intent(in) :: HType

	selectcase(HType) 
		case(1) 
			Hfunc => Hexp
			damageFunc => DamageExp
			derDamageFunc => derDamageExp
		case(2) 
			Hfunc => HexpReg
			damageFunc => DamageExpReg
			derDamageFunc => derDamageExpReg
		case(3) 
			Hfunc => HexpRegSqrt
			damageFunc => DamageExpRegSqrt
			derDamageFunc => derDamageExpRegSqrt
		case(4)		
			Hfunc => HBlanco
!~ 			damageFunc => DamageBlanco
!~ 			derDamageFunc => derDamageBlanco
		case(5)		
			Hfunc => HSanchez  
!~ 			damageFunc => DamageSanchez
			derDamageFunc => derDamageSanchez
		case default 
			Hfunc => Hexp
			damageFunc => DamageExp
			derDamageFunc => derDamageExp
	end select

end subroutine

subroutine chooseRfunc(RType)
	integer, intent(in) :: RType

	selectcase(RType) 
		case(1) 
			Rfunc => RmaxEnergy
		case(2) 
			Rfunc => RmaxEnergySqrt
		case default 
			Rfunc => RmaxEnergy
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

subroutine damageExp(d,HPar,r,q) 
	implicit none
	real*8 , intent(in) ::  r, q, HPar(:)
	real*8 , intent(out) :: d
	real*8 :: beta, dinf
	
	dinf = HPar(1)
	beta = HPar(2)
	
	d = dinf*(1.0d0 - dexp(-r/beta))

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

subroutine damageExpRegSqrt(d,HPar,r,q) 
	implicit none
	real*8 , intent(in) ::  r, q, HPar(:)
	real*8 , intent(out) :: d
	real*8 :: Gf, r0, beta, he0
	
	Gf = HPar(1)
	r0 = HPar(2)
	he0 = HPar(3)
	
	beta = 0.5d0*(-r0 + dsqrt(4.0d0*Gf/he0  - r0*r0)) 	

	d = 1.0d0 - dexp(-(r-r0)/beta)

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

subroutine derDamageExpRegSqrt(derD,HPar,r,q) 
	implicit none
	real*8 , intent(in) ::  r, q, HPar(:)
	real*8 , intent(out) :: derD
	real*8 :: Gf, r0, beta, he0
	
	Gf = HPar(1)
	r0 = HPar(2)
	he0 = HPar(2)
	 	
	beta = 0.5d0*(-r0 + dsqrt(4.0d0*Gf/he0  - r0*r0)) 	
	
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
	real*8 :: Gf, alpha, r0, he0
	
	Gf = HPar(1)
	r0 = HPar(2)
	he0 = HPar(3) 
	
	alpha = he0*r0/Gf
		 	
	H = -alpha*r0*dexp(-alpha*(r-r0))

end subroutine

subroutine derDamageSanchez(derD,HPar,r,q) 
	implicit none
	real*8 , intent(in) :: r, q, HPar(:)
	real*8 , intent(out) :: derD
	real*8 :: H
	
	call HSanchez(H,Hpar,r,q)
	
	derD = (-H*r+q)/(r*r)

end subroutine

subroutine loadDamage(Param,iShiftDamageParam,nG,inNew)
	Real*8 , intent(in) :: Param(*)
	integer , intent(in) :: iShiftDamageParam, nG
	logical , intent(in) :: inNew
	integer :: iShiftDamageParamMod
	
	if(inNew) then
		iShiftDamageParamMod = LengthParam + (ng-1)*NdamageParam + iShiftDamageParam 
	else
		iShiftDamageParamMod = (ng-1)*NdamageParam + iShiftDamageParam
	end if
	
	damage = Param(iShiftDamageParamMod + Idamage)
	rd = Param(iShiftDamageParamMod + Ird)
	qd = Param(iShiftDamageParamMod + Iqd)
	rdDot = Param(iShiftDamageParamMod + IrdDot)
	derDamage = Param(iShiftDamageParamMod + IderDamage)

end subroutine

subroutine putDamage(Param,iShiftDamageParam, nG)
	Real*8 , intent(out) :: Param(*)
	integer , intent(in) :: iShiftDamageParam, nG
	integer :: iShiftDamageParamMod

    iShiftDamageParamMod = LengthParam + (ng-1)*NdamageParam + iShiftDamageParam 
     
	Param(iShiftDamageParamMod + Idamage) = damageNew 
	Param(iShiftDamageParamMod + Ird) = rdNew
	Param(iShiftDamageParamMod + Iqd) = qdNew 
	Param(iShiftDamageParamMod + IrdDot) = rdDotNew 
	Param(iShiftDamageParamMod + IderDamage) = derDamageNew
end subroutine

subroutine copyDamage(Param,iShiftDamageParam, nG, iShiftCopy)
	Real*8 , intent(out) :: Param(*)
	integer , intent(in) :: iShiftDamageParam, iShiftCopy, nG
	integer :: iShiftIn, iShiftOut
	
	iShiftOut = iShiftDamageParam + (ng-1)*NdamageParam + iShiftCopy
	iShiftIn = iShiftDamageParam + (ng-1)*NdamageParam
	
	Param(iShiftOut + Idamage) = Param(iShiftIn + Idamage)
	Param(iShiftOut + Ird) = Param(iShiftIn + Ird)
	Param(iShiftOut + Iqd) = Param(iShiftIn + Iqd)
	Param(iShiftOut + IrdDot) = Param(iShiftIn + IrdDot)
	Param(iShiftOut + IderDamage) = Param(iShiftIn + IderDamage)
	
end subroutine

end module

