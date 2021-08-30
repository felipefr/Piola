module damageNewLib
use funcAux
use damageModelsLIb


implicit none

! for compatibility is Real*8, integer, integer. Det to me is local, not shared in common
real*8 :: dummy
integer:: LengthParam,LengthJParam
COMMON dummy,LengthParam,LengthJParam

private qBasedUpdate, directBasedUpdate
public writeDamage , damageUpdate, damageModifyTangent, writeDamageFromParam

integer , parameter :: Idamage = 1 , Ird = 2 , Iqd = 3,  IrdDot = 4, IderDamage = 5, NdamageParam = 5, iShiftDamageParam = 0

type damageState
	real*8 :: damage, rd,qd , rdDot,derDamage
	
	contains
	procedure :: loadDamage, putDamage
end type

contains

subroutine qBasedUpdate(Hfunc, Rfunc, sd , sdNew, Rpar, Hpar, energy, Htype, Rtype)
	type(damageState), intent(inout) :: sd , sdNew 
	integer, intent(in) :: Htype, Rtype
	real*8, intent(in) :: Rpar(1) , Hpar(3), energy 
	procedure (HfuncLike), pointer , intent(in) :: Hfunc
	procedure (RfuncLike), pointer , intent(in):: Rfunc
	real*8 ::  H


	call Rfunc(sdNew%rd,RPar, sd%rd, energy)
	call Hfunc(H, HPar, 0.5d0*(sd%rd+sdNew%rd) , sd%qd) !! integration errors are bigger otherwise
!~ 	call Hfunc(H, HPar, sd%rd , sd%qd) !! integration errors are bigger otherwise

	sdNew%rdDot = sdNew%rd - sd%rd
	sdNew%qd = sd%qd + H*sdNew%rdDot
	
	if(sdNew%qd < 0.0d0) sdNew%qd = 0.000001d0
	
	sdNew%damage = 1.0d0 - sdNew%qd/sdNew%rd !! integrating
		
	call Hfunc(H, HPar, sdNew%rd , sdNew%qd)
	
	sdNew%derDamage = (-H*sdNew%rd + sdNew%qd)/(sdNew%rd**2.0d0)
	
end subroutine

subroutine directBasedUpdate(damageFunc, derDamageFunc, Rfunc,  sd , sdNew, Rpar, Hpar, energy, Htype, Rtype)
	type(damageState), intent(inout) :: sd , sdNew 
	integer, intent(in) :: Htype, Rtype
	real*8, intent(in) :: Rpar(1) , Hpar(3), energy 
	procedure (RfuncLike), pointer , intent(in):: Rfunc
	procedure (damageFuncLike), pointer , intent(in):: damageFunc 
	procedure (derDamageFuncLike), pointer , intent(in) :: derDamageFunc 
	real*8 ::  H

	call Rfunc(sdNew%rd,RPar, sd%rd, energy)

	sdNew%rdDot = sdNew%rd - sd%rd
	
	call damageFunc(sdNew%damage,HPar,sdNew%rd,0.0d0) !! last is q
	call derDamageFunc(sdNew%derDamage,HPar,sdNew%rd,0.0d0) !! last is q
	
end subroutine

subroutine damageUpdate(damagePar,Param,energy, nG)
	integer, intent(in) :: nG
	real*8, intent(in) :: damagePar(:), energy 
	real*8, intent(inout) :: Param(*)
	integer :: Htype, Rtype
	real*8 :: Rpar(1) , Hpar(3) , H
	type(damageState) :: sd , sdNew 
	procedure (HfuncLike), pointer :: Hfunc => null ()
	procedure (RfuncLike), pointer :: Rfunc => null ()
	procedure (damageFuncLike), pointer :: damageFunc => null ()
	procedure (derDamageFuncLike), pointer :: derDamageFunc => null ()

	call loadDamagePar(Hpar,Rpar,Htype,Rtype,damagePar)
	
	call chooseRfunc(Rfunc,Rtype)
	call chooseHfunc(Hfunc,damageFunc,derDamageFunc,Htype) !! DerDamageFunc damageFunc as well if available

	call sd%loadDamage(Param, nG,.false.)

	if(sd%qd == 0.0d0 .and. sd%rd == 0.0d0) then 
		sd%rd = Rpar(1)
		sd%qd = Rpar(1)
	end if
	
	if(associated(Hfunc)) then
		call qBasedUpdate(Hfunc, Rfunc, sd , sdNew, Rpar, Hpar, energy, Htype, Rtype)
	else
		call directBasedUpdate(damageFunc, derDamageFunc, Rfunc, sd , sdNew, Rpar, Hpar, energy, Htype, Rtype)
	end if
	
	if(Rtype==2) sdNew%derDamage = sdNew%derDamage/sdNew%rd

	call sdNew%putDamage(Param, nG)

!~ 	call writeDamage(sd%damage,sdNew%damage,energy)
end subroutine

!~ function damageImplex(damagePar) result(damageTilde)
!~ 	Real*8 , intent(in) :: damagePar(:)
!~ 	
!~ 	real*8 :: damageTilde, qTilde, rTilde, H, HPar(3), RPar(1)
!~ 	integer :: Rtype, Htype
!~ 
!~ 	call loadDamagePar(Hpar,Rpar,Htype,Rtype,damagePar)
!~ 		
!~ 	call chooseRfunc(Rtype)
!~ 	call chooseHfunc(Htype)
!~ 
!~ 	rTilde = rd + rdDot
!~ 	call Hfunc(H, HPar, rd , qd)
!~ 	qTilde = qd + H*(rTilde - rd)
!~ 	
!~ 	if(rTilde < 0.00000001d0) then
!~ 		rTilde = Rpar(1)
!~ 		qTilde = Rpar(1)
!~ 	end if  
!~ 	
!~ 	damageTilde = 1.0d0 - qTilde/rTilde
!~ 		
!~ end function

subroutine damageUpdateEquilibrium(F,matPar,damagePar,Param,constLaw,nG)
	use constitutiveLib , only : strainEnergy

	real*8 , intent(in) :: F(:,:) , matPar(:), damagePar(:)
	real*8 , intent(inout) :: Param(*)
	integer , intent(in) :: constLaw , nG
	real*8 :: energy
	integer :: IdamageIntegration
	
	IdamageIntegration = nint(damagePar(1))
	
	energy = 0.0d0
	if(IdamageIntegration == 1 .or. IdamageIntegration == 5 ) then ! Coupled
		call strainEnergy(energy,F,matPar,constLaw)
		call damageUpdate(damagePar,Param,energy,nG)
	end if

end subroutine

subroutine damageModifyStress(Ppk,Param,nG)
			
	Real*8 , intent(inout):: Ppk(:,:)
	Real*8 , intent(in) :: Param(*)
	integer, intent(in) :: nG
	logical :: flagDamageNew
	type(damageState) :: sd
	
	flagDamageNew = .true.
	call sd%loadDamage(Param, nG,flagDamageNew)
		
	Ppk = (1.0d0-sd%damage)*Ppk
	
end subroutine

subroutine damageModifyTangent(D,Ppk,damagePar,Param,nG,IcalledFrom) !!! also modifies stress!!!!!
	
	Real*8 , intent(inout):: D(:,:,:,:),Ppk(:,:)
	Real*8 , intent(in) :: damagePar(:),Param(*)
	integer, intent(in) :: nG, IcalledFrom
	logical :: flagDamageDer, flagDamageNew
	integer :: i,j,k,l, IdamageIntegration, NdimE
	type(damageState) :: sd
	
	NdimE = size(Ppk,1)
	
	!! IcalledFrom = (1 - equilibrium, 2- canonical Problem, 3 - Hom tang, 4 - PosProc)
	IdamageIntegration = nint(damagePar(1))
	if(IdamageIntegration == 0 .and. (IcalledFrom == 3 .or. IcalledFrom == 2)  ) IdamageIntegration = 1
	
	select case(IdamageIntegration)
		case(-1)
			return !!! Damage is not defined
		case(0,4)  !! 4 is for implex
			flagDamageDer = .false.
			flagDamageNew = .false.
		case(1)	!! implicit
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
	 
	call sd%loadDamage(Param, nG,flagDamageNew)
		
!~ 	if(IdamageIntegration == 4) damage = damageImplex(damagePar)
		
	D = (1.0d0-sd%damage)*D 
	
	if(sd%rdDot>0.0d0 .and. sd%derDamage>0.0d0 .and. flagDamageDer) then 
			
		do i = 1,NdimE
		do j = 1,NdimE
		do k = 1,NdimE
		do l = 1,NdimE
				D(i,j,k,l) = D(i,j,k,l) - sd%derDamage*Ppk(i,j)*Ppk(k,l)
		end do
		end do
		end do
		end do
	end if
	
	Ppk = (1.0d0-sd%damage)*Ppk
	
end subroutine


subroutine loadDamagePar(Hpar,Rpar,Htype,Rtype,damagePar)
	integer , intent(out) :: Htype, Rtype
	real*8  , intent(out) :: Rpar(1) , Hpar(3)
	real*8  , intent(in) :: damagePar(:)
	
	!! position is reserved to the integration scheme
	Rtype =  nint(damagePar(2))
	Rpar(1) = damagePar(3)
	Htype =  nint(damagePar(4))
	HPar(1:3) = damagePar(5 : 7)

end subroutine

subroutine loadDamage(s,Param, nG,inNew)
	class(damageState) :: s
	Real*8 , intent(in) :: Param(*)
	integer , intent(in) :: nG
	logical , intent(in) :: inNew
	integer :: iShiftDamageParamMod
	
	if(inNew) then
		iShiftDamageParamMod = LengthParam + (ng-1)*NdamageParam + iShiftDamageParam 
	else
		iShiftDamageParamMod = (ng-1)*NdamageParam + iShiftDamageParam
	end if
	
	s%damage = Param(iShiftDamageParamMod + Idamage)
	s%rd = Param(iShiftDamageParamMod + Ird)
	s%qd = Param(iShiftDamageParamMod + Iqd)
	s%rdDot = Param(iShiftDamageParamMod + IrdDot)
	s%derDamage = Param(iShiftDamageParamMod + IderDamage)
	
end subroutine

subroutine putDamage(s, Param, nG)
	class(damageState) :: s
	Real*8 , intent(out) :: Param(*)
	integer , intent(in) :: nG
	integer :: iShiftDamageParamMod

    iShiftDamageParamMod = LengthParam + (ng-1)*NdamageParam + iShiftDamageParam 
     
	Param(iShiftDamageParamMod + Idamage) = s%damage
	Param(iShiftDamageParamMod + Ird) = s%rd
	Param(iShiftDamageParamMod + Iqd) = s%qd 
	Param(iShiftDamageParamMod + IrdDot) = s%rdDot 
	Param(iShiftDamageParamMod + IderDamage) = s%derDamage
	
end subroutine


subroutine writeDamageFromParam(Param,energy,nG)
	real*8 , intent(in) :: Param(*), energy
	integer , intent(in) :: nG
	type(damageState) :: sd,sdNew
	
	call sd%loadDamage(Param, nG, .false.)
	call sdNew%loadDamage(Param, nG, .true.)
	
	call writeDamage(sd%damage, sdNew%damage, energy)
	
end subroutine

subroutine writeDamage(damage,damageNew, energy)
	real*8 , intent(in) :: damage, damageNew, energy
	integer , parameter :: OUnitDamage = 15, OUnitDamageNew = 16, OUnitStrainEnergy = 17
	
	open (OUnitDamage, file='damage.txt', Access = 'append')
	open (OUnitDamageNew, file='damageNew.txt', Access = 'append') 
	open (OUnitStrainEnergy, file='strainEnergy.txt', Access = 'append') 
	
	write(OUnitDamage,*) damage
	write(OUnitDamageNew,*) damageNew
	write(OUnitStrainEnergy,*) energy
	
	close(OUnitDamage)
	close(OUnitDamageNew)
	close(OUnitStrainEnergy)
end subroutine

end module
