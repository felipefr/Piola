module damageFibLib
use funcAux
use damageNewLib

implicit none

contains

subroutine damageUpdateFib(damagePar,Param,lambda)
	real*8, intent(in) :: damagePar(:), lambda 
	real*8, intent(inout) :: Param(*)
	real*8 :: lambda_c, lambda_f, dmax
	type(damageState) :: sd, sdNew 

	lambda_c = damagePar(1)
	lambda_f = damagePar(2)
	dmax = damagePar(3)
	
	call sd%loadDamage(Param, 1, .false.)
	
!~ 	if(lambda > lambda_c) then
!~ 		sdNew%damage = 0.999999999999d0
!~ 	else
!~ 		sdNew%damage = sd%damage
!~ 	end if 

	if(lambda > lambda_f) then
		sdNew%damage = dmax
	else if(lambda > lambda_c) then
		sdNew%damage = (lambda - lambda_c)*dmax/(lambda_f - lambda_c)
	else
		sdNew%damage = sd%damage
	end if 

	sdNew%rd = 0.0d0
	sdNew%qd = 0.0d0
	sdNew%rdDot = 0.0d0
	sdNew%derDamage = 0.0d0
	
!~ 	write(0,*) lambda
	
	call sdNew%putDamage(Param, 1)

end subroutine

subroutine damageModifySfib(Sfib,Param)
			
	Real*8 , intent(inout):: Sfib(:)
	Real*8 , intent(in) :: Param(*)
	logical :: flagDamageNew
	type(damageState) :: sd
	
	flagDamageNew = .true.
	 
	call sd%loadDamage(Param, 1,flagDamageNew)
				
	Sfib = (1.0d0-sd%damage)*Sfib 

	
end subroutine

subroutine damageModifyDSfib(DSfib,Param)
			
	Real*8 , intent(inout):: DSfib(:,:)
	Real*8 , intent(in) :: Param(*)
	logical :: flagDamageNew
	type(damageState) :: sd
	
	flagDamageNew = .true.
	 
	call sd%loadDamage(Param, 1,flagDamageNew)
				
	DSfib = (1.0d0-sd%damage)*DSfib 
	
end subroutine

end module
