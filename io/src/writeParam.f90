module writeParam
use converterData
use auxConverterFuncs
use coreConverterFuncs

implicit none

private constantDamage , slopeDamage , constantDamagePlate, slopeDamagePlate
public write2Param

contains

subroutine write2Param()

	integer , parameter :: OUnitParam = 24 
	character(len=30) :: formatSetVar , formatParam
	character(len=40) :: arg
	real(8) , allocatable :: ParamElem(:)
	integer :: i , j , ip , ipp , jpp, NParam , NinputMax, NkindFunc, NodG
	integer , allocatable :: kindFunc(:) , shiftParam(:)
	real(8) , allocatable :: xG(:), input(:,:)
	real(8) :: valueParam

	if(iargc() > 1 ) then
		call getarg(2, arg)
		write(*,*) arg
		read (arg,*) NParam, NkindFunc, NinputMax, NodG, NvolGroups, NFaceGroups
		allocate(kindFunc(NkindFunc),shiftParam(NkindFunc),input(NinputMax,NkindFunc))
		
		do i = 1 , NkindFunc
			call getarg(2 + i, arg)
			read (arg,*) kindFunc(i), shiftParam(i), input(:,i)
		end do		
	end if
	
	formatSetVar = "(a,/,I5)"
	write(formatParam,"(a,I1,a)") "(", NParam, "(E24.16,2X))"

	call readSGPMesh()

	allocate(xG(Ndim))
	allocate(ParamElem(NParam*NElem))

	ParamElem = 0.0
	
	do i = 1, Nelem
		ip = (i-1)*NodEl
		call findBaricentric(xG,X,Elem(ip+1:ip+NodG),Ndim,NodG)
		ipp = (i-1)*NParam
			
		do j = 1 , NkindFunc
				
			select case(kindFunc(j))
				case(1)
					call constantDamage(valueParam,xG,input(:,j))
				case(2)
					call slopeDamage(valueParam,xG,input(:,j))
				case(3)
					call constantDamagePlate(valueParam,xG,input(:,j))
				case(4)
					call slopeDamagePlate(valueParam,xG,input(:,j))
				case(5)
					call slopeDamageGen(valueParam,xG,input(:,j))
				case(6)
					call planeDamage(valueParam,xG,input(:,j))
				case(7)
					call setRotateThetaFrame(valueParam,xG,input(:,j),ElemGroup(i))
				case(8)
					call setRVElabel(valueParam,xG,input(:,j),ElemGroup(i))
				case default 
					valueParam = 0.0
			end select
				
			jpp = ipp + shiftParam(j) + 1
			ParamElem(jpp) = valueParam
			
		end do
		
	end do 

	open (OUnitParam, file='Param000.txt')
	write(OUnitParam,formatSetVar) "*Parameter Groups" , NElem + 1  ! it has always a neutral in order to link to others elements
	write(OUnitParam,"(a)") "*Real Parameters"
	write(OUnitParam,"(I3)") 0
	write(OUnitParam,"(I3)") ( NParam , i = 1 , Nelem ) 
	write(OUnitParam, formatParam)  ParamElem
	write(OUnitParam,"(a)") "*Integer Parameters"
	write(OUnitParam,"(I6,'*',I1)") NElem + 1, 0
	close(OUnitParam)

end subroutine

subroutine constantDamage(damage,xG,inputDamage)
	real*8 , intent(in) ::  xG(:) , inputDamage(:)
	real*8 , intent(out) :: damage
	real*8 :: xa,xb,maxDamage

	xa = inputDamage(1)
	xb = inputDamage(2)
	maxDamage = inputDamage(3)
	
	if(xG(1)<xa) then
		damage = 0.0
	else if(xG(1)>xb) then 
		damage = 0.0
	else
		damage = maxDamage
	end if

end subroutine

subroutine slopeDamage(damage,xG,inputDamage)
	real*8 , intent(in) ::  xG(:) , inputDamage(:)
	real*8 , intent(out) :: damage
	real*8 :: xa,xb,xc,xd,maxDamage
	
	xa = inputDamage(1)
	xb = inputDamage(2)
	xc = inputDamage(3)
	xd = inputDamage(4)
	maxDamage = inputDamage(5)
	
	if(xG(1)<=xa .or. xG(1)>=xd) then
		damage = 0.0
	else if(xG(1)>=xb .and. xG(1)<=xc) then
		damage = maxDamage
	else if(xG(1)<xb) then
		damage = maxDamage*(xG(1)-xa)/(xb-xa)
	else ! if(xG(1)>xc) then
		damage = maxDamage*(xG(1)-xd)/(xc-xd)
	end if

end subroutine

subroutine constantDamagePlate(damage,xG,inputDamage)
	real*8 , intent(in) ::  xG(:) , inputDamage(:)
	real*8 , intent(out) :: damage
	real*8 :: xa,xb,ya,yb,maxDamage
	
	xa = inputDamage(1)
	xb = inputDamage(2)
	ya = inputDamage(3)
	yb = inputDamage(4)
	maxDamage = inputDamage(5)
	
!~ 	write(6,*) xG 
	if(xG(1)<xa) then
		damage = 0.0
	else if(xG(1)>xb) then 
		damage = 0.0
	else if(xG(2)>yb) then 
		damage = 0.0
	else if(xG(2)<ya) then 
		damage = 0.0
	else
		damage = maxDamage
	end if

end subroutine

subroutine slopeDamagePlate(damage,xG,inputDamage)
	real*8 , intent(in) ::  xG(:) , inputDamage(:)
	real*8 , intent(out) :: damage
	real*8 :: x0,y0,ra,rb,rc,rd,maxDamage,r
	
	x0 = inputDamage(1)
	y0 = inputDamage(2)
	ra = inputDamage(3)
	rb = inputDamage(4)
	rc = inputDamage(5)
	rd = inputDamage(6)
	maxDamage = inputDamage(7)
	
	r = dsqrt( (xG(1)-x0)**2.0d0 + (xG(2)-y0)**2.0d0 )
	
	if(r<=ra .or. r>=rd) then
		damage = 0.0
	else if(r>=rb .and. r<=rc) then
		damage = maxDamage
	else if(r<rb) then
		damage = maxDamage*(r-ra)/(rb-ra)
	else ! if(r>rc) then
		damage = maxDamage*(r-rd)/(rc-rd)
	end if

end subroutine

subroutine slopeDamageGen(damage,xG,inputDamage)
	real*8 , intent(in) ::  xG(:) , inputDamage(:)
	real*8 , intent(out) :: damage
	real*8 :: x0,y0,z0,ra,rb,rc,rd,maxDamage,r
	
	x0 = inputDamage(1)
	y0 = inputDamage(2)
	z0 = inputDamage(3)
	ra = inputDamage(4)
	rb = inputDamage(5)
	rc = inputDamage(6)
	rd = inputDamage(7)
	maxDamage = inputDamage(8)
	
	r = dsqrt( (xG(1)-x0)**2.0d0 + (xG(2)-y0)**2.0d0 + (xG(3)-z0)**2.0d0  )
	
	if(r<=ra .or. r>=rd) then
		damage = 0.0
	else if(r>=rb .and. r<=rc) then
		damage = maxDamage
	else if(r<rb) then
		damage = maxDamage*(r-ra)/(rb-ra)
	else ! if(r>rc) then
		damage = maxDamage*(r-rd)/(rc-rd)
	end if
	
!~ 	write(0,*) damage

end subroutine

subroutine planeDamage(damage,xG,inputDamage)
	real*8 , intent(in) ::  xG(:) , inputDamage(:)
	real*8 , intent(out) :: damage
	integer , parameter :: ndim = 3
	real*8 :: x0(ndim),normal(ndim), thickness, maxDamage , minDamage, dist
	
	x0 = inputDamage(1:3)
	normal = inputDamage(4:6)
	normal = normal/dsqrt( dot_product(normal,normal))
	thickness = inputDamage(7)
	minDamage = inputDamage(8)
	maxDamage = inputDamage(9)

	dist = dot_product(xG - x0, normal )
	dist = dabs(dist)
	
	if( dist < thickness) then
		damage = (maxDamage*(thickness - dist) + dist*minDamage)/thickness ! sign(dist) * maxdamage
	else
		damage = 0.0
	end if
	
end subroutine

subroutine setRotateThetaFrame(ThetaFrame,xG,input,ElemGroup)
	real*8 , intent(in) ::  xG(:) , input(:)
	integer , intent(in) :: ElemGroup 
	real*8 , intent(out) :: ThetaFrame

	real*8 :: a,b, exc, x, y, xref, yref
	real*8 , parameter :: tol = 1.0d-6
	
	integer :: groupTarget
	
	xref = input(1)
	yref = input(2)
	exc = input(3) ! excentricity ellipsis := b/a (b in vertical axis, a in horizontal, thetaFrame is measured from the horizontal)
	groupTarget = nint(input(4)) 
	
	x = xG(1) - xref 
	y = xG(2) - yref 
	
	if(ElemGroup == groupTarget) then
		if(dabs(y)>tol) then
			ThetaFrame = datan(-x*exc*exc/y)
		else if(x>0.0d0) then
			ThetaFrame = -PI/2.0d0
		else
			ThetaFrame = PI/2.0d0
		end if
	else
		ThetaFrame = 0.0d0
	end if
	
end subroutine

subroutine setRVElabel(RVElabel,xG,input,ElemGroup)
	real*8 , intent(in) ::  xG(:) , input(:)
	integer , intent(in) :: ElemGroup 
	real*8 , intent(out) :: RVElabel

	integer, save :: RVEcount = 0 
	integer :: groupTarget
	
	
	groupTarget = nint(input(1)) 
	
	if(ElemGroup == groupTarget) then
		RVEcount = RVEcount + 1
		RVElabel = real(RVEcount)
	end if
		
		
end subroutine


end module
