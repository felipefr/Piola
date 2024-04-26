module utils
implicit none
private IprintMat, IprintVec , printVec, printMat
public assert , stop_error , FindKeyword , myPrint , countKeyword, isVinW, crossProduct, reallocate, calcGradU, PI

interface myPrint
	module procedure :: IprintMat
	module procedure :: IprintVec
	module procedure :: printMat
	module procedure :: printVec
end interface			

real*8, parameter :: PI=4.D0*DATAN(1.D0)

contains

character(len=20) function str(k)
!   "Convert an integer to string."
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
end function

subroutine calcGradU(GradU,Sol,dPhi_G,NdimE,NodElT,iDofSol)
	integer, intent(in) :: NdimE,NodElT , iDofSol !!! Not exactly iDofT
	Real*8, intent(out) :: GradU(NdimE,NdimE)
	Real*8, intent(in) ::dPhi_G(NdimE,NodElT) , Sol(iDofSol*NodElt)
	integer :: i, j , e ,ep 
	
	GradU = 0.0d0
	do e = 1 , NodElT
		ep = (e-1)*iDofSol
		do i = 1, NdimE
			do j = 1, NdimE
				GradU(i,j) = GradU(i,j) + Sol(ep+i)*dPhi_G(j,e)
			end do
		end do
	end do	
end subroutine

subroutine assert(condition) ! If condition == .false., it aborts the program.
	logical, intent(in) :: condition

	if (.not. condition) call stop_error("Assert failed.")
end subroutine

subroutine stop_error(msg) !Aborts the program with nonzero exit code
	character(len=*) :: msg ! Message to print on stdout
	print *, msg
	stop 1
end subroutine

subroutine FindKeyword(ioUnit,str1)
! Str1: String to find, ioUnit: Unit assigned to Input File
	Character(*) str1,str2*120
	integer :: n , iError, ioUnit
	
!~ 	rewind(IoUnit)   
	iError=0
	n = Len_Trim(str1)
	
	!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	do while (.true.)
		read(ioUnit, "(1A120)",iostat=iError) str2
		if (iError.lt.0 ) then 
			write(*,*) str1
			call stop_error("Keyword not finded")
		end if 
		if ( str2(1:n) == str1(1:n) )  exit 
	end do
	
end subroutine

subroutine CountKeyword(ocurrences,ioUnit,str1)
! Str1: String to find, ioUnit: Unit assigned to Input File
	Character(*) str1,str2*120
	integer :: n , iError, ioUnit, ocurrences
	
	rewind(IoUnit)   
	iError=0
	n = Len_Trim(str1)
	
	ocurrences = 0 
	!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	do while (.true.)
		read(ioUnit, "(1A120)",iostat=iError) str2
		if (iError.lt.0 ) exit 
		if ( str2(1:n) == str1(1:n) ) ocurrences = ocurrences + 1 
	end do
	
	rewind(IoUnit) 
	
end subroutine

subroutine printVec(v)
	real*8 :: v(:)
	integer n , i
	
	n = size(v)
	write(*,*) "=============="
	do i = 1, n
		write(*,"(g11.5)") v(i)
	end do
	write(*,*) "=============="
	
end subroutine

subroutine IprintVec(v)
	integer :: v(:)
	integer n , i
	
	n = size(v)
	write(*,*) "=============="
	do i = 1, n
		write(*,"(g11.5)") v(i)
	end do
	write(*,*) "=============="
	
end subroutine

subroutine printMat(A)
	real*8 :: A(:,:)
	integer n , m, i , j
	
	n = size(A,1)
	m = size(A,2)
	
	write(*,*) "=============================="
	do i = 1, n
		write(*,"(100g15.7,2X)") ( A(i,j) , j = 1 , m)
	end do
	write(*,*) "=============================="
	
end subroutine

subroutine IprintMat(A)
	integer :: A(:,:)
	integer n , m, i , j
	
	n = size(A,1)
	m = size(A,2)
	
	write(*,*) "=============================="
	do i = 1, n
		write(*,"(100g10.5)") ( A(i,j) , j = 1 , m)
	end do
	write(*,*) "=============================="
	
end subroutine

Logical function isVinW(v,w,n,m) ! W-V, m > n assumed
	integer, intent(in) :: v(n),w(m),n,m
	integer :: i,j,k
	
	isVinW = .False.
	outer: do i = 1,n
		do j = 1, m
			if(v(i)==w(j)) cycle outer 
		end do
		return 
	end do outer
	
	isVinW=.True. 
	 
end function

subroutine crossProduct(w,normW,u,v) ! w = (u x v)/|u x v|
	real*8 , intent(out) :: w(3), normW
	real*8 , intent(in) :: u(3),v(3)
	w(1) = u(2)*v(3) - v(2)*u(3)
	w(2) = u(3)*v(1) - v(3)*u(1)
	w(3) = u(1)*v(2) - v(1)*u(2)
	
	normW = dsqrt(dot_product(w,w)) 
	w = w/normW
	
end subroutine

subroutine reallocate (array, newsz) 
    real*8, pointer :: array(:) 
    integer newsz 
    real*8, pointer :: newarray(:) 
    allocate(newarray(newsz))
	newarray(1:ubound(array,1)) = array	
    array => newarray 
end subroutine 

end module
